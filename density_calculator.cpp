#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TRandom3.h"

// ============================================================================
// CONFIGURATION CONSTANTS
// ============================================================================
const double PT_MIN = 0.2; // min pT cut [GeV/c]
const int NY_BINS = 100;
const double Y_MIN = -5.0;
const double Y_MAX = 5.0;
const double BIN_WIDTH_Y = (Y_MAX - Y_MIN) / NY_BINS;
const int NPHI_BINS = 64;
const double PHI_MIN = -TMath::Pi();
const double PHI_MAX = TMath::Pi();
const double BIN_WIDTH_PHI = (PHI_MAX - PHI_MIN) / NPHI_BINS;
const int TENSOR_SAMPLES = 200; // increased to smooth spike at (0,0)
const double BIN_AREA = BIN_WIDTH_Y * BIN_WIDTH_PHI;

// ============================================================================
// UTILITY CLASSES (unchanged)
// ============================================================================
class StatusFilter {
private:
    std::set<int> allowed_status;
    bool accept_all;
public:
    StatusFilter() : accept_all(false) {
        allowed_status.insert(1);
        allowed_status.insert(11);
    }
    void setAllowedStatus(const std::vector<int>& status_codes) {
        allowed_status.clear();
        for (int status : status_codes) allowed_status.insert(status);
        accept_all = false;
    }
    void acceptAllStatus() { accept_all = true; allowed_status.clear(); }
    bool isAllowed(int status) const {
        if (accept_all) return true;
        return allowed_status.count(status);
    }
    void printConfiguration() const {
        if (accept_all) std::cout << "Status filter: ALL accepted" << std::endl;
        else {
            std::cout << "Status filter: [";
            for (int s : allowed_status) std::cout << s << " ";
            std::cout << "]" << std::endl;
        }
    }
};

class PDGManager {
public:
    static std::vector<int> expandPDGList(const std::vector<int>& input_pdgs) {
        std::set<int> pdg_set;
        for (int pdg : input_pdgs) { pdg_set.insert(pdg); pdg_set.insert(-pdg); }
        return std::vector<int>(pdg_set.begin(), pdg_set.end());
    }
    static std::string generateSuffix(int pdg1, int pdg2) {
        std::string s1 = (pdg1 < 0) ? "m" + std::to_string(-pdg1) : std::to_string(pdg1);
        std::string s2 = (pdg2 < 0) ? "m" + std::to_string(-pdg2) : std::to_string(pdg2);
        return s1 + "_" + s2;
    }
};

class KinematicsCalculator {
public:
    struct Particle {
        double px, py, pz, E;
        int pdg;
        int event_id;
        double getPt() const { return std::sqrt(px*px + py*py); }
        double getRapidity() const { return 0.5 * std::log((E + pz) / (E - pz + 1e-12)); }
        double getPhi() const { return std::atan2(py, px); }
    };
    static Particle createParticle(double px, double py, double pz, double E, int pdg, int eid) {
        return {px, py, pz, E, pdg, eid};
    }
    static double getDeltaY(const Particle& p1, const Particle& p2) { return p1.getRapidity() - p2.getRapidity(); }
    static double getDeltaPhi(const Particle& p1, const Particle& p2) {
        double dphi = p1.getPhi() - p2.getPhi();
        while (dphi >= TMath::Pi()) dphi -= 2*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
        return dphi;
    }
    static int getCharge(int pdg) {
        int apdg = std::abs(pdg);
        int s = (pdg > 0 ? 1 : -1);
        switch (apdg) {
            case 11: case 13: case 15: return -s;
            case 12: case 14: case 16: case 22: case 111: case 130: case 2112: case 310: case 311: case 3122: case 221: case 331: case 441: case 551: case 3322: return 0;
            case 211: case 321: case 411: case 421: case 431: case 511: case 521: case 531: case 2212: case 3222: return s;
            case 3112: case 3312: case 3334: return -s;
            default: return 0;
        }
    }
    static bool isChargedHadron(int pdg) {
        int charge = getCharge(pdg);
        int apdg = std::abs(pdg);
        return (charge != 0) && (apdg != 11 && apdg != 13 && apdg != 15);
    }
};

// ============================================================================
// HISTOGRAM MANAGER – only the 5 needed folders
// ============================================================================
class HistogramManager {
private:
    TFile* output_file;
    std::vector<TDirectory*> dir_singles, dir_pairs, dir_tensors, dir_c2s, dir_a2s;
    TDirectory *dir_event;
public:
    HistogramManager(TFile* file, int n_classes) : output_file(file) {
        dir_event = output_file->mkdir("EventPlots");
        for (int c = 0; c < n_classes; ++c) {
            std::string class_name = "Class" + std::to_string(c);
            TDirectory* class_dir = output_file->mkdir(class_name.c_str());
            dir_singles.push_back(class_dir->mkdir("SingleDensities"));
            dir_pairs.push_back(class_dir->mkdir("PairDensities"));
            dir_tensors.push_back(class_dir->mkdir("TensorProducts"));
            dir_c2s.push_back(class_dir->mkdir("C2"));
            dir_a2s.push_back(class_dir->mkdir("A2"));
        }
    }
    TH1F* createTH1F(const std::string& n, const std::string& t, int nb, double min, double max) {
        return new TH1F(n.c_str(), t.c_str(), nb, min, max);
    }
    TH2F* createTH2F(const std::string& n, const std::string& t, int nbx, double xmin, double xmax, int nby, double ymin, double ymax, const std::string& x_title = "", const std::string& y_title = "") {
        TH2F* h = new TH2F(n.c_str(), t.c_str(), nbx, xmin, xmax, nby, ymin, ymax);
        if (!x_title.empty()) h->GetXaxis()->SetTitle(x_title.c_str());
        if (!y_title.empty()) h->GetYaxis()->SetTitle(y_title.c_str());
        return h;
    }
    void writeSingle(int class_id, TObject* h) { if (class_id >= 0 && class_id < (int)dir_singles.size()) { dir_singles[class_id]->cd(); h->Write(); } }
    void writePair(int class_id, TObject* h) { if (class_id >= 0 && class_id < (int)dir_pairs.size()) { dir_pairs[class_id]->cd(); h->Write(); } }
    void writeTensor(int class_id, TObject* h) { if (class_id >= 0 && class_id < (int)dir_tensors.size()) { dir_tensors[class_id]->cd(); h->Write(); } }
    void writeC2(int class_id, TObject* h) { if (class_id >= 0 && class_id < (int)dir_c2s.size()) { dir_c2s[class_id]->cd(); h->Write(); } }
    void writeA2(int class_id, TObject* h) { if (class_id >= 0 && class_id < (int)dir_a2s.size()) { dir_a2s[class_id]->cd(); h->Write(); } }
    void writeEvent(TObject* h) { dir_event->cd(); h->Write(); }
    void writeEventSummary(int total_events) {
        dir_event->cd();
        TH1F* h = new TH1F("h_total_events", "Total Number of Events", 1, 0, 1);
        h->SetBinContent(1, total_events);
        h->Write();
        delete h;
    }
};

// ============================================================================
// SINGLE PARTICLE ANALYZER (unchanged)
// ============================================================================
class SingleParticleAnalyzer {
private:
    int target_pdg;
    std::vector<TH2F*> h_density_raw;
    std::vector<TH1F*> h_multiplicity;
    std::vector<int> total_events_per_class;
    int n_classes;
public:
    SingleParticleAnalyzer(int pdg, HistogramManager* hm, int n_classes) : target_pdg(pdg), n_classes(n_classes) {
        std::string suffix = (pdg < 0) ? "m" + std::to_string(-pdg) : std::to_string(pdg);
        for (int c = 0; c < n_classes; ++c) {
            std::string name = "h_density_raw_" + suffix + "_class" + std::to_string(c);
            h_density_raw.push_back(hm->createTH2F(name, "Raw Density " + suffix + " Class " + std::to_string(c),
                                                   NY_BINS, Y_MIN, Y_MAX, NPHI_BINS, PHI_MIN, PHI_MAX, "y", "#phi"));
            std::string mname = "h_mult_raw_" + suffix + "_class" + std::to_string(c);
            h_multiplicity.push_back(hm->createTH1F(mname, "Mult " + suffix + " Class " + std::to_string(c), 1000, 0, 1000));
            total_events_per_class.push_back(0);
        }
    }
    ~SingleParticleAnalyzer() {
        for (auto h : h_density_raw) if (h) delete h;
        for (auto h : h_multiplicity) if (h) delete h;
    }
    void processParticle(const KinematicsCalculator::Particle& p, int class_id) {
        if (p.pdg != target_pdg) return;
        if (p.getPt() < PT_MIN) return;
        double y = p.getRapidity();
        if (y < Y_MIN || y >= Y_MAX) return;
        h_density_raw[class_id]->Fill(y, p.getPhi(), 1.0);
    }
    void processEvent(int count, int class_id) {
        h_multiplicity[class_id]->Fill(count);
        total_events_per_class[class_id]++;
    }
    int getPDG() const { return target_pdg; }
    int getTotalEvents(int class_id) const { return total_events_per_class[class_id]; }
    double getMeanMultiplicity(int class_id) const {
        return (total_events_per_class[class_id] > 0) ? h_multiplicity[class_id]->GetMean() : 0.0;
    }
    TH2F* getDensity(int class_id) const {
        return (class_id >= 0 && class_id < (int)h_density_raw.size()) ? h_density_raw[class_id] : nullptr;
    }
    void finalize(HistogramManager* hm) {
        for (int c = 0; c < n_classes; ++c) {
            if (total_events_per_class[c] == 0) continue;
            std::string suffix = (target_pdg < 0) ? "m" + std::to_string(-target_pdg) : std::to_string(target_pdg);
            std::string cs = "_class" + std::to_string(c);
            std::string rho_name = "h_rho1_" + suffix + cs;
            TH2F* hf = (TH2F*) h_density_raw[c]->Clone(rho_name.c_str());
            hf->Scale(1.0 / total_events_per_class[c]);
            hf->Scale(1.0 / BIN_AREA); // true density
            hm->writeSingle(c, hf);
            delete hf;
        }
    }
};

// ============================================================================
// PAIR ANALYZER – fixed normalization & C2 writing
// ============================================================================
class PairAnalyzer {
private:
    int pdg1, pdg2;
    std::vector<TH2F*> h_pair_density_raw, h_tensor_product, h_c2_cumulant;
    std::vector<TH1F*> h_pair_count;
    std::vector<int> total_events_per_class;
    SingleParticleAnalyzer *s1, *s2;
    bool is_self;
    int n_classes;
public:
    PairAnalyzer(int p1, int p2, HistogramManager* hm, SingleParticleAnalyzer* sa1, SingleParticleAnalyzer* sa2, int n_classes)
        : pdg1(p1), pdg2(p2), s1(sa1), s2(sa2), n_classes(n_classes), is_self(p1 == p2) {
        std::string suffix = PDGManager::generateSuffix(p1, p2);
        for (int c = 0; c < n_classes; ++c) {
            std::string cs = "_class" + std::to_string(c);
            h_pair_density_raw.push_back(hm->createTH2F("h_rho2_raw_" + suffix + cs, "Raw #rho2 " + suffix + cs,
                                                        NY_BINS, Y_MIN, Y_MAX, NPHI_BINS, PHI_MIN, PHI_MAX, "#Deltay", "#Delta#phi"));
            h_tensor_product.push_back(hm->createTH2F("h_tensor_raw_" + suffix + cs, "Tensor " + suffix + cs,
                                                      NY_BINS, Y_MIN, Y_MAX, NPHI_BINS, PHI_MIN, PHI_MAX, "#Deltay", "#Delta#phi"));
            h_pair_count.push_back(hm->createTH1F("h_pair_mult_" + suffix + cs, "Pair mult " + suffix + cs, 1000, 0, 1000));
            h_c2_cumulant.push_back(nullptr);
            total_events_per_class.push_back(0);
        }
    }
    ~PairAnalyzer() {
        for (auto h : h_pair_density_raw) if (h) delete h;
        for (auto h : h_tensor_product) if (h) delete h;
        for (auto h : h_pair_count) if (h) delete h;
        for (auto h : h_c2_cumulant) if (h) delete h;
    }
    void processPairs(const std::vector<KinematicsCalculator::Particle>& parts, int class_id) {
        int count = 0;
        for (size_t i = 0; i < parts.size(); ++i) {
            if (parts[i].pdg != pdg1 || parts[i].getPt() < PT_MIN) continue;
            for (size_t j = 0; j < parts.size(); ++j) {
                if (is_self && i == j) continue;
                if (!is_self && (parts[j].pdg != pdg2 || parts[j].getPt() < PT_MIN)) continue;
                double dy = KinematicsCalculator::getDeltaY(parts[i], parts[j]);
                double dphi = KinematicsCalculator::getDeltaPhi(parts[i], parts[j]);
                h_pair_density_raw[class_id]->Fill(dy, dphi, 1.0);
                count++;
            }
        }
        h_pair_count[class_id]->Fill(count);
        total_events_per_class[class_id]++;
    }
    void calculateTensor(int class_id) {
        if (!s1 || !s2 || total_events_per_class[class_id] == 0) return;
        TH2F* rho1 = s1->getDensity(class_id);
        TH2F* rho2 = s2->getDensity(class_id);
        if (!rho1 || !rho2) return;
        h_tensor_product[class_id]->Reset();
        TRandom3 rand(0);
        double w = 1.0 / TENSOR_SAMPLES;
        double dx1 = rho1->GetXaxis()->GetBinWidth(1);
        double dy1 = rho1->GetYaxis()->GetBinWidth(1);
        double dx2 = rho2->GetXaxis()->GetBinWidth(1);
        double dy2 = rho2->GetYaxis()->GetBinWidth(1);
        for (int bx1 = 1; bx1 <= rho1->GetNbinsX(); ++bx1) {
            double yc1 = rho1->GetXaxis()->GetBinCenter(bx1);
            for (int by1 = 1; by1 <= rho1->GetNbinsY(); ++by1) {
                double v1 = rho1->GetBinContent(bx1, by1);
                if (v1 <= 0) continue;
                double phic1 = rho1->GetYaxis()->GetBinCenter(by1);
                for (int bx2 = 1; bx2 <= rho2->GetNbinsX(); ++bx2) {
                    double yc2 = rho2->GetXaxis()->GetBinCenter(bx2);
                    for (int by2 = 1; by2 <= rho2->GetNbinsY(); ++by2) {
                        double v2 = rho2->GetBinContent(bx2, by2);
                        if (v2 <= 0) continue;
                        double phic2 = rho2->GetYaxis()->GetBinCenter(by2);
                        for (int s = 0; s < TENSOR_SAMPLES; ++s) {
                            double y1 = yc1 + rand.Uniform(-dx1/2, dx1/2);
                            double phi1 = phic1 + rand.Uniform(-dy1/2, dy1/2);
                            double y2 = yc2 + rand.Uniform(-dx2/2, dx2/2);
                            double phi2 = phic2 + rand.Uniform(-dy2/2, dy2/2);
                            double dy = y1 - y2;
                            double dphi = phi1 - phi2;
                            while (dphi >= TMath::Pi()) dphi -= 2 * TMath::Pi();
                            while (dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
                            h_tensor_product[class_id]->Fill(dy, dphi, v1 * v2 * w);
                        }
                    }
                }
            }
        }
    }
    void computeC2(int class_id) {
        if (total_events_per_class[class_id] == 0) return;
        if (h_c2_cumulant[class_id]) delete h_c2_cumulant[class_id];
        std::string suffix = PDGManager::generateSuffix(pdg1, pdg2);
        std::string cs = "_class" + std::to_string(class_id);
        std::string c2_name = "h_c2_" + suffix + cs;
        TH2F* rho2_per_ev = (TH2F*) h_pair_density_raw[class_id]->Clone("temp_rho2");
        rho2_per_ev->Scale(1.0 / total_events_per_class[class_id]);
        TH2F* tensor_per_ev = (TH2F*) h_tensor_product[class_id]->Clone("temp_tensor");
        h_c2_cumulant[class_id] = (TH2F*) rho2_per_ev->Clone(c2_name.c_str());
        h_c2_cumulant[class_id]->Add(tensor_per_ev, -1.0);
        h_c2_cumulant[class_id]->Scale(1.0 / BIN_AREA);
        delete rho2_per_ev;
        delete tensor_per_ev;
    }
    TH2F* getC2(int class_id) { return (class_id < (int)h_c2_cumulant.size()) ? h_c2_cumulant[class_id] : nullptr; }
    TH2F* getTensor(int class_id) { return (class_id < (int)h_tensor_product.size()) ? h_tensor_product[class_id] : nullptr; }
    void finalize(HistogramManager* hm) {
        for (int c = 0; c < n_classes; ++c) {
            if (total_events_per_class[c] == 0) continue;
            calculateTensor(c);
            computeC2(c);

            // ─────────────────────────────────────────────────────────────
            // Write C2 histogram (this was missing before)
            // ─────────────────────────────────────────────────────────────
            if (h_c2_cumulant[c]) {
                std::string name = h_c2_cumulant[c]->GetName();
                TH2F* copy_for_write = (TH2F*) h_c2_cumulant[c]->Clone(name.c_str());
                hm->writeC2(c, copy_for_write);
                std::cout << "[C2 written] " << name
                          << "  entries = " << copy_for_write->GetEntries()
                          << "  mean = " << copy_for_write->GetMean() << std::endl;
                // delete copy_for_write;   // optional
            }

            std::string suffix = PDGManager::generateSuffix(pdg1, pdg2);
            std::string cs = "_class" + std::to_string(c);

            std::string rho2_name = "h_rho2_" + suffix + cs;
            TH2F* h2 = (TH2F*) h_pair_density_raw[c]->Clone(rho2_name.c_str());
            h2->Scale(1.0 / total_events_per_class[c] / BIN_AREA);
            hm->writePair(c, h2); delete h2;

            std::string tensor_name = "h_tensor_" + suffix + cs;
            TH2F* ht = (TH2F*) h_tensor_product[c]->Clone(tensor_name.c_str());
            ht->Scale(1.0 / total_events_per_class[c] / BIN_AREA);
            hm->writeTensor(c, ht); delete ht;
        }
    }
};

// ============================================================================
// MAIN PROCESSOR – only C2 and A2
// ============================================================================
class HepMCProcessor {
private:
    StatusFilter status_filter;
    HistogramManager* hist_manager;
    std::vector<SingleParticleAnalyzer*> singles;
    std::vector<PairAnalyzer*> pairs;
    std::map<std::pair<int,int>, PairAnalyzer*> pair_map;
    TH1F* h_total_mult;
    int events_processed = 0;
    int n_classes;
    std::vector<double> mult_thresholds;
public:
    HepMCProcessor(HistogramManager* manager, int n_classes) : hist_manager(manager), n_classes(n_classes) {
        h_total_mult = hist_manager->createTH1F("h_total_event_multiplicity_raw", "Charged |y|<1 Multiplicity", 500, 0, 500);
    }
    ~HepMCProcessor() {
        for (auto s : singles) delete s;
        for (auto p : pairs) delete p;
        delete h_total_mult;
    }
    void setStatusFilter(const std::vector<int>& s) { status_filter.setAllowedStatus(s); }
    void addSingle(int pdg) { singles.push_back(new SingleParticleAnalyzer(pdg, hist_manager, n_classes)); }
    void addPair(int trigger, int assoc) {
        SingleParticleAnalyzer *s_trig = nullptr, *s_assoc = nullptr;
        for (auto s : singles) {
            if (s->getPDG() == trigger) s_trig = s;
            if (s->getPDG() == assoc) s_assoc = s;
        }
        if (s_trig && s_assoc) {
            std::pair<int,int> key{trigger, assoc};
            if (pair_map.find(key) == pair_map.end()) {
                auto* pa = new PairAnalyzer(trigger, assoc, hist_manager, s_trig, s_assoc, n_classes);
                pairs.push_back(pa);
                pair_map[key] = pa;
            }
        }
    }
    SingleParticleAnalyzer* findSingleAnalyzer(int pdg) {
        for (auto s : singles) if (s->getPDG() == pdg) return s;
        return nullptr;
    }
    PairAnalyzer* getPairAnalyzer(int trig, int assoc) {
        auto key = std::make_pair(trig, assoc);
        auto it = pair_map.find(key);
        return (it != pair_map.end()) ? it->second : nullptr;
    }
    void computeMultiplicityThresholds(const std::string& fname) {
        std::ifstream file(fname);
        if (!file.is_open()) return;
        std::string line;
        int eid = 0;
        std::vector<KinematicsCalculator::Particle> event_parts;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == 'E') {
                if (eid > 0) {
                    int charged_mult = 0;
                    for (const auto& p : event_parts) {
                        if (std::abs(p.getRapidity()) < 1.0 && KinematicsCalculator::isChargedHadron(p.pdg))
                            charged_mult++;
                    }
                    h_total_mult->Fill(charged_mult);
                    event_parts.clear();
                }
                eid++;
                continue;
            }
            if (line[0] != 'P') continue;
            std::stringstream ss(line);
            char c; int id, pdg, stat;
            double px, py, pz, E, m;
            ss >> c >> id >> pdg >> px >> py >> pz >> E >> m >> stat;
            if (status_filter.isAllowed(stat)) {
                event_parts.push_back(KinematicsCalculator::createParticle(px, py, pz, E, pdg, eid));
            }
        }
        if (!event_parts.empty()) {
            int charged_mult = 0;
            for (const auto& p : event_parts) {
                if (std::abs(p.getRapidity()) < 1.0 && KinematicsCalculator::isChargedHadron(p.pdg))
                    charged_mult++;
            }
            h_total_mult->Fill(charged_mult);
        }
        file.close();
        std::vector<double> mults;
        for (int b = 1; b <= h_total_mult->GetNbinsX(); ++b) {
            for (int i = 0; i < h_total_mult->GetBinContent(b); ++i) {
                mults.push_back(h_total_mult->GetBinCenter(b));
            }
        }
        std::sort(mults.begin(), mults.end());
        double total = mults.size();
        if (total == 0) return;
        mult_thresholds.resize(n_classes - 1);
        double frac = 1.0;
        for (int k = 0; k < n_classes - 1; ++k) {
            frac -= 1.0 / n_classes;
            size_t idx = static_cast<size_t>(total * frac);
            if (idx >= total) idx = total - 1;
            mult_thresholds[k] = mults[idx];
        }
    }
    int getClass(int mult) const {
        for (int c = 0; c < n_classes - 1; ++c) {
            if (mult > mult_thresholds[c]) return c;
        }
        return n_classes - 1;
    }
    void processFile(const std::string& fname) {
        computeMultiplicityThresholds(fname);
        std::ifstream file(fname);
        if (!file.is_open()) return;
        status_filter.printConfiguration();
        std::string line;
        int eid = 0;
        std::vector<KinematicsCalculator::Particle> event_parts;
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == 'E') {
                if (eid > 0) {
                    int charged_mult = 0;
                    for (const auto& p : event_parts) {
                        if (std::abs(p.getRapidity()) < 1.0 && KinematicsCalculator::isChargedHadron(p.pdg))
                            charged_mult++;
                    }
                    int class_id = getClass(charged_mult);
                    processEvent(event_parts, class_id);
                    event_parts.clear();
                }
                eid++;
                if (eid % 100 == 0) std::cout << "Processing event " << eid << std::endl;
                continue;
            }
            if (line[0] != 'P') continue;
            std::stringstream ss(line);
            char c; int id, pdg, stat;
            double px, py, pz, E, m;
            ss >> c >> id >> pdg >> px >> py >> pz >> E >> m >> stat;
            if (status_filter.isAllowed(stat)) {
                event_parts.push_back(KinematicsCalculator::createParticle(px, py, pz, E, pdg, eid));
            }
        }
        if (!event_parts.empty()) {
            int charged_mult = 0;
            for (const auto& p : event_parts) {
                if (std::abs(p.getRapidity()) < 1.0 && KinematicsCalculator::isChargedHadron(p.pdg))
                    charged_mult++;
            }
            int class_id = getClass(charged_mult);
            processEvent(event_parts, class_id);
        }
        events_processed = eid;
        file.close();
        std::cout << "Processed " << events_processed << " events.\n";
    }
    void processEvent(const std::vector<KinematicsCalculator::Particle>& parts, int class_id) {
        int charged_mult = 0;
        for (const auto& p : parts) {
            if (std::abs(p.getRapidity()) < 1.0 && KinematicsCalculator::isChargedHadron(p.pdg))
                charged_mult++;
        }
        h_total_mult->Fill(charged_mult);
        for (auto s : singles) {
            int count = 0;
            for (const auto& p : parts) {
                s->processParticle(p, class_id);
                if (p.pdg == s->getPDG() && p.getPt() >= PT_MIN) count++;
            }
            s->processEvent(count, class_id);
        }
        for (auto p : pairs) p->processPairs(parts, class_id);
    }
    void computeDerivedPhysics() {
        if (events_processed == 0) return;
        std::cout << "\n--- Computing C2 and A2 ---\n";
        for (auto& kv : pair_map) {
            int trig = kv.first.first;
            int assoc = kv.first.second;
            PairAnalyzer* pa = kv.second;
            for (int c = 0; c < n_classes; ++c) {
                TH2F* c2 = pa->getC2(c);
                if (!c2) continue;
                SingleParticleAnalyzer* s_assoc = nullptr;
                for (auto* s : singles) {
                    if (s->getPDG() == assoc) {
                        s_assoc = s;
                        break;
                    }
                }
                if (!s_assoc) continue;
                double N_beta = s_assoc->getMeanMultiplicity(c); // N_beta = associate
                if (N_beta <= 0) continue;
                std::string suf = PDGManager::generateSuffix(trig, assoc);
                std::string cs = "_class" + std::to_string(c);
                std::string a2_name = "h_A2_" + suf + cs;
                TH2F* h_a2 = (TH2F*) c2->Clone(a2_name.c_str());
                h_a2->Scale(1.0 / N_beta);
                hist_manager->writeA2(c, h_a2);

                std::cout << "[A2 written] " << a2_name
                          << "  N_beta = " << N_beta
                          << "  mean(C2 input) = " << c2->GetMean()
                          << "  mean(A2) = " << h_a2->GetMean() << std::endl;

                delete h_a2;
            }
        }
    }
    void finalize() {
        for (auto s : singles) s->finalize(hist_manager);
        for (auto p : pairs) p->finalize(hist_manager);
        computeDerivedPhysics();
        TH1F* hfm = (TH1F*) h_total_mult->Clone("h_total_event_multiplicity_norm");
        if (events_processed > 0) hfm->Scale(1.0 / events_processed);
        hist_manager->writeEvent(hfm);
        delete hfm;
        hist_manager->writeEventSummary(events_processed);
    }
};

// ============================================================================
// MAIN
// ============================================================================
int main(int argc, char* argv[]) {
    TH1::AddDirectory(false);
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " input.hepmc output.root n_classes PDG1 [PDG2 ...]\n";
        std::cerr << "Example: " << argv[0] << " events.hepmc output.root 3 211 321 2212\n";
        return 1;
    }
    std::string fin = argv[1];
    std::string fout = argv[2];
    int n_classes = std::stoi(argv[3]);
    std::vector<int> inputs;
    for (int i = 4; i < argc; ++i) {
        try { inputs.push_back(std::stoi(argv[i])); } catch (...) {}
    }
    if (inputs.empty()) {
        std::cerr << "No PDG codes provided.\n";
        return 1;
    }
    auto pdgs = PDGManager::expandPDGList(inputs);
    TFile* file = new TFile(fout.c_str(), "RECREATE");
    HistogramManager hm(file, n_classes);
    HepMCProcessor proc(&hm, n_classes);
    proc.setStatusFilter({1, 11});
    for (int p : pdgs) proc.addSingle(p);
    for (size_t i = 0; i < pdgs.size(); ++i) {
        for (size_t j = 0; j < pdgs.size(); ++j) {
            proc.addPair(pdgs[i], pdgs[j]);
        }
    }
    std::cout << "Setup done. Processing file: " << fin << "\n";
    proc.processFile(fin);
    proc.finalize();
    file->Close();
    delete file;
    std::cout << "Done.\n";
    return 0;
}
