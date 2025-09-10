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

class StatusFilter {
private:
    std::set<int> allowed_status;
    bool accept_all;
    
public:
    StatusFilter() : accept_all(false) {
        allowed_status.insert(1);
    }
    
    void setAllowedStatus(const std::vector<int>& status_codes) {
        allowed_status.clear();
        for (int status : status_codes) {
            allowed_status.insert(status);
        }
        accept_all = false;
    }
    
    void acceptAllStatus() {
        accept_all = true;
        allowed_status.clear();
    }
    
    bool isAllowed(int status) const {
        if (accept_all) return true;
        return allowed_status.find(status) != allowed_status.end();
    }
    
    void printConfiguration() const {
        if (accept_all) {
            std::cout << "Status filter: ALL particle status codes accepted" << std::endl;
        } else {
            std::cout << "Status filter: Accepting particles with status [";
            bool first = true;
            for (int status : allowed_status) {
                if (!first) std::cout << ",";
                std::cout << status;
                first = false;
            }
            std::cout << "]" << std::endl;
        }
    }
    
    std::vector<int> getStatusCodes() const {
        std::vector<int> codes;
        if (!accept_all) {
            for (int status : allowed_status) {
                codes.push_back(status);
            }
        }
        return codes;
    }
    
    bool isDefaultStatus() const {
        return !accept_all && allowed_status.size() == 1 && allowed_status.count(1) == 1;
    }
};

class PDGNamer {
public:
    static std::string generateSuffix(int pdg_code, const std::vector<int>& status_codes = {}, bool is_default_status = true) {
        std::string suffix = "";
        
        // Handle negative PDG codes with 'm' prefix (m = minus)
        // This avoids issues with minus signs in ROOT histogram names
        // Example: PDG -211 becomes "m211" in histogram names
        if (pdg_code < 0) {
            suffix += "m" + std::to_string(-pdg_code);
        } else {
            suffix += std::to_string(pdg_code);
        }
        
        // Only add status codes if not default (status=1 only)
        if (!is_default_status && !status_codes.empty()) {
            suffix += "_status";
            for (size_t i = 0; i < status_codes.size(); i++) {
                suffix += std::to_string(status_codes[i]);
                if (i < status_codes.size() - 1) {
                    suffix += "_";
                }
            }
        }
        
        return suffix;
    }
    
    static std::string generatePairSuffix(int pdg1, int pdg2, const std::vector<int>& status_codes = {}, bool is_default_status = true) {
        std::string suffix1 = (pdg1 < 0) ? "m" + std::to_string(-pdg1) : std::to_string(pdg1);
        std::string suffix2 = (pdg2 < 0) ? "m" + std::to_string(-pdg2) : std::to_string(pdg2);
        std::string suffix = suffix1 + "_" + suffix2;
        
        if (!is_default_status && !status_codes.empty()) {
            suffix += "_status";
            for (size_t i = 0; i < status_codes.size(); i++) {
                suffix += std::to_string(status_codes[i]);
                if (i < status_codes.size() - 1) {
                    suffix += "_";
                }
            }
        }
        
        return suffix;
    }
    
    static std::string generateTitle(int pdg_code, const std::vector<int>& status_codes = {}, bool is_default_status = true) {
        std::string title = "PDG[" + std::to_string(pdg_code) + "]";
        
        if (!is_default_status && !status_codes.empty()) {
            title += " Status[";
            for (size_t i = 0; i < status_codes.size(); i++) {
                title += std::to_string(status_codes[i]);
                if (i < status_codes.size() - 1) {
                    title += ",";
                }
            }
            title += "]";
        }
        
        return title;
    }
    
    static std::string generatePairTitle(int pdg1, int pdg2, const std::vector<int>& status_codes = {}, bool is_default_status = true) {
        std::string title = "PDG[" + std::to_string(pdg1) + "," + std::to_string(pdg2) + "]";
        
        if (!is_default_status && !status_codes.empty()) {
            title += " Status[";
            for (size_t i = 0; i < status_codes.size(); i++) {
                title += std::to_string(status_codes[i]);
                if (i < status_codes.size() - 1) {
                    title += ",";
                }
            }
            title += "]";
        }
        
        return title;
    }
    
    static void printNamingConvention() {
        std::cout << "\nHistogram Naming Convention:" << std::endl;
        std::cout << "- Positive PDG codes: 211 → '211'" << std::endl;
        std::cout << "- Negative PDG codes: -211 → 'm211' (m = minus)" << std::endl;
        std::cout << "- Default status (1): no status suffix" << std::endl;
        std::cout << "- Custom status: adds '_status[codes]' suffix" << std::endl;
    }
};

class HistogramManager {
private:
    TFile* output_file;
    
public:
    HistogramManager(TFile* file) : output_file(file) {}
    
    TH1F* createTH1F(const std::string& name, const std::string& title,
                     int nbins, double xmin, double xmax) {
        return new TH1F(name.c_str(), title.c_str(), nbins, xmin, xmax);
    }
    
    TH2F* createTH2F(const std::string& name, const std::string& title,
                     int nbinsx, double xmin, double xmax,
                     int nbinsy, double ymin, double ymax) {
        return new TH2F(name.c_str(), title.c_str(),
                       nbinsx, xmin, xmax, nbinsy, ymin, ymax);
    }
    
    void writeHistogram(TH1* hist) {
        output_file->cd();
        hist->Write();
    }
};

class KinematicsCalculator {
public:
    struct Particle {
        double px, py, pz, E;
        int pdg;
        int event_id;
        
        double getPt() const {
            return std::sqrt(px*px + py*py);
        }
        
        double getRapidity() const {
            return 0.5 * std::log((E + pz) / (E - pz + 1e-12));
        }
        
        double getPhi() const {
            return std::atan2(py, px);
        }
        
        double getEta() const {
            double p = std::sqrt(px*px + py*py + pz*pz);
            return 0.5 * std::log((p + pz) / (p - pz + 1e-12));
        }
    };
    
    static Particle createParticle(double px, double py, double pz, double E, int pdg, int event_id) {
        Particle p;
        p.px = px; p.py = py; p.pz = pz; p.E = E; p.pdg = pdg; p.event_id = event_id;
        return p;
    }
    
    static double getDeltaY(const Particle& p1, const Particle& p2) {
        return p1.getRapidity() - p2.getRapidity();
    }
    
    static double getDeltaPhi(const Particle& p1, const Particle& p2) {
        double dphi = p1.getPhi() - p2.getPhi();
        while (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();
        while (dphi < -TMath::Pi()) dphi += 2*TMath::Pi();
        return dphi;
    }
};



class SingleParticleAnalyzer {
private:
    int target_pdg;
    std::vector<int> status_codes;
    bool is_default_status;
    HistogramManager* hist_manager;
    std::string pdg_suffix;
    std::string pdg_title;
    int total_events;
    
    struct Config {
        int n_rapidity_bins = 100;
        double rapidity_min = -5.0;
        double rapidity_max = 5.0;
        int n_phi_bins = 64;
        double phi_min = -TMath::Pi();
        double phi_max = TMath::Pi();
        int n_pt_bins = 200;
        double pt_min = 0.0;
        double pt_max = 20.0;
    } config;
    
    TH2F* h_density_rapidity_phi;
    TH1F* h_pt;
    TH1F* h_rapidity;
    TH1F* h_phi;
    TH1F* h_multiplicity;
    
public:
    SingleParticleAnalyzer(int pdg_code, HistogramManager* manager, const std::vector<int>& status_list = {1}, bool default_status = true)
        : target_pdg(pdg_code),
          status_codes(status_list),
          is_default_status(default_status),
          hist_manager(manager),
          total_events(0) {
        
        pdg_suffix = PDGNamer::generateSuffix(pdg_code, status_codes, is_default_status);
        pdg_title = PDGNamer::generateTitle(pdg_code, status_codes, is_default_status);
        
        initializeHistograms();
    }
    
    void processParticle(const KinematicsCalculator::Particle& particle) {
        if (particle.pdg != target_pdg) return;
        
        double pt = particle.getPt();
        if (pt < 1e-6) return;
        
        double rapidity = particle.getRapidity();
        double phi = particle.getPhi();
        
        h_density_rapidity_phi->Fill(rapidity, phi, 1.0/pt);
        h_pt->Fill(pt, 1.0/pt);
        h_rapidity->Fill(rapidity);
        h_phi->Fill(phi);
    }
    
    void processEvent(int count) {
        h_multiplicity->Fill(count);
        total_events++;
    }
    
    TH2F* getDensityHistogram() const { return h_density_rapidity_phi; }
    
    int getPDG() const { return target_pdg; }
    
    std::vector<int> getStatusCodes() const { return status_codes; }
    
    int getTotalEvents() const { return total_events; }
    
    void finalize() {
        std::string final_name = "h_final_density_" + pdg_suffix;
        std::string final_title = "#rho_{1}(y,#phi)=dN/dp_{T}(y,#phi) " + pdg_title + ";rapidity y;#phi [rad]";
        
        TH2F* h_final_density = (TH2F*)h_density_rapidity_phi->Clone(final_name.c_str());
        h_final_density->SetTitle(final_title.c_str());
        
        // Normalize by number of events and bin size
        double dy = (config.rapidity_max - config.rapidity_min) / config.n_rapidity_bins;
        double dphi = (config.phi_max - config.phi_min) / config.n_phi_bins;
        if (total_events > 0) {
            h_final_density->Scale(1.0 / (total_events * dy * dphi));
        }
        
        hist_manager->writeHistogram(h_final_density);
        hist_manager->writeHistogram(h_density_rapidity_phi);
        hist_manager->writeHistogram(h_pt);
        hist_manager->writeHistogram(h_rapidity);
        hist_manager->writeHistogram(h_phi);
        hist_manager->writeHistogram(h_multiplicity);
        
        std::cout << "Single particle analysis for " << pdg_title
                  << ": " << h_rapidity->GetEntries() << " particles in " 
                  << total_events << " events" << std::endl;
        
        delete h_final_density;
    }
    
private:
    void initializeHistograms() {
        std::string base_name, base_title;
        
        base_name = "h_density_rapidity_phi_" + pdg_suffix;
        base_title = "#rho_{1}(y,#phi) raw " + pdg_title + ";rapidity y;#phi [rad]";
        h_density_rapidity_phi = hist_manager->createTH2F(base_name, base_title,
            config.n_rapidity_bins, config.rapidity_min, config.rapidity_max,
            config.n_phi_bins, config.phi_min, config.phi_max);
        
        base_name = "h_pt_" + pdg_suffix;
        base_title = "p_{T} distribution " + pdg_title + ";p_{T} [GeV];dN/dp_{T}";
        h_pt = hist_manager->createTH1F(base_name, base_title,
            config.n_pt_bins, config.pt_min, config.pt_max);
        
        base_name = "h_rapidity_" + pdg_suffix;
        base_title = "Rapidity distribution " + pdg_title + ";y;counts";
        h_rapidity = hist_manager->createTH1F(base_name, base_title,
            config.n_rapidity_bins, config.rapidity_min, config.rapidity_max);
        
        base_name = "h_phi_" + pdg_suffix;
        base_title = "#phi distribution " + pdg_title + ";#phi [rad];counts";
        h_phi = hist_manager->createTH1F(base_name, base_title,
            config.n_phi_bins, config.phi_min, config.phi_max);
        
        base_name = "h_multiplicity_" + pdg_suffix;
        base_title = "Event multiplicity " + pdg_title + ";N_{particles};events";
        h_multiplicity = hist_manager->createTH1F(base_name, base_title, 1000, 0, 1000);
    }
};

class PairAnalyzer {
private:
    int pdg1, pdg2;
    std::vector<int> status_codes;
    bool is_default_status;
    HistogramManager* hist_manager;
    std::string pair_suffix;
    std::string pair_title;
    int total_events;
    
    struct Config {
        int n_delta_y_bins = 100;
        double delta_y_min = -4.0;
        double delta_y_max = 4.0;
        int n_delta_phi_bins = 64;
        double delta_phi_min = -TMath::Pi();
        double delta_phi_max = TMath::Pi();
    } config;
    
    TH2F* h_pair_density;
    TH2F* h_tensor_product;
    TH1F* h_delta_y;
    TH1F* h_delta_phi;
    TH1F* h_pair_multiplicity;
    
    SingleParticleAnalyzer* analyzer1;
    SingleParticleAnalyzer* analyzer2;
    
public:
    PairAnalyzer(int particle1_pdg, int particle2_pdg, HistogramManager* manager,
                 SingleParticleAnalyzer* single1, SingleParticleAnalyzer* single2,
                 const std::vector<int>& status_list = {1}, bool default_status = true)
        : pdg1(particle1_pdg), pdg2(particle2_pdg), hist_manager(manager),
          analyzer1(single1), analyzer2(single2), status_codes(status_list), 
          is_default_status(default_status), total_events(0) {
        
        pair_suffix = PDGNamer::generatePairSuffix(pdg1, pdg2, status_codes, is_default_status);
        pair_title = PDGNamer::generatePairTitle(pdg1, pdg2, status_codes, is_default_status);
        
        initializeHistograms();
    }
    
    void processPairs(const std::vector<KinematicsCalculator::Particle>& particles) {
        int pair_count = 0;
        
        for (size_t i = 0; i < particles.size(); i++) {
            if (particles[i].pdg != pdg1) continue;
            
            for (size_t j = 0; j < particles.size(); j++) {
                if (i == j) continue;
                if (particles[j].pdg != pdg2) continue;
                
                if (particles[i].event_id != particles[j].event_id) {
                    std::cerr << "Warning: Particles from different events!" << std::endl;
                    continue;
                }
                
                const auto& p1 = particles[i];
                const auto& p2 = particles[j];
                
                double pt1 = p1.getPt();
                double pt2 = p2.getPt();
                if (pt1 < 1e-6 || pt2 < 1e-6) continue;
                
                double delta_y = KinematicsCalculator::getDeltaY(p1, p2);
                double delta_phi = KinematicsCalculator::getDeltaPhi(p1, p2);
                
                double pair_weight = 1.0 / (pt1 + pt2);
                
                h_pair_density->Fill(delta_y, delta_phi, pair_weight);
                h_delta_y->Fill(delta_y);
                h_delta_phi->Fill(delta_phi);
                
                pair_count++;
            }
        }
        
        h_pair_multiplicity->Fill(pair_count);
        total_events++;
    }
    
    void calculateTensorProduct() {
        if (!analyzer1 || !analyzer2) {
            std::cout << "Warning: Cannot calculate tensor product - missing analyzers" << std::endl;
            return;
        }
        
        TH2F* rho1 = analyzer1->getDensityHistogram();
        TH2F* rho2 = analyzer2->getDensityHistogram();
        
        if (!rho1 || !rho2) {
            std::cout << "Warning: Missing density histograms for tensor product" << std::endl;
            return;
        }
        
        // Get total events from single particle analyzers
        int events1 = analyzer1->getTotalEvents();
        int events2 = analyzer2->getTotalEvents();
        int min_events = std::min(events1, events2);
        
        if (min_events == 0) {
            std::cout << "Warning: No events processed for tensor product calculation" << std::endl;
            return;
        }
        
        // Create normalized copies for tensor product calculation
        TH2F* rho1_norm = (TH2F*)rho1->Clone("temp_rho1");
        TH2F* rho2_norm = (TH2F*)rho2->Clone("temp_rho2");
        
        // Apply same normalization as in single particle analysis
        double dy_single = 10.0 / 100;  // (5.0 - (-5.0)) / 100
        double dphi_single = (2*TMath::Pi()) / 64;
        
        rho1_norm->Scale(1.0 / (min_events * dy_single * dphi_single));
        rho2_norm->Scale(1.0 / (min_events * dy_single * dphi_single));
        
        // Calculate tensor product: rho1(y1,phi1) * rho2(y1+delta_y, phi1+delta_phi)
        for (int i = 1; i <= h_tensor_product->GetNbinsX(); i++) {
            for (int j = 1; j <= h_tensor_product->GetNbinsY(); j++) {
                double delta_y = h_tensor_product->GetXaxis()->GetBinCenter(i);
                double delta_phi = h_tensor_product->GetYaxis()->GetBinCenter(j);
                
                double tensor_sum = 0.0;
                int valid_bins = 0;
                
                // Sum over all possible (y1, phi1) positions
                for (int y1_bin = 1; y1_bin <= rho1_norm->GetNbinsX(); y1_bin++) {
                    for (int phi1_bin = 1; phi1_bin <= rho1_norm->GetNbinsY(); phi1_bin++) {
                        double y1 = rho1_norm->GetXaxis()->GetBinCenter(y1_bin);
                        double phi1 = rho1_norm->GetYaxis()->GetBinCenter(phi1_bin);
                        
                        double y2 = y1 + delta_y;
                        double phi2 = phi1 + delta_phi;
                        
                        // Handle phi periodicity
                        while (phi2 > TMath::Pi()) phi2 -= 2*TMath::Pi();
                        while (phi2 < -TMath::Pi()) phi2 += 2*TMath::Pi();
                        
                        // Find corresponding bin in rho2
                        int y2_bin = rho2_norm->GetXaxis()->FindBin(y2);
                        int phi2_bin = rho2_norm->GetYaxis()->FindBin(phi2);
                        
                        // Check if within bounds
                        if (y2_bin >= 1 && y2_bin <= rho2_norm->GetNbinsX() &&
                            phi2_bin >= 1 && phi2_bin <= rho2_norm->GetNbinsY()) {
                            
                            double rho1_val = rho1_norm->GetBinContent(y1_bin, phi1_bin);
                            double rho2_val = rho2_norm->GetBinContent(y2_bin, phi2_bin);
                            
                            tensor_sum += rho1_val * rho2_val;
                            valid_bins++;
                        }
                    }
                }
                
                if (valid_bins > 0) {
                    h_tensor_product->SetBinContent(i, j, tensor_sum / valid_bins);
                }
            }
        }
        
        delete rho1_norm;
        delete rho2_norm;
    }
    
    void finalize() {
        calculateTensorProduct();
        
        std::string final_pair_name = "h_final_pair_density_" + pair_suffix;
        std::string final_pair_title = "#rho_{2}(#Delta y,#Delta #phi)=dN_{pairs}/dp_{T} " + pair_title + ";#Delta y;#Delta #phi [rad]";
        
        TH2F* h_final_pair_density = (TH2F*)h_pair_density->Clone(final_pair_name.c_str());
        h_final_pair_density->SetTitle(final_pair_title.c_str());
        
        // Normalize by number of events and bin size
        double dy = (config.delta_y_max - config.delta_y_min) / config.n_delta_y_bins;
        double dphi = (config.delta_phi_max - config.delta_phi_min) / config.n_delta_phi_bins;
        if (total_events > 0) {
            h_final_pair_density->Scale(1.0 / (total_events * dy * dphi));
        }
        
        std::string final_tensor_name = "h_final_tensor_product_" + pair_suffix;
        std::string final_tensor_title = "#rho_{1} #otimes #rho_{1}(#Delta y,#Delta #phi) " + pair_title + ";#Delta y;#Delta #phi [rad]";
        TH2F* h_final_tensor = (TH2F*)h_tensor_product->Clone(final_tensor_name.c_str());
        h_final_tensor->SetTitle(final_tensor_title.c_str());
        
        hist_manager->writeHistogram(h_final_pair_density);
        hist_manager->writeHistogram(h_final_tensor);
        hist_manager->writeHistogram(h_pair_density);
        hist_manager->writeHistogram(h_tensor_product);
        hist_manager->writeHistogram(h_delta_y);
        hist_manager->writeHistogram(h_delta_phi);
        hist_manager->writeHistogram(h_pair_multiplicity);
        
        std::cout << "Pair analysis for " << pair_title
                  << ": " << h_pair_multiplicity->GetMean() << " avg pairs/event in "
                  << total_events << " events" << std::endl;
        
        delete h_final_pair_density;
        delete h_final_tensor;
    }
    
private:
    void initializeHistograms() {
        std::string name, title;
        
        name = "h_pair_density_" + pair_suffix;
        title = "#rho_{2}(#Delta y,#Delta #phi) raw " + pair_title + ";#Delta y;#Delta #phi [rad]";
        h_pair_density = hist_manager->createTH2F(name, title,
            config.n_delta_y_bins, config.delta_y_min, config.delta_y_max,
            config.n_delta_phi_bins, config.delta_phi_min, config.delta_phi_max);
        
        name = "h_tensor_product_" + pair_suffix;
        title = "#rho_{1} #otimes #rho_{1}(#Delta y,#Delta #phi) " + pair_title + ";#Delta y;#Delta #phi [rad]";
        h_tensor_product = hist_manager->createTH2F(name, title,
            config.n_delta_y_bins, config.delta_y_min, config.delta_y_max,
            config.n_delta_phi_bins, config.delta_phi_min, config.delta_phi_max);
        
        name = "h_delta_y_" + pair_suffix;
        title = "#Delta y distribution " + pair_title + ";#Delta y;counts";
        h_delta_y = hist_manager->createTH1F(name, title,
            config.n_delta_y_bins, config.delta_y_min, config.delta_y_max);
        
        name = "h_delta_phi_" + pair_suffix;
        title = "#Delta #phi distribution " + pair_title + ";#Delta #phi [rad];counts";
        h_delta_phi = hist_manager->createTH1F(name, title,
            config.n_delta_phi_bins, config.delta_phi_min, config.delta_phi_max);
        
        name = "h_pair_multiplicity_" + pair_suffix;
        title = "Pair multiplicity " + pair_title + ";N_{pairs};events";
        h_pair_multiplicity = hist_manager->createTH1F(name, title, 1000, 0, 1000);
    }
};


class HepMCProcessor {
private:
    std::vector<SingleParticleAnalyzer*> single_analyzers;
    std::vector<PairAnalyzer*> pair_analyzers;
    StatusFilter status_filter;
    
public:
    void setStatusFilter(const std::vector<int>& status_codes) {
        status_filter.setAllowedStatus(status_codes);
    }
    
    void acceptAllStatus() {
        status_filter.acceptAllStatus();
    }
    
    void addSingleParticleAnalyzer(SingleParticleAnalyzer* analyzer) {
        single_analyzers.push_back(analyzer);
    }
    
    void addPairAnalyzer(PairAnalyzer* analyzer) {
        pair_analyzers.push_back(analyzer);
    }
    
    std::vector<int> getStatusCodes() const {
        return status_filter.getStatusCodes();
    }
    
    bool isDefaultStatus() const {
        return status_filter.isDefaultStatus();
    }
    
    SingleParticleAnalyzer* findSingleAnalyzer(int pdg) {
        for (auto* analyzer : single_analyzers) {
            if (analyzer->getPDG() == pdg) {
                return analyzer;
            }
        }
        return nullptr;
    }
    
    void processFile(const std::string& filename) {
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Error opening input file: " << filename << std::endl;
            return;
        }
        
        status_filter.printConfiguration();
        
        std::string line;
        int event_count = 0;
        std::vector<KinematicsCalculator::Particle> event_particles;
        std::map<int, int> event_particle_counts;
        
        while (std::getline(infile, line)) {
            if (line.empty()) continue;
            
            if (line[0] == 'E') {
                if (event_count > 0) {
                    processEvent(event_particles, event_particle_counts, event_count);
                }
                
                event_count++;
                event_particles.clear();
                event_particle_counts.clear();
                
                if (event_count % 1000 == 0) {
                    std::cout << "Processing event " << event_count << std::endl;
                }
                continue;
            }
            
            if (line[0] != 'P') continue;
            
            std::istringstream iss(line);
            char tag;
            int id, pdg, status;
            double px, py, pz, E, m;
            
            iss >> tag >> id >> pdg >> px >> py >> pz >> E >> m >> status;
            
            if (!status_filter.isAllowed(status)) continue;
            
            KinematicsCalculator::Particle particle =
                KinematicsCalculator::createParticle(px, py, pz, E, pdg, event_count);
            
            event_particles.push_back(particle);
            event_particle_counts[pdg]++;
        }
        
        if (event_count > 0) {
            processEvent(event_particles, event_particle_counts, event_count);
        }
        
        infile.close();
        std::cout << "Processed " << event_count << " events" << std::endl;
    }
    
    void finalize() {
        for (auto* analyzer : single_analyzers) {
            analyzer->finalize();
        }
        for (auto* analyzer : pair_analyzers) {
            analyzer->finalize();
        }
    }
    
private:
    void processEvent(const std::vector<KinematicsCalculator::Particle>& particles,
                     const std::map<int, int>& particle_counts, int current_event_id) {
        
        // Process single particles
        for (auto* analyzer : single_analyzers) {
            int count = 0;
            for (const auto& particle : particles) {
                analyzer->processParticle(particle);
                if (particle.pdg == analyzer->getPDG()) {
                    count++;
                }
            }
            analyzer->processEvent(count);
        }
        
        // Process pairs
        for (auto* analyzer : pair_analyzers) {
            analyzer->processPairs(particles);
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " input.hepmc output.root PDG1 PDG2 ... [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --status STATUS1 STATUS2  : Particle status filter (default: 1)" << std::endl;
        std::cerr << "  --all-status              : Accept all particle status codes" << std::endl;
        std::cerr << "\nExample:" << std::endl;
        std::cerr << "  " << argv[0] << " events.hepmc output.root 211 -211 2212" << std::endl;
        std::cerr << "  This will create:" << std::endl;
        std::cerr << "  - Single densities: 211, -211, 2212" << std::endl;
        std::cerr << "  - Pair densities: 211→-211, -211→211, 211→2212, 2212→211, -211→2212, 2212→-211" << std::endl;
        PDGNamer::printNamingConvention();
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
    std::vector<int> pdg_codes;
    std::vector<int> status_codes;
    bool use_all_status = false;
    bool status_specified = false;
    
    // Parse PDG codes and options
    for (int i = 3; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--status") {
            i++;
            while (i < argc && argv[i][0] != '-') {
                try {
                    int status = std::stoi(argv[i]);
                    status_codes.push_back(status);
                } catch (const std::exception& e) {
                    std::cerr << "Invalid status code: " << argv[i] << std::endl;
                    return 1;
                }
                i++;
            }
            i--;
            status_specified = true;
        }
        else if (arg == "--all-status") {
            use_all_status = true;
            status_specified = true;
        }
        else {
            // Try to parse as PDG code
            try {
                int pdg = std::stoi(arg);
                pdg_codes.push_back(pdg);
                std::cout << "Added PDG code: " << pdg << std::endl;
            } catch (const std::exception& e) {
                std::cerr << "Invalid PDG code: " << arg << std::endl;
                return 1;
            }
        }
    }
    
    if (pdg_codes.empty()) {
        std::cerr << "No PDG codes specified!" << std::endl;
        return 1;
    }
    
    TFile* root_file = new TFile(output_file.c_str(), "RECREATE");
    HistogramManager hist_manager(root_file);
    HepMCProcessor processor;
    
    if (use_all_status) {
        processor.acceptAllStatus();
    } else if (status_specified) {
        processor.setStatusFilter(status_codes);
    }
    
    std::vector<int> actual_status_codes = processor.getStatusCodes();
    bool is_default_status = processor.isDefaultStatus();
    
    // Create single particle analyzers for each PDG
    std::vector<SingleParticleAnalyzer*> single_analyzers;
    std::cout << "\nCreating single particle analyzers:" << std::endl;
    for (int pdg : pdg_codes) {
        auto* analyzer = new SingleParticleAnalyzer(pdg, &hist_manager, actual_status_codes, is_default_status);
        single_analyzers.push_back(analyzer);
        processor.addSingleParticleAnalyzer(analyzer);
        std::cout << "  - PDG " << pdg << std::endl;
    }
    
    // Create all possible pair combinations (including XY and YX)
    std::vector<PairAnalyzer*> pair_analyzers;
    std::cout << "\nCreating pair analyzers:" << std::endl;
    for (size_t i = 0; i < pdg_codes.size(); i++) {
        for (size_t j = 0; j < pdg_codes.size(); j++) {
            if (i == j) continue; // Skip self-pairs
            
            int pdg1 = pdg_codes[i];
            int pdg2 = pdg_codes[j];
            
            // Find corresponding single analyzers
            SingleParticleAnalyzer* analyzer1 = processor.findSingleAnalyzer(pdg1);
            SingleParticleAnalyzer* analyzer2 = processor.findSingleAnalyzer(pdg2);
            
            if (!analyzer1 || !analyzer2) {
                std::cerr << "Error: Could not find single analyzers for pair " << pdg1 << "," << pdg2 << std::endl;
                continue;
            }
            
            auto* pair_analyzer = new PairAnalyzer(pdg1, pdg2, &hist_manager, analyzer1, analyzer2, actual_status_codes, is_default_status);
            pair_analyzers.push_back(pair_analyzer);
            processor.addPairAnalyzer(pair_analyzer);
            std::cout << "  - Pair " << pdg1 << " → " << pdg2 << std::endl;
        }
    }
    
    std::cout << "\n" << std::string(50, '=') << std::endl;
    std::cout << "ANALYSIS SUMMARY:" << std::endl;
    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Output file: " << output_file << std::endl;
    std::cout << "PDG codes: ";
    for (size_t i = 0; i < pdg_codes.size(); i++) {
        std::cout << pdg_codes[i];
        if (i < pdg_codes.size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    std::cout << "Single particle analyzers: " << single_analyzers.size() << std::endl;
    std::cout << "Pair analyzers: " << pair_analyzers.size() << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    // Process the file
    processor.processFile(input_file);
    processor.finalize();
    
    // Cleanup
    for (auto* analyzer : single_analyzers) delete analyzer;
    for (auto* analyzer : pair_analyzers) delete analyzer;
    
    root_file->Close();
    delete root_file;
    
    std::cout << "\n" << std::string(50, '=') << std::endl;
    std::cout << "ANALYSIS COMPLETED SUCCESSFULLY!" << std::endl;
    std::cout << "Results saved to: " << output_file << std::endl;
    std::cout << "\nHistogram naming convention:" << std::endl;
    
    if (is_default_status) {
        std::cout << "- Single particle densities: h_final_density_[PDG]" << std::endl;
        std::cout << "- Pair densities (true): h_final_pair_density_[PDG1]_[PDG2]" << std::endl;
        std::cout << "- Tensor products: h_final_tensor_product_[PDG1]_[PDG2]" << std::endl;
        std::cout << "\nNote: 'm' prefix indicates negative PDG codes (e.g., -211 → m211)" << std::endl;
    } else {
        std::cout << "- Single particle densities: h_final_density_[PDG]_status[STATUS]" << std::endl;
        std::cout << "- Pair densities (true): h_final_pair_density_[PDG1]_[PDG2]_status[STATUS]" << std::endl;
        std::cout << "- Tensor products: h_final_tensor_product_[PDG1]_[PDG2]_status[STATUS]" << std::endl;
    }
    
    std::cout << "\nGenerated pairs for PDG codes [";
    for (size_t i = 0; i < pdg_codes.size(); i++) {
        std::cout << pdg_codes[i];
        if (i < pdg_codes.size() - 1) std::cout << ", ";
    }
    std::cout << "]:" << std::endl;
    
    for (size_t i = 0; i < pdg_codes.size(); i++) {
        for (size_t j = 0; j < pdg_codes.size(); j++) {
            if (i != j) {
                std::cout << "  " << pdg_codes[i] << " → " << pdg_codes[j] << std::endl;
            }
        }
    }
    
    std::cout << "\nTo view results:" << std::endl;
    std::cout << "root -l " << output_file << std::endl;
    std::cout << ".ls" << std::endl;
    std::cout << "TBrowser b" << std::endl;
    std::cout << std::string(50, '=') << std::endl;
    
    return 0;
}


