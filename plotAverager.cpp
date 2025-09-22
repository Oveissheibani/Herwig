#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <cstring>

#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TClass.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"

std::vector<std::string> findRootFiles(const std::string& baseDir) {
    std::vector<std::string> rootFiles;
    
    std::string findCommand = "find " + baseDir + " -maxdepth 2 -name '*.root' -type f";
    FILE* pipe = popen(findCommand.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error running find command" << std::endl;
        return rootFiles;
    }
    
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        std::string filepath(buffer);
        if (!filepath.empty() && filepath[filepath.length()-1] == '\n') {
            filepath.erase(filepath.length()-1);
        }
        
        std::string relativePath = filepath.substr(baseDir.length() + 1);
        if (relativePath.find('/') != std::string::npos) {
            rootFiles.push_back(filepath);
            std::cout << "Found ROOT file: " << filepath << std::endl;
        }
    }
    pclose(pipe);
    
    return rootFiles;
}

std::vector<std::string> parseplotNames(const std::string& plotNamesStr) {
    std::vector<std::string> plotNames;
    std::stringstream ss(plotNamesStr);
    std::string item;
    
    while (std::getline(ss, item, ',')) {
        item.erase(0, item.find_first_not_of(" \t"));
        item.erase(item.find_last_not_of(" \t") + 1);
        if (!item.empty()) {
            plotNames.push_back(item);
        }
    }
    
    return plotNames;
}

bool createOutputDirectory(const std::string& outputPath) {
    std::string mkdirCommand = "mkdir -p " + outputPath;
    int result = system(mkdirCommand.c_str());
    if (result == 0) {
        std::cout << "Created output directory: " << outputPath << std::endl;
        return true;
    } else {
        std::cerr << "Error creating output directory: " << outputPath << std::endl;
        return false;
    }
}

TH1* searchInDirectory(TDirectory* dir, const std::string& plotName) {
    TIter next(dir->GetListOfKeys());
    TKey* key;
    while ((key = static_cast<TKey*>(next()))) {
        TObject* keyObj = key->ReadObj();
        
        if (keyObj->InheritsFrom("TDirectory")) {
            TH1* found = searchInDirectory(static_cast<TDirectory*>(keyObj), plotName);
            if (found) return found;
        } else if (keyObj->InheritsFrom("TH1")) {
            if (std::string(keyObj->GetName()) == plotName) {
                TH1* foundHist = static_cast<TH1*>(keyObj->Clone());
                foundHist->SetDirectory(0);
                return foundHist;
            }
        }
    }
    return nullptr;
}

TH1* extractHistogram(const std::string& filename, const std::string& plotName) {
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return nullptr;
    }
    
    TH1* hist = nullptr;
    
    TObject* obj = file->Get(plotName.c_str());
    if (obj && obj->InheritsFrom("TH1")) {
        hist = static_cast<TH1*>(obj->Clone());
        hist->SetDirectory(0);
    } else {
        hist = searchInDirectory(file, plotName);
    }
    
    file->Close();
    delete file;
    
    if (!hist) {
        std::cerr << "Histogram '" << plotName << "' not found in file: " << filename << std::endl;
    }
    
    return hist;
}

TH1* computeAverage(const std::vector<TH1*>& histograms, const std::string& plotName) {
    if (histograms.empty()) {
        std::cerr << "No histograms to average!" << std::endl;
        return nullptr;
    }
    
    TH1* avgHist = static_cast<TH1*>(histograms[0]->Clone());
    avgHist->SetName((plotName + "_average").c_str());
    avgHist->SetTitle((plotName + " - Average of " + std::to_string(histograms.size()) + " files").c_str());
    avgHist->Reset();
    
    int nBins = avgHist->GetNbinsX();
    int nHists = histograms.size();
    
    bool is2D = avgHist->InheritsFrom("TH2");
    int nBinsY = is2D ? avgHist->GetNbinsY() : 1;
    
    for (int i = 1; i <= nBins; i++) {
        for (int j = 1; j <= nBinsY; j++) {
            std::vector<double> binValues;
            
            for (const auto& hist : histograms) {
                double binContent = is2D ? hist->GetBinContent(i, j) : hist->GetBinContent(i);
                binValues.push_back(binContent);
            }
            
            double sum = 0.0;
            for (double val : binValues) {
                sum += val;
            }
            double mean = sum / nHists;
            
            double sumSquares = 0.0;
            for (double val : binValues) {
                sumSquares += (val - mean) * (val - mean);
            }
            double stdDev = (nHists > 1) ? std::sqrt(sumSquares / (nHists - 1)) : 0.0;
            double stdError = stdDev / std::sqrt(nHists);
            
            if (is2D) {
                avgHist->SetBinContent(i, j, mean);
                avgHist->SetBinError(i, j, stdError);
            } else {
                avgHist->SetBinContent(i, mean);
                avgHist->SetBinError(i, stdError);
            }
        }
    }
    
    return avgHist;
}

int processSinglePlot(const std::vector<std::string>& rootFiles, const std::string& plotName, 
                      TFile* outputFile, const std::string& outputPath) {
    std::cout << "\nProcessing plot: " << plotName << std::endl;
    
    std::vector<TH1*> histograms;
    for (const auto& filename : rootFiles) {
        TH1* hist = extractHistogram(filename, plotName);
        if (hist) {
            histograms.push_back(hist);
        }
    }
    
    if (histograms.empty()) {
        std::cerr << "No histograms found with name: " << plotName << std::endl;
        return 0;
    }
    
    std::cout << "Successfully loaded " << histograms.size() << " histograms for " << plotName << std::endl;
    
    TH1* avgHistogram = computeAverage(histograms, plotName);
    if (!avgHistogram) {
        std::cerr << "Failed to compute average histogram for: " << plotName << std::endl;
        for (auto* hist : histograms) {
            delete hist;
        }
        return 0;
    }
    
    if (outputFile && !outputFile->IsZombie()) {
        outputFile->cd();
        avgHistogram->Write();
        
        for (size_t i = 0; i < histograms.size(); i++) {
            std::string name = plotName + "_file" + std::to_string(i);
            histograms[i]->SetName(name.c_str());
            histograms[i]->Write();
        }
    }
    
    gStyle->SetOptStat(1111);
    TCanvas* canvas = new TCanvas(("canvas_" + plotName).c_str(), 
                                 ("Average Plot: " + plotName).c_str(), 800, 600);
    
    avgHistogram->SetLineColor(kBlue);
    avgHistogram->SetMarkerColor(kBlue);
    avgHistogram->SetMarkerStyle(20);
    avgHistogram->Draw("E1");
    
    std::string pngPath = outputPath + "/" + plotName + "_average.png";
    std::string pdfPath = outputPath + "/" + plotName + "_average.pdf";
    
    canvas->SaveAs(pngPath.c_str());
    canvas->SaveAs(pdfPath.c_str());
    
    std::cout << "Plots saved to: " << outputPath << "/" << plotName << "_average.png and .pdf" << std::endl;
    
    std::cout << "=== Summary for " << plotName << " ===" << std::endl;
    std::cout << "Number of files processed: " << histograms.size() << std::endl;
    std::cout << "Average histogram entries: " << avgHistogram->GetEntries() << std::endl;
    std::cout << "Average histogram mean: " << avgHistogram->GetMean() << std::endl;
    std::cout << "Average histogram RMS: " << avgHistogram->GetRMS() << std::endl;
    
    delete canvas;
    delete avgHistogram;
    for (auto* hist : histograms) {
        delete hist;
    }
    
    return 1;
}

int averagePlots(const std::string& baseDirectory, const std::string& plotNamesStr, 
                 const std::string& outputPath = "./output", const std::string& outputFileName = "multi_plot_averaged.root") {
    
    std::cout << "=== ROOT Multi-Plot Averaging Script ===" << std::endl;
    std::cout << "Base directory: " << baseDirectory << std::endl;
    std::cout << "Plot names: " << plotNamesStr << std::endl;
    std::cout << "Output path: " << outputPath << std::endl;
    std::cout << "Output file: " << outputFileName << std::endl;
    std::cout << "========================================" << std::endl;
    
    if (!createOutputDirectory(outputPath)) {
        return 1;
    }
    
    std::vector<std::string> plotNames = parseplotNames(plotNamesStr);
    if (plotNames.empty()) {
        std::cerr << "No valid plot names provided!" << std::endl;
        return 1;
    }
    
    std::cout << "Will process " << plotNames.size() << " plots:" << std::endl;
    for (const auto& name : plotNames) {
        std::cout << "  - " << name << std::endl;
    }
    
    std::vector<std::string> rootFiles = findRootFiles(baseDirectory);
    
    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found in subdirectories of: " << baseDirectory << std::endl;
        return 1;
    }
    
    std::cout << "Found " << rootFiles.size() << " ROOT files" << std::endl;
    
    std::string fullOutputPath = outputPath + "/" + outputFileName;
    TFile* output = TFile::Open(fullOutputPath.c_str(), "RECREATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Error creating output file: " << fullOutputPath << std::endl;
        return 1;
    }
    
    int successCount = 0;
    for (const auto& plotName : plotNames) {
        successCount += processSinglePlot(rootFiles, plotName, output, outputPath);
    }
    
    output->Close();
    delete output;
    
    std::cout << "\n=== Final Summary ===" << std::endl;
    std::cout << "Total plots requested: " << plotNames.size() << std::endl;
    std::cout << "Successfully processed: " << successCount << std::endl;
    std::cout << "Results saved to: " << fullOutputPath << std::endl;
    std::cout << "Plot images saved to: " << outputPath << std::endl;
    
    return (successCount > 0) ? 0 : 1;
}

void runAveraging() {
    std::string baseDir = "./output";
    std::string plotNames = "h_pt_211_status1,h_eta_211_status1,h_phi_211_status1";
    std::string outputPath = "./averaged_results";
    std::string outputFileName = "multi_averaged_results.root";
    
    averagePlots(baseDir, plotNames, outputPath, outputFileName);
}

#ifndef __CINT__
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: " << argv[0] << " <base_directory> <plot_names> [output_path] [output_filename]" << std::endl;
        std::cout << "Plot names should be comma-separated (no spaces around commas)" << std::endl;
        std::cout << "Default output_path: ./output" << std::endl;
        std::cout << "Default output_filename: multi_plot_averaged.root" << std::endl;
        std::cout << "Examples:" << std::endl;
        std::cout << "  " << argv[0] << " ./output \"h_pt_211_status1,h_eta_211_status1\"" << std::endl;
        std::cout << "  " << argv[0] << " ./output \"h_pt_211_status1\" ./results" << std::endl;
        std::cout << "  " << argv[0] << " ./output \"h_pt_211_status1\" ./results my_analysis.root" << std::endl;
        return 1;
    }
    
    std::string baseDir = argv[1];
    std::string plotNames = argv[2];
    std::string outputPath = (argc > 3) ? argv[3] : "./output";
    std::string outputFileName = (argc > 4) ? argv[4] : "multi_plot_averaged.root";
    
    return averagePlots(baseDir, plotNames, outputPath, outputFileName);
}
#endif
