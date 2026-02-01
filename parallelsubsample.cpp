// parallelsubsample.cpp
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TDirectory.h"

std::vector<std::string> findRootFiles(const std::string& baseDir) {
    std::vector<std::string> rootFiles;
    std::string findCommand = "find " + baseDir + " -maxdepth 2 -name '*.root' -type f";
    std::cout << "Executing: " << findCommand << "\n";
    FILE* pipe = popen(findCommand.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Could not execute find command\n";
        return rootFiles;
    }
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        std::string filepath(buffer);
        if (!filepath.empty() && filepath.back() == '\n') {
            filepath.pop_back();
        }
        std::string relative = filepath.substr(baseDir.length() + 1);
        if (relative.find('/') != std::string::npos) {
            rootFiles.push_back(filepath);
        }
    }
    pclose(pipe);
    std::cout << "Found " << rootFiles.size() << " ROOT files:\n";
    for (size_t i = 0; i < std::min(rootFiles.size(), size_t(10)); ++i) {
        std::cout << " " << rootFiles[i] << "\n";
    }
    if (rootFiles.size() > 10) std::cout << " ... and " << (rootFiles.size() - 10) << " more\n";
    return rootFiles;
}

std::vector<std::string> findSubDirectories(const std::string& baseDir) {
    std::vector<std::string> subdirs;
    std::string findCommand = "find " + baseDir + " -maxdepth 2 -name '*.root' -type f -exec dirname {} \\; | sort -u";
    std::cout << "Finding unique subdirectories...\n";
    FILE* pipe = popen(findCommand.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Could not execute find subdirs command\n";
        return subdirs;
    }
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        std::string line(buffer);
        if (!line.empty() && line.back() == '\n') line.pop_back();
        if (!line.empty()) {
            subdirs.push_back(line);
        }
    }
    pclose(pipe);
    std::cout << "Found " << subdirs.size() << " unique subdirectories with ROOT files.\n";
    return subdirs;
}

std::vector<std::string> findRootFilesInDir(const std::string& dirPath) {
    std::vector<std::string> files;
    std::string cmd = "find \"" + dirPath + "\" -maxdepth 1 -name '*.root' -type f";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cerr << "Error: Could not execute find in dir " << dirPath << "\n";
        return files;
    }
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), pipe) != NULL) {
        std::string filepath(buffer);
        if (!filepath.empty() && filepath.back() == '\n') {
            filepath.pop_back();
        }
        if (!filepath.empty()) {
            files.push_back(filepath);
        }
    }
    pclose(pipe);
    return files;
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
    std::cout << "Plots to process (" << plotNames.size() << "):\n";
    for (const auto& name : plotNames) {
        std::cout << " - " << name << "\n";
    }
    return plotNames;
}

bool createOutputDirectory(const std::string& outputPath) {
    std::string mkdirCommand = "mkdir -p \"" + outputPath + "\"";
    std::cout << "Creating output directory: " << outputPath << "\n";
    int result = system(mkdirCommand.c_str());
    if (result != 0) std::cerr << "Failed to create directory\n";
    return (result == 0);
}

TH1* extractHistogram(const std::string& filename, const std::string& plotName,
                      const std::string& subDirPath) {
    std::cout << " Extracting '" << plotName << "' from " << filename << " ... ";
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cout << "FAILED (cannot open)\n";
        return nullptr;
    }
    TH1* hist = nullptr;
    std::string fullPath = subDirPath.empty() ? plotName : subDirPath + "/" + plotName;
    TObject* obj = file->Get(fullPath.c_str());
    if (obj && obj->InheritsFrom("TH1")) {
        hist = static_cast<TH1*>(obj->Clone());
        hist->SetDirectory(nullptr);
        std::cout << "OK (" << (subDirPath.empty() ? "top" : subDirPath) << ")\n";
        delete obj;
        file->Close(); delete file;
        return hist;
    }
    if (obj) delete obj;
    if (!subDirPath.empty()) {
        obj = file->Get(plotName.c_str());
        if (obj && obj->InheritsFrom("TH1")) {
            hist = static_cast<TH1*>(obj->Clone());
            hist->SetDirectory(nullptr);
            std::cout << "OK (top level fallback)\n";
            delete obj;
        }
        if (obj) delete obj;
    }
    if (!hist) std::cout << "NOT FOUND\n";
    file->Close();
    delete file;
    return hist;
}

TH1* computeAverage(const std::vector<TH1*>& histograms, const std::string& plotName) {
    if (histograms.empty()) return nullptr;
    std::cout << "Computing average for " << plotName << " (" << histograms.size() << " hists)\n";
    TH1* avg = static_cast<TH1*>(histograms[0]->Clone());
    avg->SetName(plotName.c_str());  // Keep original plot name for iterative averaging
    std::string title = plotName + " - Average of " + std::to_string(histograms.size()) + " files";
    avg->SetTitle(title.c_str());
    avg->Reset();
    bool is2D = avg->InheritsFrom("TH2");
    int nBinsX = avg->GetNbinsX();
    int nBinsY = is2D ? avg->GetNbinsY() : 1;
    for (int ix = 1; ix <= nBinsX; ++ix) {
        for (int iy = 1; iy <= nBinsY; ++iy) {
            double sum = 0.0;
            for (const TH1* h : histograms) {
                sum += is2D ? h->GetBinContent(ix, iy) : h->GetBinContent(ix);
            }
            double mean = sum / histograms.size();
            double sumsq = 0.0;
            for (const TH1* h : histograms) {
                double v = is2D ? h->GetBinContent(ix, iy) : h->GetBinContent(ix);
                sumsq += (v - mean) * (v - mean);
            }
            double stddev = (histograms.size() > 1) ? std::sqrt(sumsq / (histograms.size() - 1)) : 0.0;
            double err = stddev / std::sqrt(static_cast<double>(histograms.size()));
            if (is2D) {
                avg->SetBinContent(ix, iy, mean);
                avg->SetBinError(ix, iy, err);
            } else {
                avg->SetBinContent(ix, mean);
                avg->SetBinError(ix, err);
            }
        }
    }
    std::cout << " Average computed (integral = " << avg->Integral() << ")\n";
    return avg;
}

int processSinglePlot(const std::vector<std::string>& rootFiles,
                      const std::string& plotName,
                      const std::string& subDirPath,
                      TFile* outputFile) {
    std::cout << "\n┌──────────────────────────────────────\n";
    std::cout << "│ Processing: " << plotName << "\n";
    std::cout << "│ Directory: " << (subDirPath.empty() ? "<top level>" : subDirPath) << "\n";
    std::cout << "└──────────────────────────────────────\n";
    std::vector<TH1*> histograms;
    int found = 0;
    for (size_t i = 0; i < rootFiles.size(); ++i) {
        TH1* h = extractHistogram(rootFiles[i], plotName, subDirPath);
        if (h) {
            histograms.push_back(h);
            ++found;
        }
        if ((i + 1) % 5 == 0 || i + 1 == rootFiles.size()) {
            std::cout << " Progress: " << (i + 1) << "/" << rootFiles.size()
                      << " files checked | found so far: " << found << "\n";
        }
    }
    if (histograms.empty()) {
        std::cout << "→ No histograms found for '" << plotName << "'\n";
        return 0;
    }
    std::cout << "→ Found " << histograms.size() << " valid histograms\n";
    TH1* average = computeAverage(histograms, plotName);
    if (!average) {
        std::cout << "→ Failed to compute average\n";
        for (auto* h : histograms) delete h;
        return 0;
    }
    if (outputFile && !outputFile->IsZombie()) {
        outputFile->cd();
        average->Write("", TObject::kOverwrite);
        outputFile->Flush();
        std::cout << "→ Written & flushed: " << average->GetName() << "\n";
    }
    for (auto* h : histograms) delete h;
    return 1;
}

int averagePlots(const std::string& baseDirectory,
                 const std::string& plotNamesStr,
                 const std::string& subDirPath,
                 const std::string& outputPath = "./output",
                 const std::string& outputFileName = "averaged_plots.root",
                 int numJobs = 1,
                 int jobId = 0) {
    std::cout << "\n=== Starting plot averaging ===\n";
    std::cout << "Base directory: " << baseDirectory << "\n";
    std::cout << "Target subdir: " << (subDirPath.empty() ? "<top level>" : subDirPath) << "\n";
    std::cout << "Output: " << outputPath << "/" << outputFileName << "\n";
    std::cout << "Mode: " << (numJobs > 1 ? "Parallel (Job " + std::to_string(jobId) + "/" + std::to_string(numJobs) + ")" : "Single") << "\n\n";

    if (!createOutputDirectory(outputPath)) return 1;

    auto plotNames = parseplotNames(plotNamesStr);
    if (plotNames.empty()) {
        std::cerr << "No valid plot names provided\n";
        return 1;
    }

    std::vector<std::string> rootFiles;
    if (numJobs <= 1) {
        rootFiles = findRootFiles(baseDirectory);
    } else {
        auto subdirs = findSubDirectories(baseDirectory);
        if (subdirs.empty()) {
            std::cerr << "No subdirectories found\n";
            return 1;
        }
        std::sort(subdirs.begin(), subdirs.end());
        int total = subdirs.size();
        int chunk_size = (total + numJobs - 1) / numJobs;
        int start = jobId * chunk_size;
        int end = std::min(start + chunk_size, total);
        std::vector<std::string> assigned_subdirs(subdirs.begin() + start, subdirs.begin() + end);
        std::cout << "Job " << jobId << "/" << numJobs << " assigned " << assigned_subdirs.size() 
                  << " subdirs (" << start << " to " << (end-1) << ")\n";
        for (const auto& sd : assigned_subdirs) {
            auto local_files = findRootFilesInDir(sd);
            rootFiles.insert(rootFiles.end(), local_files.begin(), local_files.end());
        }
        std::cout << "Collected " << rootFiles.size() << " ROOT files for this job\n";
    }

    if (rootFiles.empty()) {
        std::cerr << "No ROOT files found for this job\n";
        return 1;
    }

    std::string fullPath = outputPath + "/" + outputFileName;
    std::cout << "Opening output file: " << fullPath << " (UPDATE mode)\n";
    TFile* output = TFile::Open(fullPath.c_str(), "UPDATE");
    if (!output || output->IsZombie()) {
        std::cerr << "Failed to open/create output file\n";
        return 1;
    }

    int success = 0;
    for (const auto& plot : plotNames) {
        success += processSinglePlot(rootFiles, plot, subDirPath, output);
    }

    output->Flush();
    output->Close();
    delete output;

    std::cout << "\nFinished processing " << success << " / " << plotNames.size() << " plots.\n";
    std::cout << "=== Done ===\n";
    return (success > 0) ? 0 : 1;
}

#ifndef __CINT__
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <base_directory> \"plot1,plot2,...\" <subdirectory> [output_path] [output_filename] [num_jobs] [job_id]\n";
        std::cerr << "  For parallel: provide num_jobs and job_id (0-based)\n";
        std::cerr << "  Example subdirectory: \"Class1/SymmetricBalanceFunctions_Bs\" or \"\" for top level\n";
        return 1;
    }
    std::string baseDir = argv[1];
    std::string plots = argv[2];
    std::string subDir = argv[3];
    std::string outPath = (argc > 4) ? argv[4] : "./output";
    std::string outFile = (argc > 5) ? argv[5] : "averaged_plots.root";
    int num_jobs = (argc > 6) ? std::stoi(argv[6]) : 1;
    int job_id = (argc > 7) ? std::stoi(argv[7]) : 0;

    return averagePlots(baseDir, plots, subDir, outPath, outFile, num_jobs, job_id);
}
#endif
