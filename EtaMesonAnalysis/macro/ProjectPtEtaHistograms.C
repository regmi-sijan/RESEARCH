#include <TFile.h>
#include <TH3.h>
#include <TH1.h>
#include <TString.h>
#include <TSystem.h>
#include <vector>
#include <iostream>

void ProjectPtEtaHistograms(bool doProjectPt = true, bool doProjectEta = false)
{
    const std::string input_file = "combined.root";
    const std::string output_file = "proj_" + input_file;

    // Hardcoded pT ranges
    std::vector<double> ptLows  = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0};
    std::vector<double> ptHighs = {4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 18.0};

    // Hardcoded eta ranges
    std::vector<double> etaLows  = {-1.1, -0.7, -0.3, 0.0, 0.3, 0.7};
    std::vector<double> etaHighs = {-0.7, -0.3, 0.0, 0.3, 0.7, 1.1};

    if (ptLows.size() != ptHighs.size() || etaLows.size() != etaHighs.size()) {
        std::cerr << "Error: pT and/or eta range vectors don't match in size!" << std::endl;
        return;
    }

    // Open input and output files
    TFile *file = new TFile(input_file.c_str(), "READ");
    TFile *outFile = new TFile(output_file.c_str(), "RECREATE");

    TH3F *hist_f = (TH3F*)file->Get("pairInvMassPtEta");
    TH3F *hist_b = (TH3F*)file->Get("pairInvMassPtEtaBkgd");

    if (!hist_f || !hist_b) {
        std::cerr << "Error: Could not find one or both histograms!" << std::endl;
        return;
    }

    // ==== Project on pT ranges ====
    if (doProjectPt) {
        for (size_t i = 0; i < ptLows.size(); ++i) {
            int yMin_f = hist_f->GetYaxis()->FindBin(ptLows[i]);
            int yMax_f = hist_f->GetYaxis()->FindBin(ptHighs[i]);
            int zMin_f = hist_f->GetZaxis()->FindBin(hist_f->GetZaxis()->GetXmin());
            int zMax_f = hist_f->GetZaxis()->FindBin(hist_f->GetZaxis()->GetXmax());

            int yMin_b = hist_b->GetYaxis()->FindBin(ptLows[i]);
            int yMax_b = hist_b->GetYaxis()->FindBin(ptHighs[i]);
            int zMin_b = hist_b->GetZaxis()->FindBin(hist_b->GetZaxis()->GetXmin());
            int zMax_b = hist_b->GetZaxis()->FindBin(hist_b->GetZaxis()->GetXmax());

            TH1* h_f = hist_f->ProjectionX(Form("f%zu", i), yMin_f, yMax_f-1, zMin_f, zMax_f, "C");
            TH1* h_b = hist_b->ProjectionX(Form("b%zu", i), yMin_b, yMax_b-1, zMin_b, zMax_b, "C");

            //TH1* h_f = hist_f->ProjectionX(Form("pT_f_%zu", i), yMin_f, yMax_f-1, zMin_f, zMax_f, "C");
            //TH1* h_b = hist_b->ProjectionX(Form("pT_b_%zu", i), yMin_b, yMax_b-1, zMin_b, zMax_b, "C");
            
						h_f->SetTitle(Form("Mass projection for pT [%g, %g]", ptLows[i], ptHighs[i]));
            h_b->SetTitle(Form("Mass projection for pT [%g, %g]", ptLows[i], ptHighs[i]));

            h_f->Write();
            h_b->Write();
        }
    }

    // ==== Project on eta ranges ====
    if (doProjectEta) {
        for (size_t i = 0; i < etaLows.size(); ++i) {
            int yMin_f = hist_f->GetYaxis()->FindBin(hist_f->GetYaxis()->GetXmin());
            int yMax_f = hist_f->GetYaxis()->FindBin(hist_f->GetYaxis()->GetXmax());
            int zMin_f = hist_f->GetZaxis()->FindBin(etaLows[i]);
            int zMax_f = hist_f->GetZaxis()->FindBin(etaHighs[i]);

            int yMin_b = hist_b->GetYaxis()->FindBin(hist_b->GetYaxis()->GetXmin());
            int yMax_b = hist_b->GetYaxis()->FindBin(hist_b->GetYaxis()->GetXmax());
            int zMin_b = hist_b->GetZaxis()->FindBin(etaLows[i]);
            int zMax_b = hist_b->GetZaxis()->FindBin(etaHighs[i]);

            TH1* h_f = hist_f->ProjectionX(Form("eta_f_%zu", i), yMin_f, yMax_f, zMin_f, zMax_f-1, "C");
            TH1* h_b = hist_b->ProjectionX(Form("eta_b_%zu", i), yMin_b, yMax_b, zMin_b, zMax_b-1, "C");

            h_f->SetTitle(Form("Mass projection for eta [%g, %g]", etaLows[i], etaHighs[i]));
            h_b->SetTitle(Form("Mass projection for eta [%g, %g]", etaLows[i], etaHighs[i]));

            h_f->Write();
            h_b->Write();
        }
    }

    outFile->Close();
    file->Close();
    delete file;
    delete outFile;

    std::cout << "Projection complete!" << std::endl;
    gSystem->Exit(0);
}


