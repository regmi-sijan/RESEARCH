#include "TFile.h"
#include "TH3.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TStyle.h"
#include "TH1D.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

void ExpPol2Fit() {
    TFile *file = TFile::Open("combine_run24_pp_photrig_evtmix.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: cannot open input ROOT file!" << std::endl;
        return;
    }

    TH3F *h3 = dynamic_cast<TH3F*>(file->Get("pairInvMassPtEta"));
    if (!h3) {
        std::cerr << "Error: histogram 'pairInvMassPtEta' not found!" << std::endl;
        return;
    }

    std::ifstream csv("ExcludedRanges.csv");
    if (!csv.is_open()) {
        std::cerr << "Error: could not open outliers.csv!" << std::endl;
        return;
    }

    std::string line;
    std::getline(csv, line); // skip header

    int index = 0;
    while (std::getline(csv, line)) {
        std::stringstream ss(line);
        double ptMin, ptMax, leftOutlierX, rightOutlierX, fitMean;
        char delim;

        ss >> ptMin >> delim >> ptMax >> delim >> leftOutlierX >> delim >> rightOutlierX >> delim >> fitMean;

        int ptBinMin = h3->GetYaxis()->FindBin(ptMin + 0.0001);
        int ptBinMax = h3->GetYaxis()->FindBin(ptMax - 0.0001);

        TH1D *h_proj = h3->ProjectionX(Form("h_proj_%d", index), ptBinMin, ptBinMax, 1, h3->GetZaxis()->GetNbins());
        h_proj->SetTitle(Form("InvMass (%.0f < pT < %.0f); Mass [GeV/c^{2}]; Counts", ptMin, ptMax));

        int leftOutBin = h_proj->FindBin(leftOutlierX);
        int rightOutBin = h_proj->FindBin(rightOutlierX);

        // Define fitting range: 10 bins left of leftOutBin and 20 bins right of rightOutBin
        int fitStartBin = std::max(1, leftOutBin - 10);
        int fitEndBin = std::min(h_proj->GetNbinsX(), rightOutBin + 20);

        // Clone and mask excluded region
        TH1D *h_fit = (TH1D*)h_proj->Clone(Form("h_fit_%d", index));
        for (int b = 1; b <= h_proj->GetNbinsX(); ++b) {
            if (b < fitStartBin || b > fitEndBin || (b >= leftOutBin && b <= rightOutBin)) {
                h_fit->SetBinContent(b, 0);
                h_fit->SetBinError(b, 0);
            }
        }

        // Fit: exp + pol2
        TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-[1]*x) + [2] + [3]*x + [4]*x*x",
                               h_proj->GetXaxis()->GetXmin(), h_proj->GetXaxis()->GetXmax());
        fitFunc->SetParameters(h_fit->GetMaximum(), 1.0, 0.0, 0.0, 0.0);
        h_fit->Fit(fitFunc, "RQ0");

        // Extrapolate into excluded region
        TH1D* h_extrapolated = (TH1D*)h_proj->Clone(Form("h_extrapolated_%d", index));
        h_extrapolated->Reset();
        for (int b = leftOutBin; b <= rightOutBin; ++b) {
            double x = h_extrapolated->GetBinCenter(b);
            double val = fitFunc->Eval(x);
            h_extrapolated->SetBinContent(b, val);
        }

        h_extrapolated->SetFillColor(kRed);
        h_extrapolated->SetLineColor(kRed);
        h_extrapolated->SetFillStyle(1001);

        // === Draw: hist → fit → extrapolated → lines → hist again
        TCanvas *c = new TCanvas(Form("c_%d", index), "c", 800, 600);
        h_proj->SetLineColor(kBlack);
        h_proj->Draw("E");

        fitFunc->SetLineColor(kBlue);
        fitFunc->SetLineWidth(2);
        fitFunc->Draw("SAME");

        h_extrapolated->Draw("HIST SAME");

        TLine *l1 = new TLine(leftOutlierX, 0, leftOutlierX, h_proj->GetMaximum());
        TLine *l2 = new TLine(rightOutlierX, 0, rightOutlierX, h_proj->GetMaximum());
        l1->SetLineStyle(2); l2->SetLineStyle(2);
        l1->SetLineColor(kGray + 2); l2->SetLineColor(kGray + 2);
        l1->Draw("SAME");
        l2->Draw("SAME");

        h_proj->Draw("E SAME");  // Final draw to bring histogram on top

        c->SaveAs(Form("invMassFit_expPol2_global_%d.png", index));

        // Cleanup
        delete c;
        delete h_proj;
        delete h_fit;
        delete h_extrapolated;
        delete fitFunc;
        delete l1;
        delete l2;

        ++index;
    }

    csv.close();
    file->Close();
}

