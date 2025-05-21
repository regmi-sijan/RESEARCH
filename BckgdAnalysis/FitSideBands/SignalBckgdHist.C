#include "TFile.h"
#include "TH3.h"
#include "TH1D.h"
#include "TF1.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

void SignalBckgdHist() {

    // --- Open input ROOT file and 3D histogram
    TFile *inFile = TFile::Open("combine_run24_pp_photrig_evtmix.root");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: cannot open input ROOT file!" << std::endl;
        return;
    }
    TH3F *h3 = dynamic_cast<TH3F*>(inFile->Get("pairInvMassPtEta"));
    if (!h3) {
        std::cerr << "Error: histogram 'pairInvMassPtEta' not found!" << std::endl;
        inFile->Close();
        return;
    }

    // --- Create output ROOT file
    TFile *outFile = TFile::Open("expPol2Fits.root", "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: cannot open output ROOT file!" << std::endl;
        inFile->Close();
        return;
    }

    // --- Open CSV of excluded ranges
    std::ifstream csvIn("ExcludedRanges.csv");
    if (!csvIn.is_open()) {
        std::cerr << "Error: could not open ExcludedRanges.csv!" << std::endl;
        inFile->Close();
        outFile->Close();
        return;
    }
    std::string line;
    std::getline(csvIn, line); // skip header

    // --- Prepare output CSV for results
    std::ofstream csvOut("Results.csv");
    csvOut << "ptMin,ptMax,mean,sigma,yield" << std::endl;

    int index = 0;
    while (std::getline(csvIn, line)) {
        std::stringstream ss(line);
        double ptMin, ptMax, leftOutX, rightOutX, fitMean;
        char delim;
        ss >> ptMin >> delim
           >> ptMax >> delim
           >> leftOutX >> delim
           >> rightOutX >> delim
           >> fitMean;

        // --- Determine pT-bin range
        int ptBinMin = h3->GetYaxis()->FindBin(ptMin + 1e-4);
        int ptBinMax = h3->GetYaxis()->FindBin(ptMax - 1e-4);

        // --- Project onto mass axis and name fg_<ptMin>_<ptMax>
        TString fgName = Form("fg_%.0f_%.0f", ptMin, ptMax);
        TH1D *h_fg = h3->ProjectionX(
            fgName,
            ptBinMin, ptBinMax,
            1, h3->GetZaxis()->GetNbins()
        );
        h_fg->SetTitle(Form("InvMass (%.0f<pT<%.0f); Mass [GeV]; Counts", ptMin, ptMax));

        // --- Identify excluded-region bins and fitting window
        int leftBin  = h_fg->FindBin(leftOutX);
        int rightBin = h_fg->FindBin(rightOutX);
        int fitStart = std::max(1,                   leftBin - 10);
        int fitEnd   = std::min(h_fg->GetNbinsX(),   rightBin + 20);

        // --- Clone & mask for background fit
        TH1D *h_fit = (TH1D*)h_fg->Clone(Form("h_fit_%02d", index));
        for (int b = 1; b <= h_fg->GetNbinsX(); ++b) {
            if (b < fitStart || b > fitEnd || (b >= leftBin && b <= rightBin)) {
                h_fit->SetBinContent(b, 0);
                h_fit->SetBinError(b,   0);
            }
        }

        // --- Fit function: exponential + quadratic
        TF1 *f_bg = new TF1(
            Form("f_bg_%02d", index),
            "[0]*exp(-[1]*x) + [2] + [3]*x + [4]*x*x",
            h_fg->GetXaxis()->GetXmin(),
            h_fg->GetXaxis()->GetXmax()
        );
        f_bg->SetParameters(h_fit->GetMaximum(), 1.0, 0.0, 0.0, 0.0);
        h_fit->Fit(f_bg, "RQ0");

        // --- Build background histogram bg_<ptMin>_<ptMax>
        TString bgName = Form("bg_%.0f_%.0f", ptMin, ptMax);
        TH1D *h_bg = (TH1D*)h_fg->Clone(bgName);
        h_bg->Reset();
        for (int b = fitStart; b <= fitEnd; ++b) {
            double x = h_bg->GetBinCenter(b);
            double y = f_bg->Eval(x);
            h_bg->SetBinContent(b, y);
        }

        // --- Subtract background from foreground only in fit window
        TH1D *h_sig = (TH1D*)h_fg->Clone(Form("sig_%.0f_%.0f", ptMin, ptMax));
        for (int b = 1; b <= h_sig->GetNbinsX(); ++b) {
            if (b >= fitStart && b <= fitEnd) {
                double sig = h_fg->GetBinContent(b) - h_bg->GetBinContent(b);
                h_sig->SetBinContent(b, sig);
            } else {
                h_sig->SetBinContent(b, 0);
                h_sig->SetBinError(b,   0);
            }
        }

        // --- Find global peak between 0.45 and 0.65 GeV
        int binLow  = h_sig->FindBin(0.45);
        int binHigh = h_sig->FindBin(0.65);
        double maxVal = -1e9;
        int peakBin = binLow;
        for (int b = binLow; b <= binHigh; ++b) {
            if (h_sig->GetBinContent(b) > maxVal) {
                maxVal = h_sig->GetBinContent(b);
                peakBin = b;
            }
        }
        
        // --- Fit Gaussian ±10 bins around peak
        int gaussStartBin = std::max(1, peakBin - 5);
        int gaussEndBin   = std::min(h_sig->GetNbinsX(), peakBin + 5);
        double xLo = h_sig->GetBinCenter(gaussStartBin);
        double xHi = h_sig->GetBinCenter(gaussEndBin);

        TF1 *f_gaus = new TF1(
            Form("f_gaus_%02d", index),
            "gaus",
            xLo, xHi
        );
        f_gaus->SetParameters(maxVal, h_sig->GetBinCenter(peakBin), 0.01);
        h_sig->Fit(f_gaus, "RQ0");

        double mean  = f_gaus->GetParameter(1);
        double sigma = f_gaus->GetParameter(2);

        // --- Integrate ±3σ around peak
        int intBinLo = h_sig->FindBin(mean - 3*sigma);
        int intBinHi = h_sig->FindBin(mean + 3*sigma);
        double yield = h_sig->Integral(intBinLo, intBinHi);

        // --- Write histograms to output file
        outFile->cd();
        h_fg->Write();
        h_bg->Write();
        h_sig->Write();

        // --- Append results to CSV
        csvOut << ptMin << "," << ptMax << ","
               << mean  << "," << sigma  << "," << yield << std::endl;

        // --- Cleanup
        delete h_fg;
        delete h_fit;
        delete h_bg;
        delete h_sig;
        delete f_bg;
        delete f_gaus;

        ++index;
    }

    // --- Close files
    csvIn.close();
    csvOut.close();
    outFile->Close();
    inFile->Close();

    std::cout << "Processed " << index << " pT bins."
              << " Results in expPol2Fits.root and Results.csv\n";
}

