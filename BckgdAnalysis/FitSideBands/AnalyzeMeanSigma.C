#include <TFile.h>
#include <TH3F.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TBox.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>


// this code will give us tentative range where eta-meson peak will be located
// we will get mean and sigma for eta meson peak (not dedicated towards pi0 now but can be used with simple modification and testing)
// output *.csv file have the range to extrapolate, but is soft range and must fine-tune manually before using in subsequent analysis

// we are supposed to check green color range in the *.png images
// they are exactly same as the range we save to extrapolate in .csv file
// we can manually adjust .csv file looking at the images it produce by this code
// such check will be first level check

// all the necessary parameters
static const std::string InputFile = "combine_run24_pp_photrig_evtmix.root";
static const std::string HistName  = "pairInvMassPtEta";
static const std::string OutputCSV = "ExcludedRanges.csv";
static const double seedMass = 0.55;
static const double seedSigma = 0.02;

// Define pair pT ranges here: {pT_min, pT_max}
static const std::vector<std::pair<double,double>> pTRanges = {
    {2.0, 4.0}, {4.0, 6.0}, {6.0, 9.0}, {9.0,12.0}, {12.0,15.0}, {15.0,22.0}
	};

// helper function to perform dynamic-window Gaussian fit around a seed mean
// this will be useful to get the best range for gaussian fit
bool PerformDynamicGaussianFit(TH1D* hist, int seedBin, double seedMean, double seedSigma, double& outMean, double& outSigma)
{
    const int nbins = hist->GetNbinsX();
    double bestChi2 = 1e9;
    bool success = false;
		
		// using window of 3 to 8 (this is  why we call this as dynamic window gauss fit
    for (int win = 3; win <= 8; ++win) {
        
				// lower bin and higher bin 
				int bL = std::max(1, seedBin - win);
        int bH = std::min(nbins, seedBin + win);
        
				if (hist->Integral(bL, bH) < 10) continue;
				
				// lower and higher x values in hist
        double xL = hist->GetXaxis()->GetBinLowEdge(bL);
        double xH = hist->GetXaxis()->GetBinUpEdge(bH);
        
				TF1 fit("fit", "gaus", xL, xH);
        fit.SetParameters(hist->GetBinContent(seedBin), seedMean, seedSigma);
        
				if (hist->Fit(&fit, "RQ0") != 0) continue;

        double m = fit.GetParameter(1);
        double s = fit.GetParameter(2);
        double chi2ndf = fit.GetChisquare() / fit.GetNDF();
        
				// to select best fit
				if (chi2ndf < bestChi2 && s > 0) {
            bestChi2 = chi2ndf;
            outMean = m;
            outSigma = s;
            success = true;
        }
    }

    return success;
}

// Compute median absolute deviation
static double ComputeMAD(std::vector<double>& vals) {
    
		if (vals.empty()) return 0;
    
		std::sort(vals.begin(), vals.end());
    
		double med = vals[vals.size()/2];
    
		std::vector<double> absdev;
    absdev.reserve(vals.size());
    
		for (double v : vals) absdev.push_back(std::fabs(v - med));
    
		std::sort(absdev.begin(), absdev.end());
    
		return absdev[absdev.size()/2];
}

// main function to analyze mean and sigma that will give us tentative range for extrapolation and fitting side bands
void AnalyzeMeanSigma()
{
		std::ofstream csvOut(OutputCSV);
    if (!csvOut) {
        std::cerr << "Error: cannot open outliers.csv for writing.\n";
        return;
    }
    csvOut << "pTmin,pTmax,ExcludeLeft,ExcludeRight,fitMean\n";

    TFile* file = TFile::Open(InputFile.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: can't open file.\n";
        return;
    }
    TH3F* h3 = dynamic_cast<TH3F*>(file->Get(HistName.c_str()));
    if (!h3) {
        std::cerr << "Histogram not found: " << HistName << std::endl;
        file->Close();
        return;
    }

    auto* massAxis = h3->GetXaxis();
    auto* ptAxis   = h3->GetYaxis();
    auto* etaAxis  = h3->GetZaxis();
		
		// now we will be looking at different ranges
    for (size_t i = 0; i < pTRanges.size(); ++i) {
        
				float p1 = pTRanges[i].first, p2 = pTRanges[i].second;

        int b1 = ptAxis->FindBin(p1), b2 = ptAxis->FindBin(p2) - 1;
        
				TH1D* hProj = h3->ProjectionX(Form("mass_pT_%zu", i), b1, b2, 1, etaAxis->GetNbins());
        
				if (hProj->GetEntries() <= 10) { delete hProj; continue; }

        // Perform two-pass Gaussian fit
        int seedBin = massAxis->FindBin(seedMass);
        double m1 = 0, s1 = 0;
        
				// first pass will help to roughly locate the peak location
				PerformDynamicGaussianFit(hProj, seedBin, seedMass, 0.02, m1, s1);
        
				int midBin = hProj->GetXaxis()->FindBin(m1);
        double finalMean = 0, finalSigma = 0;
        
				// second pass will be useful for more precise analysis for mean and sigma
				PerformDynamicGaussianFit(hProj, midBin, m1, s1, finalMean, finalSigma);

        // Build final fit function
        TF1* finalFit = new TF1(Form("finalFit_%zu", i), "gaus", finalMean - 3*finalSigma, finalMean + 3*finalSigma);
        
				int cBin = hProj->GetXaxis()->FindBin(finalMean);  
				double avgH = (hProj->GetBinContent(cBin-1) + hProj->GetBinContent(cBin) + hProj->GetBinContent(cBin+1)) / 3.0;

        finalFit->SetParameters(avgH, finalMean, finalSigma);
        finalFit->SetLineColor(kRed);
        finalFit->SetLineWidth(2);

        // Calculate deviations and MAD within ±3σ
        int leftLimit  = hProj->GetXaxis()->FindBin(finalMean - 3*finalSigma);
        int rightLimit = hProj->GetXaxis()->FindBin(finalMean + 3*finalSigma);
        
				std::vector<double> devsL, devsR;
        
				for (int b = cBin-1; b >= leftLimit; --b) {
            double x = hProj->GetXaxis()->GetBinCenter(b);
            devsL.push_back(std::fabs(hProj->GetBinContent(b) - finalFit->Eval(x)));
        }
        
				for (int b = cBin+1; b <= rightLimit; ++b) {
            double x = hProj->GetXaxis()->GetBinCenter(b);
            devsR.push_back(std::fabs(hProj->GetBinContent(b) - finalFit->Eval(x)));
        }

        double madL = ComputeMAD(devsL);
        double madR = ComputeMAD(devsR);

        // Identify ignored/extrapolated bins based on deviation > MAD
        int outL_bin = -1;
        for (int b = cBin-1; b >= leftLimit; --b) {
            double x = hProj->GetXaxis()->GetBinCenter(b);
            double dev = std::fabs(hProj->GetBinContent(b) - finalFit->Eval(x));
            if (dev > madL) { outL_bin = b; break; }
        }
        
				int outR_bin = -1;
        for (int b = cBin+1; b <= rightLimit; ++b) {
            double x = hProj->GetXaxis()->GetBinCenter(b);
            double dev = std::fabs(hProj->GetBinContent(b) - finalFit->Eval(x));
            if (dev > madR) { outR_bin = b; break; }
        }
        
				double outL_x = (outL_bin>0 ? hProj->GetXaxis()->GetBinCenter(outL_bin) : -1);
        double outR_x = (outR_bin>0 ? hProj->GetXaxis()->GetBinCenter(outR_bin) : -1);

        // Log and record
        std::cout << "Extrapolated range :: L=" << outL_x << ", R=" << outR_x << ", Mean=" << finalMean << std::endl;
        
				// save in csv file
				csvOut << p1 << "," << p2 << "," << outL_x << "," << outR_x << "," << finalMean << "\n";

        // Draw histogram first
        TCanvas* c = new TCanvas(Form("c_%zu", i), Form("pT [%.1f,%.1f] GeV", p1, p2), 800, 600);
        hProj->Draw("E");
        
        // Fill outlier region in green
        if (outL_bin > 0 && outR_bin > outL_bin) {
            for (int b = outL_bin; b <= outR_bin; ++b) {
                double x1 = hProj->GetXaxis()->GetBinLowEdge(b);
                double x2 = hProj->GetXaxis()->GetBinLowEdge(b+1);
                double y  = hProj->GetBinContent(b);
                TBox* box = new TBox(x1, 0, x2, y);
                box->SetFillColorAlpha(kGreen, 0.3);
                box->SetLineColor(kGreen);
                box->Draw("same");
            }
        }
        
        // Draw fit function next
        finalFit->Draw("SAME");
        
        // Draw histogram last to be on top
        hProj->Draw("E SAME");

        // Annotate plot
        TLatex lbl; lbl.SetNDC(); lbl.SetTextSize(0.04);
        lbl.DrawLatex(0.15, 0.85, Form("Mean=%.4f", finalMean));
        lbl.DrawLatex(0.15, 0.80, Form("Sigma=%.4f", finalSigma));
        c->SaveAs(Form("pT_range_%.1f_%.1f.png", p1, p2));

        // Clean up
        delete finalFit;
        delete c;
        delete hProj;
    }

    file->Close();
    csvOut.close();
    std::cout << "Done, results in outliers.csv" << std::endl;
    gSystem->Exit(0);
}

