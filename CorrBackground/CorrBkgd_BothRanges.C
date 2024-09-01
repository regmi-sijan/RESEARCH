#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TF1.h>
#include <TLine.h>
#include <TSystem.h>
#include <TString.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <cstdio>
#include <thread> // For std::this_thread::sleep_for
#include <chrono> // For std::chrono::seconds
#include <tuple>
#include <algorithm>

// Constant parameters defined at the beginning
const char* INPUT_FILE_NAME = "new_input.root";
const char* OUTPUT_FILE_NAME = "b_output.root";
const char* HISTOGRAM_F_NAME = "DelR_pairpT_f";
const char* HISTOGRAM_B_NAME = "DelR_pairpT_b";
const char* INV_MASS_HIST_NAME = "pairInvMassPtdelR";
const char* INV_MASS_BKGD_HIST_NAME = "pairInvMassPtdelRBkgd";

// this parms are used to fit pol to find the peak location 
const double FIT_delR_MIN = 0.025;     
const double FIT_delR_MAX = 0.9;      

// this is exclusively designed to have initial parms for fitting (after excluding bins)
const double EXTRAP_delR_MIN = 0.02;  
const double EXTRAP_delR_MAX = 1.0;   

// working pairpT range for our analysis
const double pairpT_MIN = 1.8;  
const double pairpT_MAX = 18.0;


// Predefined ranges for pairpT projection
const double pairpT_RANGES[][2] = {
    {0.0, 2.0},
    {2.0, 2.5},
    {2.5, 3.0},
    {3.0, 3.5},
    {3.5, 4.0},
    {4.0, 4.5},
    {4.5, 5.0},
    {5.0, 5.5},
    {5.5, 6.5},
    {6.5, 7.5},
    {7.5, 10.0},
    {10.0, 18.0}
};

// size of the tuples defined earlier
const int NUM_RANGES = sizeof(pairpT_RANGES) / sizeof(pairpT_RANGES[0]);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this section is to modify the range of histogram so that it will be ready for fitting
std::tuple<double, double> ModifyEtaRange(double pT_min, double pT_max) {

		// we will make sure to return 2 points for truncating eta peak
		double a_1 = 0.35; 
		double a_2 = 0.65;

		std::vector<double> pT_1 = {2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 5.5, 6.0, 7.0, 8.0};
		std::vector<double> pT_2 = {2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5, 6.0, 7.0, 8.0, 9.0};

		std::vector<double> r1 = {0.35, 0.25, 0.24, 0.20, 0.17, 0.15, 0.13, 0.12, 0.11, 0.10, 0.09};
		std::vector<double> r2 = {0.65, 0.55, 0.5, 0.45, 0.40 , 0.35, 0.30, 0.26, 0.25, 0.24, 0.22};
			
		for (size_t i = 0; i < pT_1.size()-1; ++i) {
			if ((pT_min >= pT_1[i]) && (pT_min < pT_1[i+1])){
				a_1 = r1[i];
				a_2 = r2[i];				
			}
		}
		// define parms for larger range
		if (pT_min >= 8.0){
			a_1 = 0.09;
			a_2 = 0.22;
		}

	 // get minimum among a_1 and rangeMax and make it new a_1
	 //a_1 = std::min(a_1, rangeMax);
 return std::make_tuple(a_1, a_2);		
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::tuple<double, double> ModifyPi0Range(double pT_min, double pT_max) {

		// we will make sure to return 2 points for truncating eta peak
		double a_1 = 0.09; 
		double a_2 = 0.22;

		std::vector<double> pT_1 = {2.0,  2.5,  3.0,  3.5,  4.0,  4.5,  5.0, 5.5, 6.0, 7.0, 8.0};
		std::vector<double> pT_2 = {2.5,  3.0,  3.5,  4.0,  4.5,  5.0,  5.5, 6.0, 7.0, 8.0, 9.0};

		std::vector<double> r1 = {0.06, 0.05, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};
		std::vector<double> r2 = {0.18, 0.18, 0.14, 0.13, 0.13 , 0.13, 0.10, 0.09, 0.085, 0.08, 0.075};
			
		for (size_t i = 0; i < pT_1.size()-1; ++i) {
			if ((pT_min >= pT_1[i]) && (pT_min < pT_1[i+1])){
				a_1 = r1[i];
				a_2 = r2[i];				
			}
		}
		// define parms for larger range
		if (pT_min >= 8.0){
			a_1 = 0.03;
			a_2 = 0.07;
		}

	 // get minimum among a_1 and rangeMax and make it new a_1
	 //a_1 = std::min(a_1, rangeMax);
 return std::make_tuple(a_1, a_2);		
}	

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to extrapolate ignored bins using rangeMin and rangeMax
TH1D* extrapolateIgnoredBins(TH1D* proj_f, double a1, double a2, double pairpT_min, double pairpT_max) {

    // Find the bin numbers corresponding to rangeMin and rangeMax for pi0 peak
		int TotalBins = proj_f->GetNbinsX();	
		int FinalBin = proj_f->FindBin(EXTRAP_delR_MAX);
    
		// get bins corresponding to each a1, a2
		// a1 is the min for ignoring and a2 is the max
		int a1_bin = proj_f->FindBin(a1);	
		int a2_bin = proj_f->FindBin(a2);	
		
		// readjusting the lower bin where we will start ignoring things
    int newRangeMinBin = a1_bin-1;

		// adjust the min bin to start if we encounter the negative range
		if (newRangeMinBin < 0){newRangeMinBin = 1; }
		  
		// Prepare for fitting polynomials using the updated newRangeMinBin and bins from rangeMax to 0.8
    TF1* polFits[6];
    polFits[0] = new TF1("pol0", "pol0", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);
    polFits[1] = new TF1("pol1", "pol1", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);
    polFits[2] = new TF1("pol2", "pol2", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);
    polFits[3] = new TF1("pol3", "pol3", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);
    polFits[4] = new TF1("pol4", "pol4", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);
    polFits[5] = new TF1("pol5", "pol5", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);

    double bestChi2NDF_pol = 1e30;
    TF1* bestPolFit = nullptr;

		// fitting diff pol in given range after ignoring bins
    for (int i = 0; i < 6; ++i) {
        TH1F* fitHist = (TH1F*)proj_f->Clone();
        
				// Set bins in the ignored range to 0 for fitting
        for (int j = a1_bin; j <= a2_bin; ++j) {
					fitHist->SetBinContent(j, 0);
					fitHist->SetBinError(j, 1e10); // Set a very large error to exclude the bin from fitting
				}

				// this is the line that does the fitting
        fitHist->Fit(polFits[i], "RQ", "", proj_f->GetBinCenter(newRangeMinBin), EXTRAP_delR_MAX);
        
				// chi2 check to select which polynomial is better
				double chi2NDF = polFits[i]->GetChisquare() / polFits[i]->GetNDF();
        if (chi2NDF < bestChi2NDF_pol) {
          bestChi2NDF_pol = chi2NDF;
          bestPolFit = polFits[i];
        }

        delete fitHist;
    }

    // Calculate the standard deviation of residuals in the fitting range
		// the idea is to include the fluctuations in the extrapolated region based on the other parts of histogram
    std::vector<double> residuals;
    for (int i = a2_bin; i <= FinalBin; ++i) {
        double observedValue = proj_f->GetBinContent(i);
        double fittedValue = bestPolFit->Eval(proj_f->GetBinCenter(i));
        double residual = observedValue - fittedValue;
        residuals.push_back(residual);
    }

    double sum = 0.0;
    double sq_sum = 0.0;
    for (double r : residuals) {
        sum += r;
        sq_sum += r * r;
    }
    double mean_residual = sum / residuals.size();
    double stdev_residual = std::sqrt(sq_sum / residuals.size() - mean_residual * mean_residual);

    // Use the best polynomial fit to extrapolate the ignored region with added fluctuations
		for (int j = a1_bin; j <= a2_bin; ++j) {
				
			double x = proj_f->GetBinCenter(j);
			double extrapolatedValue = bestPolFit->Eval(x);

			// Add Gaussian fluctuations to the extrapolated value
			TRandom3 randGen(static_cast<unsigned int>(time(0)));
				
			// now adding fluctuations
			double fluctuatedValue = randGen.Gaus(extrapolatedValue, stdev_residual / std::sqrt(residuals.size()));
			//double fluctuatedValue =  extrapolatedValue; 

			if (fluctuatedValue < 0) {fluctuatedValue = 0;}

			// Update the modified histogram with the fluctuated value
			proj_f->SetBinContent(j, fluctuatedValue);
		}

    // Clean up dynamically allocated objects
    for (int i = 0; i < 6; ++i) {
        delete polFits[i];
    }

	return proj_f; // return the modified histogram
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Now we are trying to apply weights to get things done
TH2F* ApplyWeightInvMass(TH1D* delR_proj_m, TH1D* delR_proj_b, TH2F* mass_delR_b){

    // Ensure histograms are not null
    if (!delR_proj_m || !delR_proj_b || !mass_delR_b) {
        std::cerr << "Error: One or more input histograms are null." << std::endl;
        return nullptr;
    }

    // Clone mass_delR_b to create a weighted histogram
    TH2F* w_mass_delR_b = (TH2F*)mass_delR_b->Clone("w_mass_delR_b");

		// now divide modified signal by background (delR_proj_m by delR_proj_b) for weight
		delR_proj_m->Divide(delR_proj_b);

    // Loop over the x and y-axis in 2D inv mass hist
    int _bindelR = w_mass_delR_b->GetNbinsX();
    int _binpT = w_mass_delR_b->GetNbinsY();
		
		// x-axis is mass and y axis is piarpT
    for (int i = 1; i <= _bindelR; ++i) {
			
			// get bin content corresponding to delR axis
			double _delR_ = w_mass_delR_b->GetXaxis()->GetBinCenter(i);

			// get weight from 1D weight
			double _w = delR_proj_m->GetBinContent(delR_proj_m->FindBin(_delR_));
			
			for (int j = 1; j <= _binpT; ++j) { // getting into another axis of pairpT in 2D hist
				
				// get bincontent and delR from 2D Hist
				double _v = w_mass_delR_b->GetBinContent(i , j);
				
				// update the histogram
				w_mass_delR_b->SetBinContent(i, j, _v*_w);
			}
		}
  return w_mass_delR_b;  // Return the invariant mass histogram
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to modify the histograms (merged with main)
void CorrBkgd_BothRanges() {
    // Open an output file for modified histograms
    TFile* outputFile = new TFile(OUTPUT_FILE_NAME, "RECREATE");

		// Open the input ROOT file
    TFile* inputFile = TFile::Open(INPUT_FILE_NAME, "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open the input file." << std::endl;
        return;
    }

    // Retrieve the TH2F and TH3F histogram from the file
    TH2F* delR_pairpT_f = (TH2F*)inputFile->Get(HISTOGRAM_F_NAME);
    TH2F* delR_pairpT_b = (TH2F*)inputFile->Get(HISTOGRAM_B_NAME);

    TH3F* masshist_f = (TH3F*)inputFile->Get(INV_MASS_HIST_NAME);
    TH3F* masshist_b = (TH3F*)inputFile->Get(INV_MASS_BKGD_HIST_NAME);

    // Check if the histogram was successfully loaded
    if (!delR_pairpT_f || !delR_pairpT_b || !masshist_f || !masshist_b) {
        std::cerr << "Error: Could not find one or more required histogram in the file." << std::endl;
        inputFile->Close();
        return;
    }

		// Iterate over each predefined pairpT range
    for (int rangeIdx = 0; rangeIdx < NUM_RANGES; ++rangeIdx) {
				
        double pairpT_min = pairpT_RANGES[rangeIdx][0];
        double pairpT_max = pairpT_RANGES[rangeIdx][1];

        // Project the histograms for the current pairpT range
				int f_min_bin = delR_pairpT_f->GetYaxis()->FindBin(pairpT_min);
				int f_max_bin = delR_pairpT_f->GetYaxis()->FindBin(pairpT_max);
				
				int b_min_bin = delR_pairpT_b->GetYaxis()->FindBin(pairpT_min);
				int b_max_bin = delR_pairpT_b->GetYaxis()->FindBin(pairpT_max);

				TH1D* delR_proj_f = delR_pairpT_f->ProjectionX(Form("delR_proj_f_%d", rangeIdx), f_min_bin, f_max_bin);
        TH1D* delR_proj_b = delR_pairpT_b->ProjectionX(Form("delR_proj_b_%d", rangeIdx), b_min_bin, b_max_bin);

        // Extrapolate ignored bins
				// First we will extrapolate eta-peak region and then pi0 peak region
	
				// remove pi0 peak
				auto [min_range_pi0, max_range_pi0] = ModifyPi0Range(pairpT_min, pairpT_max);
				delR_proj_f = extrapolateIgnoredBins(delR_proj_f, min_range_pi0, max_range_pi0, pairpT_min, pairpT_max);			
				
				// removing eta-peak
				// this will work for pT higher than 2.5
				if (pairpT_min >= 3.0) {

					auto [min_range_eta, max_range_eta] = ModifyEtaRange(pairpT_min, pairpT_max);
					delR_proj_f = extrapolateIgnoredBins(delR_proj_f, min_range_eta, max_range_eta, pairpT_min, pairpT_max);
				}
			
				// Now we will use weight to the invariant mass histogram and dynamically save all the histogram required
				
				// first let us project the mass 3D hist to the required pairpT range
				masshist_f->GetYaxis()->SetRangeUser(pairpT_min, pairpT_max);
				masshist_b->GetYaxis()->SetRangeUser(pairpT_min, pairpT_max);

				// Project onto the x and z axes
				TH2F* mass_delR_f = (TH2F*)masshist_f->Project3D("xz");
				TH2F* mass_delR_b = (TH2F*)masshist_b->Project3D("xz");
							
				// Update the weight in inv mass histogram (2D) and return modified 2D Hist
				TH2F* mass_delR_mod = ApplyWeightInvMass(delR_proj_f, delR_proj_b, mass_delR_b);

				// Project the 2D histogram onto the x-axis with the dynamic name (for both foreground and background)
				TH1* _fg = mass_delR_f->ProjectionY(Form("f_%.0f_%.0f", 10*pairpT_min, 10*pairpT_max));
				TH1* _bg = mass_delR_mod->ProjectionY(Form("b_%.0f_%.0f", 10*pairpT_min, 10*pairpT_max));

				// Set the title for both histogram as same
				_fg->SetTitle("Pair Mass Histogram");
				_bg->SetTitle("Pair Mass Histogram");

				_fg->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
				_fg->SetYTitle("Counts");

				_bg->SetXTitle("M_{#gamma#gamma} [GeV/c^{2}]");
				_bg->SetYTitle("Counts");

				outputFile->cd(); // preparing output file to write histogram as requirement
				_fg->Write(); // write the foreground in the output file
				_bg->Write(); // write the background in the output file		
	
		} // this marks the completeness of all desired range of pairpT

    outputFile->Close();
    inputFile->Close();

		std::cout << "All Done." << std::endl;
		gSystem->Exit(0);
}


