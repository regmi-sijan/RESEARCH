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
#include <TSpectrum.h>
#include <iostream>
#include "TH1.h"
#include "TMath.h"
#include <memory>
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <limits>
#include "TH1D.h"
#include <cstddef>
#include <TFitResultPtr.h>
#include <TROOT.h>

// Constant parameters defined at the beginning
const char* INPUT_FILE_NAME = "new_run24_pp_combined.root";
const char* REF_FILE_NAME = "tt7_pitoo_dR_reco2.root";
const char* OUTPUT_FILE_NAME = "test_run24_output.root";
const char* HISTOGRAM_F_NAME = "DelR_pairpT_f";
const char* HISTOGRAM_B_NAME = "DelR_pairpT_b";
const char* INV_MASS_HIST_NAME = "pairInvMassPtdelR";
const char* INV_MASS_BKGD_HIST_NAME = "pairInvMassPtdelRBkgd";
 

// following range is used for both identifying peak (caused exclusively by pi0) and also for extrapolating (after excluding few/many bins)
const double RANGE_delR_MIN = 0.03;  
const double RANGE_delR_MAX = 1.0;

// range to scale bg to match fg
const double MASS_RANGE_MIN = 1.2;  
const double MASS_RANGE_MAX = 1.8;

// min pairpT to start to apply weighitng method
const double MIN_PairpT = 2.0;  


// Predefined pT ranges for projection
const double pairpT_RANGES[][2] = {
  {2.0, 3.0},
  {3.0, 4.0},
  {4.0, 5.0},
  {5.0, 7.0},
  {7.0, 9.0},
  {9.0, 11.0},
  {11.0, 15.0}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is also the helper function that computes the root mean square error
// we use RMSE for simple calculation (we can modify that later depending on our requirement)
// Compute the root–mean–square error between a TGraphErrors and a fit TF1.
// this is helper function to extrapolate function
double calculateRMSE(const TH1D* hist, const TF1* fitFunc, 
		     int ignoreStartBin, int ignoreEndBin, 
		     int fitStartBin, int fitEndBin)
{
  double sum2 = 0.0;
  int count = 0;

  // just incase if we have null pointer
  if (!hist || !fitFunc) {
    std::cerr << "Error: null pointer passed to calculateRMSE." << std::endl;
    return 0.0;
  }

  // Loop over the fitting region
  for (int bin = fitStartBin; bin <= fitEndBin; ++bin) {
        
    // Skip ignored bins
    if (bin >= ignoreStartBin && bin <= ignoreEndBin) continue;

    double x = hist->GetBinCenter(bin);
    double y = hist->GetBinContent(bin);

    double expected = fitFunc->Eval(x);
    double resid = y - expected;
    sum2 += resid * resid;
    ++count;
  }

  if (count == 0) {
    std::cerr << "Warning: No bins left after ignoring specified bins." << std::endl;
    return 0.0;
  }

  //return std::sqrt(sum2 / count);
  return sum2;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper function to find maxima, minima and bin closet to certain value 
// first is to pass vector and then second, third will start-end index and final one is the threshold
// this is helper function to FindRanges function
int FindClosestBinToThreshold(const std::vector<double>& _value,
			      int start_index, int end_index, double threshold) 
{
  auto iter = std::min_element(_value.begin() + start_index, _value.begin() + end_index,
			       [threshold](double a, double b) {
				 return std::abs(a - threshold) < std::abs(b - threshold);
			       }
			       );

  // return -1 if we do not find any match
  if (iter == _value.begin() + end_index) {return -1;}
		
  // else get the valid value
  return std::distance(_value.begin(), iter);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to analyze the histogram based on pairpT_min and threshold
// this FindRanges is also semi-helper function and we are using this in extrapolating function  
std::vector<double> FindRanges(TH1D* ref_hist_proj, double pairpT) {
		
  // checking if histogram is valid
  if(!ref_hist_proj || ref_hist_proj == nullptr){std::cout << "No valid projected histogram" << std::endl;}
    
  std::vector<double> result(4, -1.0);  // [a0, a1, a2, a3] , four points for pi0 and eta peak

  // determine ref_delR based on pairpT
  double ref_delR = -100.0;
  if (pairpT <= 3.0) ref_delR = 0.2;
  else if (pairpT <= 3.5) ref_delR = 0.16;
  else if (pairpT <= 4.0) ref_delR = 0.15;
  else if (pairpT <= 4.5) ref_delR = 0.13;
  else if (pairpT <= 5.0) ref_delR = 0.12;
  else if (pairpT <= 5.5) ref_delR = 0.11;
	else if (pairpT <= 7.0) ref_delR = 0.10;
  else if (pairpT <= 8.0) ref_delR = 0.095;
  else ref_delR = 0.09;

  // Prepare containers to save data/info
  std::vector<double> _delR, _value;

  // Get number of bins from histogram
  int NBins = ref_hist_proj->GetNbinsX();

  // Safety check
  if (!ref_hist_proj || NBins == 0) {
    std::cerr << "Histogram is null or has no bins!" << std::endl;
  }

  // Reserve memory
  _delR.reserve(NBins);
  _value.reserve(NBins);

  // Fill vectors with bin center and content using emplace_back
  for (int i = 1; i <= NBins; ++i) {
    _delR.emplace_back(ref_hist_proj->GetXaxis()->GetBinCenter(i));
    _value.emplace_back(ref_hist_proj->GetBinContent(i));
  }

  // find the global peak
  int pi0_peak_bin = FindClosestBinToThreshold(_value, 0, NBins-2, 1e15); // set very large to find peak

  // Find the closest value (in bin-index) to ref_delR in _delR
  int delR_index = FindClosestBinToThreshold(_delR, 0, _delR.size(), ref_delR); 
	
  if (delR_index < 0 || !delR_index){std::cout << "Something is wrong in delR" << std::endl;} 

  // find location of eta-meson peak
  // the idea to find pi0 and eta peak are same except the range we want to find
  int etaloc_index = FindClosestBinToThreshold(_value, delR_index, NBins-2, 1e15); // set very large value to find peak

  // finding to minima
  int minima_index = FindClosestBinToThreshold(_value, pi0_peak_bin, etaloc_index, 0.0);
  double minima_value = _value[minima_index];

  // setting up threshold
  double threshold = 1.5 * minima_value;

  // for a0_bin
  int a0_bin = FindClosestBinToThreshold(_value, 0, pi0_peak_bin, threshold);

  // for a1_bin
  int a1_bin = FindClosestBinToThreshold(_value, pi0_peak_bin, minima_index, threshold);

  // for a2_bin
  int a2_bin = FindClosestBinToThreshold(_value, minima_index, etaloc_index, threshold);

  // for a3 and a3_bin
  int a3_bin = FindClosestBinToThreshold(_value, etaloc_index, NBins-2, threshold);
 
  // cross check things
  if (a0_bin < 0 || a1_bin < 0 || a2_bin < 0 || a3_bin < 0){
    std::cout << "One of the a's parm is less than zero. please have a look" << std::endl;
    std::cout << "Bins of a0, a1, a2, a3 = " << a0_bin << " , " << a1_bin << " , " << a2_bin << " , " << a3_bin << std::endl;
  }
	
  // looking for a1 and a2 bin
  if (a1_bin > a2_bin){
    std::cout << "a2 is larger than a3, so double check this." << std::endl;
    a1_bin = (a1_bin+a2_bin)/2 - 2;
    a2_bin = a1_bin + 4;
  }

  ////////////
  // --- Save to result vector
  result[0] = _delR[a0_bin];
  result[1] = _delR[a1_bin];
  result[2] = _delR[a2_bin];
  result[3] = _delR[a3_bin];
	
  //std::cout << _delR[a0_bin] << " , " <<  _delR[a1_bin] << " , " <<  _delR[a2_bin] << " , " <<  _delR[a3_bin] << std::endl;

  return result;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TH1D* extrapolateIgnoredBins(TH1D* proj_f, const std::vector<double>& ranges, bool etamode) {
    
  // cross-check the size of vector "ranges"
  if (ranges.size() != 4) {
    std::cerr << "Error: ranges vector must contain exactly 4 values!" << std::endl;
  }
		
  // getting all values for pi0 and eta meson
  double a0 = ranges[0];
  double a1 = ranges[1];
  double a2 = ranges[2];
  double a3 = ranges[3];

  int TotalBins = proj_f->GetNbinsX();
  int FinalBin = proj_f->FindBin(RANGE_delR_MAX);
		
  // getting bin-numbers for fitting and ignoring for peak(s)
  // while fitting we will set to bins (we ignore) to zero with very high error and use all fitting range to fit
  int ignoreStartBin, ignoreEndBin, fitStartBin, fitEndBin;

  if (etamode) {
    ignoreStartBin = proj_f->FindBin(a2);
    ignoreEndBin = proj_f->FindBin(a3);
    
    fitStartBin = ignoreStartBin - 5; // No fitting before a1 in eta mode
    fitEndBin = ignoreEndBin + 5; // we can adjust that to few bins before if we need in eta mode
    
		//fitEndBin = std::min(ignoreEndBin + 5, FinalBin-5); // we can adjust that to few bins before if we need in eta mode
    //fitStartBin = std::max(proj_f->FindBin(a1)+5, ignoreStartBin-3); // No fitting before a1 in eta mode
  
	} 
  else {
    ignoreStartBin = proj_f->FindBin(a0);
    ignoreEndBin = proj_f->FindBin(a1);
    
    fitStartBin = std::max(3, proj_f->FindBin(a0) - 5); // leave 5 bins before a0, if possible otherwise take only 2 at least
    fitEndBin = proj_f->FindBin(a2);
  }

  // starting point for fitting
  const int n_pol = 3;
  TF1* polFits[n_pol];
  double bestRMSE = 1e30;
  TF1* bestPolFit = nullptr;

  double fitMinX = proj_f->GetBinCenter(fitStartBin);
  double fitMaxX = proj_f->GetBinCenter(fitEndBin);
	
  // testing with multiple polynomials
  for (int i = 0; i < n_pol; ++i) {
    polFits[i] = new TF1(Form("pol%d", i), Form("pol%d", i), fitMinX, fitMaxX);
    TH1D* fitHist = (TH1D*)proj_f->Clone();

    for (int j = 1; j <= TotalBins; ++j) {
      if (j >= ignoreStartBin && j <= ignoreEndBin) {
				fitHist->SetBinContent(j, 0);
				fitHist->SetBinError(j, 1e10);
      }
    }

/*
    fitHist->Fit(polFits[i], "RQ", "", fitMinX, fitMaxX);
    double rmse = calculateRMSE(fitHist, polFits[i], ignoreStartBin, ignoreEndBin, fitStartBin, fitEndBin);
		

		if (rmse < bestRMSE) {
      bestRMSE = rmse;
      bestPolFit = polFits[i];
    }
	
		TFitResultPtr r = fitHist->Fit(polFits[i], "SQN", "", fitMinX, fitMaxX);
		double chi2 = r->Chi2();
		double ndf = r->Ndf();
		double rmse = (ndf>0 ? chi2/ndf : 1e30);

*/

    // this is the line that does the fitting
    fitHist->Fit(polFits[i], "RQ", "", fitMinX, fitMaxX);
        
    // chi2 check to select which polynomial is better
    double chi2NDF = polFits[i]->GetChisquare() / polFits[i]->GetNDF();
    if (chi2NDF < bestRMSE) {
      bestRMSE = chi2NDF;
      bestPolFit = polFits[i];
    }


    delete fitHist;
  }

  // this is to add random fluctuations to make extrapolation look more natural
  // this is to take noise outside of fitting range so that we do not have any artificially good looking extrapolation
  // normally there is not high fluctualtion effect in our data due to this
  










TRandom3 randGen(0); // 0 = ROOT uses automatically the best random seed
  std::vector<double> residuals;
	
  for (int i = ignoreEndBin + 5; i < FinalBin; ++i) {
    double observed = proj_f->GetBinContent(i);
    double expected = bestPolFit->Eval(proj_f->GetBinCenter(i));
    residuals.push_back(observed - expected);
  }

  double mean = 0, sq_sum = 0;
  for (double r : residuals) {
    mean += r;
    sq_sum += r * r;
  }
  mean /= residuals.size();
  double stdev = std::sqrt(sq_sum / residuals.size() - mean * mean);

  for (int j = ignoreStartBin; j <= ignoreEndBin; ++j) {
    double x = proj_f->GetBinCenter(j);
    double value = bestPolFit->Eval(x);
    double fluctuated = randGen.Gaus(value, stdev / std::sqrt(residuals.size()));
    proj_f->SetBinContent(j, fluctuated < 0 ? 0 : fluctuated);
  }

  for (int i = 0; i < n_pol; ++i) delete polFits[i];
  return proj_f;
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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to modify the histograms (merged with main)
void HistBothRangeV1() {
	
	gROOT->SetBatch(kTRUE);

  // Open an output file for modified histograms
  TFile* outputFile = new TFile(OUTPUT_FILE_NAME, "RECREATE");

  // Open the input ROOT file
  TFile* inputFile = TFile::Open(INPUT_FILE_NAME, "READ");

  TFile *ref_file = TFile::Open(REF_FILE_NAME);
    
  if (!ref_file || ref_file->IsZombie() || !inputFile || inputFile->IsZombie()) 
    {
      std::cerr << "Error opening file!" << std::endl;
      return;
    }

  // Get the TH3F histogram from the file
  TH3F *ref_hist1 = dynamic_cast<TH3F*>(ref_file->Get("ptp_pt1_dang"));
  if (!ref_hist1) 
    {
      std::cerr << "Histogram 'ptp_pt1_dang' not found!" << std::endl;
      ref_file->Close();
      return;
    }

  // Project the 3D histogram onto the x and z axes
  // original pT axis changed to Y-axis and delR changed to X-axis
  TH2F *ref_hist = (TH2F*)ref_hist1->Project3D("xz");
  //ref_hist->SetMinimum(20.0);
	
  // Retrieve the TH2F and TH3F histogram from the file
  TH3F* masshist_f = (TH3F*)inputFile->Get(INV_MASS_HIST_NAME);
  TH3F* masshist_b = (TH3F*)inputFile->Get(INV_MASS_BKGD_HIST_NAME);
	  
  // Retrieve the TH2F and TH3F histogram from the file
  TH2F* delR_pairpT_f = (TH2F*)inputFile->Get(HISTOGRAM_F_NAME);
  TH2F* delR_pairpT_b = (TH2F*)inputFile->Get(HISTOGRAM_B_NAME);

  // Check if the histogram was successfully loaded
  if (!delR_pairpT_f || !delR_pairpT_b || !masshist_f || !masshist_b) 
    {
      std::cerr << "Error: Could not find one or more required histogram in the file." << std::endl;
      inputFile->Close();
      return;
    }

  // Iterate over each bin in pairpT
  int numBins = delR_pairpT_f->GetYaxis()->GetNbins(); // get total number of bins


  // Loop over predefined pT ranges
  for (auto &range : pairpT_RANGES) {

    // getting min and max pairpT range 
    double range_min = range[0];
    double range_max = range[1];

		//std::cout << "pairpT range is = " << range_min << " , " << range_max << std::endl;
		
    // starting and ending bin for this pairpT range
    int MainStartBin = delR_pairpT_f->GetYaxis()->FindBin(range_min);
    int MainEndBin  = delR_pairpT_f->GetYaxis()->FindBin(range_max); // do not subtract 1 here we will track this later
    
    if (MainEndBin < MainStartBin) {
      std::cout << "Our MainEndBin is smaller than MainStartBin" << std::endl;
      std::exit(1);
    }
		
    // for fg and bg
    TH1D* _fg = nullptr;
    TH1D* _bg = nullptr;

    // we will start to slice each pairpT range into variable sized multiple pieces with enough number of statistics
    int bin = MainStartBin; // bin will be useful in while loop
    int binStart = MainStartBin; // binStart is to know which bin is essential for start within while loop

    while (bin < MainEndBin) { // always have less than in this while loop
      
      double totalCounts = 0.0;
      int binStart = bin;
      int binEnd = binStart;

      // accumulate until threshold or end of this pT subrange
      // we will keep this loop until we get enough number of stats
      while (totalCounts < 50000 && binEnd < MainEndBin) { // also have less than (only) in this loop as well
				++binEnd;
        totalCounts = delR_pairpT_f->Integral(0, -1, binStart, binEnd);
      }
			
      // now we will check counts for remaining bins within the range
      double checkCounts = delR_pairpT_f->Integral(0, -1, binEnd, MainEndBin);
			
      // if we have less statistics in remaining bins then we merge them with current range
      if(checkCounts < 50000){binEnd = MainEndBin;}

      // looking for min and maximum pairpT range for this window 
      double temp_pairpT_min = delR_pairpT_f->GetYaxis()->GetBinCenter(binStart);
      double temp_pairpT_max = delR_pairpT_f->GetYaxis()->GetBinCenter(binEnd);

      // Find corresponding bins in histograms for the range
      int f_min_bin = binStart;
      int f_max_bin = binEnd;

      int b_min_bin = delR_pairpT_b->GetYaxis()->FindBin(temp_pairpT_min);
      int b_max_bin = delR_pairpT_b->GetYaxis()->FindBin(temp_pairpT_max);

      int ref_min_bin = ref_hist->GetYaxis()->FindBin(temp_pairpT_min);
      int ref_max_bin = ref_hist->GetYaxis()->FindBin(temp_pairpT_max);

      // Project histograms for the current range
      TH1D* delR_proj_f = delR_pairpT_f->ProjectionX(Form("delR_proj_f_%d", binStart), f_min_bin, f_max_bin - 1);
      TH1D* delR_proj_b = delR_pairpT_b->ProjectionX(Form("delR_proj_b_%d", binStart), b_min_bin, b_max_bin - 1);
			
      // to have at least one bin in ref_hist
      if(abs(ref_max_bin - ref_min_bin) < 2){ref_max_bin = ref_min_bin +2;}

      TH1D* ref_hist_proj = ref_hist->ProjectionX(Form("ref_hist_%d", binStart), ref_min_bin, ref_max_bin - 1);
      //ref_hist_proj->SetMinimum(10);
      ref_hist_proj->Smooth(5);

      // Removing peaks
      TH1D* delR_proj_m = extrapolateIgnoredBins(delR_proj_f, FindRanges(ref_hist_proj, temp_pairpT_min), true); // Remove eta peak
      //TH1D* delR_proj_m = extrapolateIgnoredBins(delR_proj_m1, FindRanges(ref_hist_proj, temp_pairpT_min), false); // Remove pi0 peak
      
			// Process mass histograms, this seems to be working
      masshist_f->GetYaxis()->SetRangeUser(temp_pairpT_min, temp_pairpT_max);
      masshist_b->GetYaxis()->SetRangeUser(temp_pairpT_min, temp_pairpT_max);
				
      // projecting in x-z plane
      TH2F* mass_delR_f = (TH2F*)masshist_f->Project3D("xz");
      TH2F* mass_delR_b = (TH2F*)masshist_b->Project3D("xz");
			
      // removing projection in y-axis, setting that free
      masshist_f->GetYaxis()->SetRange();
      masshist_b->GetYaxis()->SetRange();

      // Apply weights, this will modify mass_delR_b to mass_delR_mod but we are still left with scaling
      TH2F* mass_delR_mod = ApplyWeightInvMass(delR_proj_m, delR_proj_b, mass_delR_b);

      /*
      // --- Add this block to display it ---
      TCanvas* cRefProj = new TCanvas("cRefProj", "Ref ΔR Projection", 800, 600);
      mass_delR_f->Draw();      // default histogram draw
      cRefProj->Update();

      // Option A: pause so you can view it interactively
      std::this_thread::sleep_for(std::chrono::seconds(1));      
      */ 

      TH1* _fg_temp = mass_delR_f->ProjectionY(Form("f_%.0f", 10*temp_pairpT_min));
      TH1* _bg_temp = mass_delR_mod->ProjectionY(Form("b_%.0f", 10*temp_pairpT_min));

      // Scale background to foreground
      int min_th = _fg_temp->FindBin(MASS_RANGE_MIN);
      int max_th = _fg_temp->FindBin(MASS_RANGE_MAX);

      double _fg_sum = 0.0, _bg_sum = 0.0;
      for (int i = min_th; i <= max_th; ++i){
				_fg_sum += _fg_temp->GetBinContent(i);
				_bg_sum += _bg_temp->GetBinContent(i);
			}

      double scale_factor = (_bg_sum > 0) ? (_fg_sum / _bg_sum) : 1.0;
      if (_bg_sum > 0) _bg_temp->Scale(scale_factor);


      // now update to the main fg and bg    
      if (!_fg) {
				
				// First slice: clone both histograms so _fg/_bg inherit the correct binning & axes
				_fg = (TH1D*)_fg_temp->Clone(Form("f_%.0f", 10*temp_pairpT_min));
				//_fg->SetDirectory(nullptr);    // detach from any file so ROOT won’t auto-own it
				//_fg->Sumw2();                  // enable proper error propagation

				_bg = (TH1D*)_bg_temp->Clone(Form("b_%.0f", 10*temp_pairpT_min));
				//_bg->SetDirectory(nullptr);
				//_bg->Sumw2();
      }
			
      else {
				// Subsequent slices: add each into the accumulator
				_fg->Add(_fg_temp);
				_bg->Add(_bg_temp);
      }

      // Clean up the temporaries each time
      //delete _fg_temp;
      //delete _bg_temp;

      // Move binStart to the next unprocessed bin
      bin = binEnd;
		
    }

    // write the output file
    // Save histograms
    outputFile->cd();
    _fg->Write();
    _bg->Write();
  }
  
  // closing output file and 
  outputFile->Close();
  inputFile->Close();

  std::cout << "All Done." << std::endl;
  gSystem->Exit(0);

}



