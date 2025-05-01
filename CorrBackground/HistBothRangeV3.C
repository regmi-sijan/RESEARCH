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
#include <TSpectrum.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TSpline.h>
#include <TGraph.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <string>
#include <cstdio>
#include <thread>
#include <chrono>
#include <tuple>
#include <algorithm>
#include <memory>
#include <limits>
#include <cstddef>
#include "TH1.h"
#include "TMath.h"
#include <TFitResultPtr.h>


// Constant parameters defined at the beginning
const char* INPUT_FILE_NAME = "combine_run24_pp_photrig_evtmix.root";
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
  {9.0, 12.0},
  {12.0, 18.0}
};

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

// this region need a lot more tuning and manual interaction to make it work for different pT ranges
// following are the parameters required for tuning
double _alpha = 1.0;
double _exponent = 0.2;

// Original minmod function
double minmod(double a, double b) {
  return (a * b <= 0) ? 0.0 : ((std::abs(a) < std::abs(b)) ? a : b);
}

// this is where we are supposed to make changes/update
double adaptiveSlope(double slope, double secant, double x, double x0, double x1) {
  double t = (x - x0) / (x1 - x0);
  t = std::clamp(t, 0.0, 1.0);

  // Exponential decay: more weight at start, tighter later
  double weight = _exponent * std::exp(-_alpha * t);  // max 0.5 at t=0, near 0 at t=1

  double limited = minmod(slope, secant);
  return weight * slope + (1.0 - weight) * limited;
}


// Main extrapolation function
TH1D* extrapolateIgnoredBins(TH1D* proj_f, const std::vector<double>& ranges, bool etamode) {
  if (!proj_f || ranges.size() != 4) {
    std::cerr << "Invalid histogram or range vector." << std::endl;
    return nullptr;
  }

  double a0 = ranges[0], a1 = ranges[1], a2 = ranges[2], a3 = ranges[3];
  int ignoreStartBin, ignoreEndBin;

  if (etamode) {
    ignoreStartBin = proj_f->FindBin(a2);
    ignoreEndBin   = proj_f->FindBin(a3) + 10;
  } else {
    ignoreStartBin = proj_f->FindBin(a0);
    ignoreEndBin   = proj_f->FindBin(a1);
  }

  int totalBins = proj_f->GetNbinsX();
  ignoreStartBin = std::clamp(ignoreStartBin, 1, totalBins);
  ignoreEndBin   = std::clamp(ignoreEndBin, 1, totalBins);

  double x0 = proj_f->GetBinCenter(ignoreStartBin);
  double x1 = proj_f->GetBinCenter(ignoreEndBin);
  if (x1 < x0) {
    std::swap(x0, x1);
    std::swap(ignoreStartBin, ignoreEndBin);
  }

  if (std::abs(x1 - x0) < 1e-3) {
    std::cerr << "Interpolation region too small." << std::endl;
    return proj_f;
  }

  // Evaluate values and naive secant slope
  double y0 = proj_f->Interpolate(x0);
  double y1 = proj_f->Interpolate(x1);
  double h = x1 - x0;
  double secant = (y1 - y0) / h;

  // Estimate slopes using centered differences
  double dx = h / 20.0;
  double dy0 = (proj_f->Interpolate(x0 + dx) - proj_f->Interpolate(x0 - dx)) / (2 * dx);
  double dy1 = (proj_f->Interpolate(x1 + dx) - proj_f->Interpolate(x1 - dx)) / (2 * dx);

  // Apply adaptive limiter with progressive constraint
  double m0 = adaptiveSlope(dy0, secant, x0, x0, x1);
  double m1 = adaptiveSlope(dy1, secant, x1, x0, x1);

  // Monotonic Hermite interpolation
  auto pchipInterp = [&](double x) {
    double t = (x - x0) / h;
    double t2 = t * t;
    double t3 = t2 * t;

    double h00 = 2*t3 - 3*t2 + 1;
    double h10 = t3 - 2*t2 + t;
    double h01 = -2*t3 + 3*t2;
    double h11 = t3 - t2;

    return h00 * y0 + h * h10 * m0 + h01 * y1 + h * h11 * m1;
  };

  // Fill ignored bins using interpolated values
  for (int bin = ignoreStartBin; bin <= ignoreEndBin; ++bin) {
    double x = proj_f->GetBinCenter(bin);
    double y = pchipInterp(x);
    proj_f->SetBinContent(bin, y < 0 ? 0 : y);
  }

  return proj_f;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  TH1D* weight_hist = (TH1D*)delR_proj_m->Clone("weight_hist");
  weight_hist->Divide(delR_proj_b);

  // Loop over the x and y-axis in 2D inv mass hist
  int _bindelR = w_mass_delR_b->GetNbinsX();
  int _binpT = w_mass_delR_b->GetNbinsY();
		
  // x-axis is mass and y axis is piarpT
  for (int i = 1; i <= _bindelR; ++i) {
			
    // get bin content corresponding to delR axis
    double _delR_ = w_mass_delR_b->GetXaxis()->GetBinCenter(i);

    // get weight from 1D weight
    double _w = weight_hist->GetBinContent(weight_hist->FindBin(_delR_));
			
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
void HistBothRangeV3() {
	
  //gROOT->SetBatch(kTRUE);

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
    double pT_min = range[0];
    double pT_max = range[1];

    // Find corresponding bins in histograms for the range
    int f_min_bin = delR_pairpT_f->GetYaxis()->FindBin(pT_min);
    int f_max_bin = delR_pairpT_f->GetYaxis()->FindBin(pT_max);

    int b_min_bin = delR_pairpT_b->GetYaxis()->FindBin(pT_min);
    int b_max_bin = delR_pairpT_b->GetYaxis()->FindBin(pT_max);

    int ref_min_bin = ref_hist->GetYaxis()->FindBin(pT_min);
    int ref_max_bin = ref_hist->GetYaxis()->FindBin(pT_max);

    // Project histograms for the current range
    TH1D* delR_proj_f = delR_pairpT_f->ProjectionX(Form("delR_proj_f_%f", 10*pT_min), f_min_bin, f_max_bin - 1);
    TH1D* delR_proj_b = delR_pairpT_b->ProjectionX(Form("delR_proj_b_%f", 10*pT_min), b_min_bin, b_max_bin - 1);
			

    TH1D* ref_hist_proj = ref_hist->ProjectionX(Form("ref_hist_%f", 10*pT_min), ref_min_bin, ref_max_bin - 1);
    //ref_hist_proj->SetMinimum(10);
    ref_hist_proj->Smooth(5);

    // Removing peaks, currently dedicated to eta-meson only
    TH1D* delR_proj_m = extrapolateIgnoredBins(delR_proj_f, FindRanges(ref_hist_proj, pT_min), true); // Remove eta peak
    //TH1D* delR_proj_m = extrapolateIgnoredBins(delR_proj_m1, FindRanges(ref_hist_proj, pT_min), false); // Remove pi0 peak
 

    // Process mass histograms, this seems to be working
    // may be work on SetRangeUser() stuff
    masshist_f->GetYaxis()->SetRangeUser(pT_min, pT_max);
    masshist_b->GetYaxis()->SetRangeUser(pT_min, pT_max);
				
    // projecting in x-z plane
    TH2F* mass_delR_f = (TH2F*)masshist_f->Project3D("xz");
    TH2F* mass_delR_b = (TH2F*)masshist_b->Project3D("xz");
			
    //std::cout << "After set-range = " << mass_delR_f->Integral() << " , " << mass_delR_b->Integral() << std::endl; 
    
    // removing projection in y-axis, setting that free
    masshist_f->GetYaxis()->SetRange();
    masshist_b->GetYaxis()->SetRange();

    // Apply weights, this will modify mass_delR_b to mass_delR_mod but we are still left with scaling
    TH2F* mass_delR_mod = ApplyWeightInvMass(delR_proj_m, delR_proj_b, mass_delR_b);


    TH1* _fg = mass_delR_f->ProjectionY(Form("f_%.0f", 10*pT_min));
    TH1* _bg = mass_delR_mod->ProjectionY(Form("b_%.0f", 10*pT_min));

    // Scale background to foreground
    int min_th = _fg->FindBin(MASS_RANGE_MIN);
    int max_th = _fg->FindBin(MASS_RANGE_MAX);

    double _fg_sum = 0.0, _bg_sum = 0.0;
    for (int i = min_th; i <= max_th; ++i){
      _fg_sum += _fg->GetBinContent(i);
      _bg_sum += _bg->GetBinContent(i);
    }

    double scale_factor = (_bg_sum > 0) ? (_fg_sum / _bg_sum) : 1.0;
    if (_bg_sum > 0) _bg->Scale(scale_factor);


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



