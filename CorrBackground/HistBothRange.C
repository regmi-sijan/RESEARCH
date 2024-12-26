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

using namespace std;

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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to calculate RMSE from histogram and fitted function
double calculateRMSE(TH1* fitHist, const TF1* fitFunc, int a1_bin, int a2_bin, int f_min, int f_max) 
{
  double sumSquaredErrors = 0.0;
  int nPoints = 1;

  // Loop over each bin in the histogram
  for (int i = 1; i <= fitHist->GetNbinsX(); ++i) 
    {
      if (( i > a1_bin ) && (i < a2_bin)) continue; // ignore the bins that we will extrapolate
      if((i < f_min) || (i > f_max)) continue; // ignoring outside bins than for fitting
    
      double observed = fitHist->GetBinContent(i); // Observed value (bin content)
      double binCenter = fitHist->GetBinCenter(i); // Center of the bin
      double predicted = fitFunc->Eval(binCenter); // Predicted value from the fit function

      sumSquaredErrors += (std::pow((observed - predicted), 2))/binCenter; // Sum of squared errors
      nPoints += 1; // updating the total number of data points
    }

  // Calculate RMSE
  double rmse = std::sqrt(sumSquaredErrors / nPoints); // Return RMSE
  return rmse;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Function to analyze the histogram based on pairpT_min and threshold
std::vector<double> FindRanges(TH1D* ref_hist_proj, double pairpT) {

  // Vectors to store the results
  std::vector<double> result(5, -1.0);  // [a0, a1, a2, a3, a4]

  // determine the threshold to look for maxima
  // then we will be looking for minima and maxima
  double ref_delR = -100.0;
 
  if(pairpT <= 3.0){ref_delR = 2.0;}
  else if((pairpT > 3.0)&&(pairpT <= 4.5)){ref_delR = 1.4;}
  else if((pairpT > 4.5)&&(pairpT <= 6.0)){ref_delR = 1.2;}
  else if (pairpT > 6.0){ref_delR = 0.95;}

  int TotalBins = ref_hist_proj->GetNbinsX();
		 
  // Find the global peak (corresponding to pi0)
  ref_hist_proj->GetXaxis()->SetRangeUser(0.0, ref_delR);

  double pi0_peak = ref_hist_proj->GetMaximum();
  int pi0_peak_bin = ref_hist_proj->GetMaximumBin();

  // now lets find the peak corresponding to eta meson
  ref_hist_proj->GetXaxis()->SetRangeUser(ref_delR, 3.0);
  double eta_peak = ref_hist_proj->GetMaximum();
  int eta_peak_bin = ref_hist_proj->GetMaximumBin();

  // looking for peak minimum
  ref_hist_proj->GetXaxis()->SetRangeUser(pi0_peak, eta_peak);
  double min_peak = ref_hist_proj->GetMinimum();
  int min_peak_bin = ref_hist_proj->GetMinimumBin();


  // setting threshold to identify a1 to a4
  // this can be modified as required in future
  double pi0_thres_right = min_peak + 0.1 * (pi0_peak - min_peak);
  double eta_thres_left = min_peak + 0.1 * (eta_peak - min_peak);

  // initializating a0 to a4
  double a0 = -1; // a0 is starting of fitting of eta-meson peak		
  double a1 = -1;
  double a2 = -1;
  double a3 = -1;
  double a4 = -1;
		
  // Iterate from left to the peak to find where content is around pi0_thres_val
  // setting 1 percent of pi0 peak as lower threshold
  for (int bin = 1; bin <= pi0_peak_bin; ++bin) {
    if (ref_hist_proj->GetBinContent(bin) > (0.05 * pi0_peak)) {
      a1 = ref_hist_proj->GetBinCenter(bin);
      break;
    }
  }
		
  // to find a2
  // setting 1 percent of height of pi0 peak from the min value
  // "5" bins will be an offset
  for (int bin = pi0_peak_bin; bin <= min_peak_bin + 5; ++bin) {
    if (ref_hist_proj->GetBinContent(bin) <  pi0_thres_right) {
      a2 = ref_hist_proj->GetBinCenter(bin);
      break;
    }
  }
		
  // to find a3
  // setting 1 percent of height of eta peak from the min value
  // setting "5" as offset
  for (int bin = eta_peak_bin; bin > min_peak_bin - 5; --bin) {
    if (ref_hist_proj->GetBinContent(bin) <  eta_thres_left) {
      a3 = ref_hist_proj->GetBinCenter(bin);
      break;
    }
  }

  // to find a4
  for (int bin = eta_peak_bin; bin < TotalBins; ++bin) {
    if (ref_hist_proj->GetBinContent(bin) < (0.05 * eta_peak)) {
      a4 = ref_hist_proj->GetBinCenter(bin);
      break;
    }
  }

  if (a4 < 0){
    std::cout << "We do not get any a4 form histogram. Now we will have some manual adjustment to it." << std::endl;
    a4 = 0.55;
  } 

  if ((a1 < 0) || (a2 < 0) || (a3 < 0)) {
    std::cout << "Something is wrong with a3. No a3 found in hist. Please check it." << std::endl;
  }

  if (a2 > a3){
    //std::swap(a2, a3);
    std::cout << "a2 is larger than a3 which is abnormal situation." << std::endl;
  }
		
  // to identify a0
  //int a0_bin = (pi0_peak_bin + min_peak_bin )/2;
  //a0 = ref_hist_proj->GetBinCenter(a0_bin);
  //a0 = ref_hist_proj->GetBinCenter(pi0_peak_bin);
		
  a0 = a2 + (a3-a2)/3;
		
  // updating the output vector
  result[0] = a0;
  result[1] = a1;
  result[2] = a2;
  result[3] = a3;
  result[4] = a4;

  std::cout << "a0 to a4 = " << " , " << a0  << " , " << a1 << " , " << a2 << " , " << a3 << " , " << a4 << std::endl;

  // Reset the histogram range
  ref_hist_proj->GetXaxis()->SetRange();  // Reset to full range

  return result; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Function to extrapolate ignored bins using rangeMin and rangeMax
TH1D* extrapolateIgnoredBins(TH1D* proj_f, std::vector<double> ranges, bool etamode) {

  // Find the bin numbers corresponding to rangeMin and rangeMax for pi0 peak
  int TotalBins = proj_f->GetNbinsX();	
  int FinalBin = proj_f->FindBin(RANGE_delR_MAX);

  // Ensure that rangeMax has exactly 3 values
  // first value is the starting point, second value is the begining of ignored region and the third one is the end of ignored region
  if (ranges.size() != 5) 
    {
      std::cerr << "Error: rangeMax vector must contain exactly 3 values!" << std::endl;
      return nullptr;
    }

  // set the starting of the fitting range
  int newRangeMinBin = proj_f->FindBin(RANGE_delR_MIN);

  // Extract the three necessary values from rangeMax
  double a0 = ranges[0]; // we will update a0 only for eta meson as of now, a0 for pi0 is taken care by "RANGE_delR_MAX"
  double a1 = ranges[1];
  double a2 = ranges[2];
  double a3 = ranges[3];
  double a4 = ranges[4];

  // initialize the parameters
  int a0_bin = -1;
  int a1_bin = -1;
  int a2_bin = -1;

  // get those bin ranges depending on if we are running in eta mode or not	
  if(etamode)
    {
      a0_bin = proj_f->FindBin(a0);
      a1_bin = proj_f->FindBin(a3);
      a2_bin = proj_f->FindBin(a4);
      newRangeMinBin = a0_bin;
    }
  else
    {
      a1_bin = proj_f->FindBin(a1);
      a2_bin = proj_f->FindBin(a2);
			
      if (a1_bin  <= newRangeMinBin)
	{
	  a1_bin = newRangeMinBin + 1;
	}
    }		
		  
  // Prepare for fitting polynomials using the updated newRangeMinBin and bins from rangeMax to 0.8
  TF1* polFits[5];
  std::cout << "newRangeBinMin = " << proj_f->GetBinCenter(newRangeMinBin) << std::endl;
  polFits[0] = new TF1("pol0", "pol0", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);
  polFits[1] = new TF1("pol1", "pol1", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);
  polFits[2] = new TF1("pol2", "pol2", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);
  polFits[3] = new TF1("pol3", "pol3", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);
  polFits[4] = new TF1("pol4", "pol4", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);
  //polFits[5] = new TF1("pol5", "pol5", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);

  double bestRMSE = 1e30;
  TF1* bestPolFit = nullptr;

  // fitting diff pol in given range after ignoring bins
  for (int i = 0; i < 5; ++i) 
    {
      TH1D* fitHist = (TH1D*)proj_f->Clone();

      // Set bins in the ignored range to 0 for fitting
      for (int j = 1; j <= TotalBins; ++j) 
	{
					
	  if ((j >= a1_bin) && (j <= a2_bin))
	    {
	      fitHist->SetBinContent(j, 0);
	      fitHist->SetBinError(j, 1e10); // Set a very large error to exclude the bin from fitting
	    }
	}
				
      // this is the line that does the fitting
      fitHist->Fit(polFits[i], "RQ", "", proj_f->GetBinCenter(newRangeMinBin), RANGE_delR_MAX);

      int f_min = newRangeMinBin;
      int f_max = proj_f->FindBin(RANGE_delR_MAX);
        
      // rmse check to select which polynomial is better
      double rmse = calculateRMSE(fitHist, polFits[i], a1_bin, a2_bin, f_min, f_max);
      std::cout << rmse << std::endl;
				
      if (rmse < bestRMSE) 
	{
	  bestRMSE = rmse;
	  bestPolFit = polFits[i];
	}

      delete fitHist;
    }

  // Calculate the standard deviation of residuals in the fitting range
  // the idea is to include the fluctuations in the extrapolated region based on the other parts of histogram
  // "5" below is the offset value, can be changed if needed with any non-negative integer
  std::vector<double> residuals;
  for (int i = a2_bin+5; i < FinalBin; ++i) 
    {
      double observedValue = proj_f->GetBinContent(i);
      double fittedValue = bestPolFit->Eval(proj_f->GetBinCenter(i));
      double residual = observedValue - fittedValue;
      residuals.push_back(residual);
    }

  double sum = 0.0;
  double sq_sum = 0.0;
  for (double r : residuals) 
    {
      sum += r;
      sq_sum += r * r;
    }
  double mean_residual = sum / residuals.size();
  double stdev_residual = std::sqrt(sq_sum / residuals.size() - mean_residual * mean_residual);

  // Use the best polynomial fit to extrapolate the ignored region with added fluctuations
  for (int j = 1; j <= TotalBins; ++j) 
    {
      if ((j > a1_bin-1) && (j < a2_bin+1))
	{
		
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
    }
	
  // Clean up dynamically allocated objects
  for (int i = 0; i < 5; ++i) 
    {
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


// Function to modify the histograms (merged with main)
void HistBothRange() {

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
  TH2F *ref_hist = (TH2F*)ref_hist1->Project3D("xz");
  //ref_hist->SetMinimum(10.0);

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
  int binStart = 1; // Initialize starting bin

  while (binStart <= numBins) 
    {
      double pairpT_min = delR_pairpT_f->GetYaxis()->GetBinCenter(binStart);
      double pairpT_max = pairpT_min; // Initially set max edge to the starting bin (we will later update it)
      double totalCounts = 0.0;

      int binEnd = binStart;

      // Dynamically combine bins until the count threshold (40K) is met
      while (totalCounts < 40000 && binEnd <= numBins) 
	{
	  ++binEnd; // Expand the range and do further analysis
      
	  // new maximum range of pairpT
	  pairpT_max = delR_pairpT_f->GetYaxis()->GetBinCenter(binEnd);  
      
	  // getting total counts
	  totalCounts = delR_pairpT_f->Integral(binStart, binEnd, 0, -1);
	}

      // Find corresponding bins in histograms for the range
      int f_min_bin = delR_pairpT_f->GetYaxis()->FindBin(pairpT_min);
      int f_max_bin = delR_pairpT_f->GetYaxis()->FindBin(pairpT_max);

      int b_min_bin = delR_pairpT_b->GetYaxis()->FindBin(pairpT_min);
      int b_max_bin = delR_pairpT_b->GetYaxis()->FindBin(pairpT_max);

      int ref_min_bin = ref_hist->GetYaxis()->FindBin(pairpT_min);
      int ref_max_bin = ref_hist->GetYaxis()->FindBin(pairpT_max);

      // Validate reference range
      if ((ref_min_bin < 1 || ref_max_bin < 1) || ref_min_bin > ref_max_bin) continue;

      // Project histograms for the current range
      TH1D* delR_proj_f = delR_pairpT_f->ProjectionX(Form("delR_proj_f_%d", binStart), f_min_bin, f_max_bin - 1);
      TH1D* delR_proj_b = delR_pairpT_b->ProjectionX(Form("delR_proj_b_%d", binStart), b_min_bin, b_max_bin - 1);

      TH1D* ref_hist_proj = ref_hist->ProjectionX(Form("ref_hist_%d", binStart), ref_min_bin, ref_max_bin - 1);
      ref_hist_proj->SetMinimum(10);
      ref_hist_proj->Smooth(5);

			cout << "do we see this" << std::endl; 

      // Removing peaks
      TH1D* delR_proj_m1 = extrapolateIgnoredBins(delR_proj_f, FindRanges(ref_hist_proj, pairpT_min), true); // Remove eta peak
      TH1D* delR_proj_m = extrapolateIgnoredBins(delR_proj_m1, FindRanges(ref_hist_proj, pairpT_min), false); // Remove pi0 peak

      // Process mass histograms
      masshist_f->GetYaxis()->SetRangeUser(pairpT_min, pairpT_max);
      masshist_b->GetYaxis()->SetRangeUser(pairpT_min, pairpT_max);

      TH2F* mass_delR_f = (TH2F*)masshist_f->Project3D("xz");
      TH2F* mass_delR_b = (TH2F*)masshist_b->Project3D("xz");

      // Apply weights and project
      TH2F* mass_delR_mod = ApplyWeightInvMass(delR_proj_m, delR_proj_b, mass_delR_b);

      TH1* _fg = mass_delR_f->ProjectionY(Form("f_%.0f_%.0f", 10 * pairpT_min, 10 * pairpT_max));
      TH1* _bg = mass_delR_mod->ProjectionY(Form("b_%.0f_%.0f", 10 * pairpT_min, 10 * pairpT_max));

      // Scale background to foreground
      int min_th = _fg->FindBin(MASS_RANGE_MIN);
      int max_th = _fg->FindBin(MASS_RANGE_MAX);

      double _fg_sum = 0.0, _bg_sum = 0.0;
      for (int i = min_th; i <= max_th; ++i) 
	{
	  _fg_sum += _fg->GetBinContent(i);
	  _bg_sum += _bg->GetBinContent(i);
	}

      double scale_factor = (_bg_sum > 0) ? (_fg_sum / _bg_sum) : 1.0;
      if (_bg_sum > 0) _bg->Scale(scale_factor);

      // Save histograms
      outputFile->cd();
      _fg->Write();
      _bg->Write();

      // Move binStart to the next unprocessed bin
      binStart = binEnd + 1;
    }

  outputFile->Close();
  inputFile->Close();

  std::cout << "All Done." << std::endl;
  gSystem->Exit(0);

}

