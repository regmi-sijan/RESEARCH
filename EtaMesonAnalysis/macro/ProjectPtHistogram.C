#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TString.h>
#include <vector>


void ProjectPtHistogram()
{
   // input and output files
   const std::string input_file = "_cent_combined.root";
	
   const std::string output_file = "all_pairpT_projected_" + input_file;
	
   std::vector<double> Lows = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
   std::vector<double> Highs  = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	
   // Ensure the vectors have the same size
   if (Lows.size() != Highs.size()) {
     std::cerr << "Error: The vectors of lower and upper limits must have the same size." << std::endl;
     return;
   }

  // Load the 3D histogram from a file
  TFile *file = new TFile(input_file.c_str(), "READ");

  TH3F *my3DHist_f = (TH3F*)file->Get("pairInvMassPtEta"); // which histogram we want to process
  TH3F *my3DHist_b = (TH3F*)file->Get("pairInvMassPtEtaBkgd"); // which histogram we want to process

  if ((!my3DHist_f) || (!my3DHist_b)) {
    std::cerr << "Either of Histogram not found!" << std::endl;
    return;
  }

  // Open a ROOT file to save all projected histograms
  TFile *outFile = new TFile(output_file.c_str(), "RECREATE");

  // Loop through the ranges
  for (size_t i = 0; i < Lows.size(); ++i) {
		
    // For "foreground"
    // Define the X range for projection
    int yf_min = my3DHist_f->GetYaxis()->FindBin(Lows[i]); // xLow is the lower bound of the X range
    int yf_max = my3DHist_f->GetYaxis()->FindBin(Highs[i]); // xHigh is the upper bound of the X range

    int zf_min = my3DHist_f->GetZaxis()->FindBin(my3DHist_f->GetZaxis()->GetXmin());
    int zf_max = my3DHist_f->GetZaxis()->FindBin(my3DHist_f->GetZaxis()->GetXmax());

    // projecting 3D Hist to 1 D hist
    TH1* projYZ_f = my3DHist_f->ProjectionX(Form("f_%zu", i), yf_min, yf_max, zf_min, zf_max, "C");
    projYZ_f->SetTitle("InvMassHistForPairpT");
													
    // Save the projected histogram to the open ROOT file
    projYZ_f->Write(Form("f_%zu", i));

		
    // For "background"
    // Define the X range for projection
    int yb_min = my3DHist_b->GetYaxis()->FindBin(Lows[i]); // xLow is the lower bound of the X range
    int yb_max = my3DHist_b->GetYaxis()->FindBin(Highs[i]); // xHigh is the upper bound of the X range

    int zb_min = my3DHist_b->GetZaxis()->FindBin(my3DHist_b->GetZaxis()->GetXmin());
    int zb_max = my3DHist_b->GetZaxis()->FindBin(my3DHist_b->GetZaxis()->GetXmax());

    // projecting 3D Hist to 1 D hist
    TH1* projYZ_b = my3DHist_b->ProjectionX(Form("b_%zu", i), yb_min, yb_max, zb_min, zb_max, "C");
    projYZ_b->SetTitle("InvMassHistForPairpT");
													
    // Save the projected histogram to the open ROOT file
    projYZ_b->Write(Form("b_%zu", i));

  }

  // Clean up
  outFile->Close();
  file->Close();
  delete file;
  delete outFile;
		
  std::cout << "All Done" << std::endl;
  gSystem->Exit(0);
}


