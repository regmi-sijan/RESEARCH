#include <fstream>
#include <string>
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include <vector>
#include <calib_emc_eta/CaloCalibEmc_eta.h>
#include "TStopwatch.h"

// load so object
R__LOAD_LIBRARY(libcalibCaloEmc_eta.so)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// set of cut parameters and trigger vectors
// photon level cuts
double pt1_cut    = 1.0;
double pt2_cut    = 0.9;
double chi2_cut   = 4.0;
  
// pair-wise cuts
double alpha_cut  = 0.6;
double delR_cut   = 1.1;
  
// event level cuts
int max_nCluster_cut = 100; 
int min_nCluster_cut = 2;
double zver_cut   = 200.0;

// cuts for event mixing
double diff_zver_cut = 10.0;
int bckgnd_evnts = 1000;
int diff_nClus_cut = 15; // alternative to Centrality cut

// cuts for position swapping
bool doSwapBkg = false;
int n_Swap = 2;
double delR_bkg = 0.1;
double delE_bkg = 0.1;

// information on triggers
std::vector<int> desired_triggers = {10, 11, 12, 13, 14}; // for min-bias
// when we include all photon triggers {36, 37, 38, 24, 25, 26, 27, 28, 29, 30, 31};


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// converting cut parameter to vectors
std::vector<double> CutParms = {
  pt1_cut,
  pt2_cut,
  chi2_cut,
  alpha_cut,
  delR_cut,
  static_cast<double> (max_nCluster_cut),
  static_cast<double> (min_nCluster_cut),
  zver_cut
};

// converting event mixing cut parameter to vectors
std::vector<double> EvtMixCutParms = {
  diff_zver_cut,
  static_cast<double> (bckgnd_evnts),
  static_cast<double> (diff_nClus_cut)
};

// converting position swapping cut parameter to vectors
std::vector<double> PosSwapParms = {
  static_cast<double> (doSwapBkg),
  static_cast<double> (n_Swap),
  delR_bkg,
  delE_bkg
};

// helper function for GetTChain
TTree* GetTChain(TString ifile="")
{
  ifstream inFile;
  inFile.open(ifile, ios::in); // open to read
  if (!inFile)
    {
      cerr << "Unable to open file:\n ";
      exit(1);
    }

  TChain *etatree = new TChain("_eventTree");
  string root_file;
  int lines=0;
  while (std::getline(inFile, root_file))
    {
      etatree->Add(root_file.c_str());
      lines += 1;
    }
  printf("total lines: %d\n",lines);
  inFile.close();
	
  return (TTree *) etatree;
}

// RunMode is defined to run in either of three mode
enum RunMode {
  kSignal = 0,
  kEventMixing  = 1,
  kPositionSwapping = 2
};

// main loop (macro) code for analysis
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void RunEtaLoop(int n_file = 1)
{
  // start of timer
  TStopwatch timer;
  timer.Start();

  int mode = 1; // or use "0" "1" or "2"
	
  // setting number of events we want to run
  int nevents = -1;
	
  // getting input files (file-list) and output file
  const std::string in_fname = "Run24_pp_Ntuples_";

  const std::string macro_loc = "/sphenix/user/sregmi/WORKING_AREA/ETA_MESON/macro";

  // check this and down than this to confirm input filelist
  const std::string ifile= macro_loc + "/Run24_pp_filelist/" + in_fname + std::to_string(n_file); // input file-list
	
  const std::string ofile= macro_loc + "/RESULTS/run24_pp_photrig_evtmix_" + std::to_string(n_file) + ".root";

  std::cout << " input file list " << ifile << std::endl;
  std::cout << " output file list " << ofile << std::endl;
	

  TTree* intree1 =  GetTChain(ifile.c_str());

  // calling out function(s) object
  CaloCalibEmc_eta calo_obj("CaloCalibEmc_eta", ofile.c_str());
  
  calo_obj.InitRun(0);
	
  // running in different mode based on our input
  switch (mode) {
  case kSignal:
    calo_obj.Loop(desired_triggers, CutParms, nevents, ifile.c_str(), intree1); // creates foreground
    break;

  case kEventMixing:
    calo_obj.Loop_Event_Mixing(EvtMixCutParms, desired_triggers, CutParms, nevents, ifile.c_str(), intree1);
    //calo_obj.Loop_Event_Mixing_Old(EvtMixCutParms, desired_triggers, CutParms, nevents, ifile.c_str(), intree1);
    break;

  case kPositionSwapping:
    calo_obj.Loop_Position_Swapping(PosSwapParms, desired_triggers, CutParms, nevents, ifile.c_str(), intree1);
    break;

  default:
    std::cerr << "Unknown mode: " << mode << "\n";
    std::exit(1);
  }
	
  calo_obj.End(0);
  std::cout << "All Done" << std::endl;
	
  // print elapsed time
  timer.Stop();
  std::cout << "  Wall‑clock time : " << timer.RealTime() << std::endl;
  std::cout << "  CPU time        : " << timer.CpuTime()  << std::endl;

  gSystem->Exit(0);

}

