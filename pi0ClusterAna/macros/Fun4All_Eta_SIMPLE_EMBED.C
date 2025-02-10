#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include <DisplayOn.C>
#include <G4Setup_sPHENIX.C>
#include <G4_Mbd.C>
#include <G4_CaloTrigger.C>
#include <G4_Centrality.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_KFParticle.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
#include <Trkr_Reco.C>
#include <Trkr_Eval.C>
#include <Trkr_QA.C>

#include <Trkr_Diagnostics.C>
#include <G4_User.C>
#include <QA.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <g4centrality/PHG4CentralityReco.h>
 
//#include <calib_emc_eta/CaloCalibEmc_eta.h>
#include <pi0clusterana/pi0ClusterAna.h>

#include <iostream>
#include <string> 

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
//R__LOAD_LIBRARY(libcaloCalibDBFile.so)

//R__LOAD_LIBRARY(libcentrality_io.so)
//R__LOAD_LIBRARY(libg4centrality.so)

//R__LOAD_LIBRARY(libcalibCaloEmc_eta.so)
R__LOAD_LIBRARY(libpi0ClusterAna.so)

// For HepMC Hijing
// try inputFile = /sphenix/sim/sim01/sphnxpro/sHijing_HepMC/sHijing_0-12fm.dat

// this macro can be used to read DST file and use information for simple particle generation using different modes and embedded them
int Fun4All_Eta_SIMPLE_EMBED(
			     string inputFile = "DST_TRUTH_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-000000.root"
			     //string inputFileList = "/sphenix/user/sregmi/WORKING_AREA/truth_par_ana/EtaDecayAna/macros/FileList/dst_truth_000100"
			     )
{
  const int nEvents = 2;
	
  //const int seg_num = 0;
  const int skip = 0;
  const string &outdir = ".";
	
  // list of all input files
  //string inputFile0 = "DST_CALO_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-";
  //string inputFile1 = "DST_GLOBAL_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-";
  //string inputFile2 = "DST_TRUTH_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-";  
  //string inputFile3 = "DST_CALO_CLUSTER_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-";
  //string inputFile4 = "DST_MBD_EPD_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-";
  //string inputFile2 = "DST_BBC_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000019-";

  // input filelist
  std::string SegNum = inputFile.substr(inputFile.length() - 10, 5); 
  // list of all output files
  //string outputFile0 = "/sphenix/user/sregmi/WORKING_AREA/truth_par_ana/pi0ClusterAna_ref/macros/RESULTS/test_Fun4All_ntuples_cluster_run19_multi_embedding_";
  string outputFile1 = "/sphenix/user/sregmi/WORKING_AREA/truth_par_ana/pi0ClusterAna/RESULTS/Fun4All_ntuples_truth_run19_eta_embedding_" + SegNum + ".root";
 
 
  std::cout << "InputFile List = " << inputFile << std::endl;
  std::cout << "Output file = " << outputFile1 << std::endl;
  
  /*
  // setting all segments as per their number
  int ynum_int = 1000000+ seg_num;
  TString yn_tstr = "";
  yn_tstr += ynum_int;
  yn_tstr.Remove(0,1);

  // final version of all input files
  //inputFile0 += yn_tstr.Data();
  //inputFile1 += yn_tstr.Data();
  inputFile2 += yn_tstr.Data();
  //inputFile3 += yn_tstr.Data();
  //inputFile4 += yn_tstr.Data();
	
  //inputFile0 += ".root";
  //inputFile1 += ".root";
  inputFile2 += ".root";
  //inputFile3 += ".root";
  //inputFile4 += ".root";
  
  // final version of all output files 
  //outputFile0 += yn_tstr.Data();
  outputFile1 += yn_tstr.Data();
	
  //outputFile0 += ".root";
  outputFile1 += ".root";
	 
  cout << "running over these files" << endl;
  //cout << inputFile0 << endl;
  //cout << inputFile1 << endl;
  cout << inputFile2 << endl;
  //cout << inputFile3 << endl;
  //cout << inputFile4 << endl;

  cout << "==============================================================================================" << endl;

  cout << "our output files will be as followings" << endl;
  //cout << outputFile0 << endl;
  cout << outputFile1 << endl;
  */
  //===========================================================================================================
  //===========================================================================================================

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);


  //Opt to print all random seed used for debugging reproducibility. Comment out to reduce stdout prints.
  PHRandomSeed::Verbosity(0);

  // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();

  //===============
  // conditions DB flags
  //===============
  Enable::CDB = true;
  // global tag
  rc->set_StringFlag("CDB_GLOBALTAG","MDC2");
  // // 64 bit timestamp
  rc->set_uint64Flag("TIMESTAMP",6);

  //===============
  // Input options
  //===============
  // verbosity setting (applies to all input managers)
  Input::VERBOSITY = 1;

  Input::READHITS = false; // true;
  //INPUTREADHITS::filename[0] = inputFile0;
  //INPUTREADHITS::filename[1] = inputFile1;

  Input::EMBED = true;
  //INPUTEMBED::filename[0] = inputFile0;  //0;
  //INPUTEMBED::filename[1] = inputFile1;
  //INPUTEMBED::filename[2] = inputFile2;
  //INPUTEMBED::filename[3] = inputFile3;
  
  INPUTEMBED::filename[0] = inputFile;
  //INPUTEMBED::listfile[0] = inputFileList; // done for input file list
  
  // for simple particle generator
  Input::SIMPLE = true;
  Input::SIMPLE_VERBOSITY = 0;
  Input::SIMPLE_NUMBER = 1; // 12 if there are 12 particles to embedded

  /*
  // while using particle gun
  Input::GUN = true;
  Input::GUN_NUMBER = 1; // if you need "N" number of them
  Input::GUN_VERBOSITY = 1;
  */


  InputInit();

  //--------------
  // Set generator specific options
  //--------------
  // can only be set after InputInit() is called

  // Simple Input generator:
  // if you run more than one of these Input::SIMPLE_NUMBER > 1
  // add the settings for other with [1], next with [2]...
	
  if (Input::SIMPLE)
    {
      // make sure eta meson (221) is embedded at the end (this is to track it right at first) to figure out the era of embedding
      INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("eta", 1);
		
      /*
	INPUTGENERATOR::SimpleEventGenerator[0]->add_particles(130,1); // K0 long
	INPUTGENERATOR::SimpleEventGenerator[1]->add_particles(211,1); // pi+ ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[2]->add_particles(310,1); // K0 short ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[3]->add_particles(321,1); // K+ ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[4]->add_particles(411,1); // D+ ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[5]->add_particles(421,1); // D0 ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[6]->add_particles(431,1); // D+ s ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[7]->add_particles(-321,1); // K- ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[8]->add_particles(-411,1); // D- ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[9]->add_particles(-421,1); // anti D0 ;; we can use pdg number of particle instead of name
	INPUTGENERATOR::SimpleEventGenerator[10]->add_particles(-431,1); // anti D s ;; we can use pdg number of particle instead of name		
	INPUTGENERATOR::SimpleEventGenerator[11]->add_particles(221,1); // eta meson ;; we can use pdg number of particle instead of name
      */

      if (Input::HEPMC || Input::EMBED)
	{
	  INPUTGENERATOR::SimpleEventGenerator[0]->set_reuse_existing_vertex(true);
	  INPUTGENERATOR::SimpleEventGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);

	  /*
	    INPUTGENERATOR::SimpleEventGenerator[0]->set_reuse_existing_vertex(true);    
	    INPUTGENERATOR::SimpleEventGenerator[1]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[2]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[3]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[4]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[5]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[6]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[7]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[8]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[9]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[10]->set_reuse_existing_vertex(true);
	    INPUTGENERATOR::SimpleEventGenerator[11]->set_reuse_existing_vertex(true);
			
      
	    INPUTGENERATOR::SimpleEventGenerator[0]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[1]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[2]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[3]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[4]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[5]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[6]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[7]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[8]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[9]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[10]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	    INPUTGENERATOR::SimpleEventGenerator[11]->set_existing_vertex_offset_vector(0.0, 0.0, 0.0);
	  */

	}
      
      else
	{
	  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
	  INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0.01, 0.01, 5.);

	  /*
	    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
						
	    INPUTGENERATOR::SimpleEventGenerator[1]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[2]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[3]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[4]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[5]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[6]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[7]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[8]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[9]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[10]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
	    INPUTGENERATOR::SimpleEventGenerator[11]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
			
      			
	    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);			
	    INPUTGENERATOR::SimpleEventGenerator[1]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[2]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[3]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[4]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[5]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[6]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[7]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[8]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[9]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[10]->set_vertex_distribution_mean(0., 0., 0.);
	    INPUTGENERATOR::SimpleEventGenerator[11]->set_vertex_distribution_mean(0., 0., 0.);

	    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[1]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[2]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[3]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[4]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[5]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[6]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[7]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[8]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[9]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[10]->set_vertex_distribution_width(0.01, 0.01, 5.);
	    INPUTGENERATOR::SimpleEventGenerator[11]->set_vertex_distribution_width(0.01, 0.01, 5.);
	  */

	}
				
      INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1.1, 1.1);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
      INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(5, 20);
      //INPUTGENERATOR::SimpleEventGenerator[0]->set_power_law_n(-6.5);
 
      /*
	INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1.1, 1.1);		
	INPUTGENERATOR::SimpleEventGenerator[1]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[2]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[3]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[4]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[5]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[6]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[7]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[8]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[9]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[10]->set_eta_range(-1.1, 1.1);
	INPUTGENERATOR::SimpleEventGenerator[11]->set_eta_range(-1.1, 1.1);
		

	INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);		
	INPUTGENERATOR::SimpleEventGenerator[1]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[2]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[3]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[4]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[5]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[6]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[7]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[8]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[9]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[10]->set_phi_range(-M_PI, M_PI);
	INPUTGENERATOR::SimpleEventGenerator[11]->set_phi_range(-M_PI, M_PI);
		

	INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[1]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[2]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[3]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[4]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[5]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[6]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[7]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[8]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[9]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[10]->set_pt_range(1, 20);
	INPUTGENERATOR::SimpleEventGenerator[11]->set_pt_range(1, 20);
		
			
	INPUTGENERATOR::SimpleEventGenerator[0]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[1]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[2]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[3]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[4]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[5]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[6]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[7]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[8]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[9]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[10]->set_power_law_n(-6.5);
	INPUTGENERATOR::SimpleEventGenerator[11]->set_power_law_n(-6.5);
      */

    }

  /*
  // particle gun
  // if you run more than one of these Input::GUN_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  if (Input::GUN)
  {
  INPUTGENERATOR::Gun[0]->AddParticle("310", 0, 1, 0);
  INPUTGENERATOR::Gun[0]->set_vtx(0, 0, 0);
  }
  */


  if (Input::PILEUPRATE > 0)
    {
      //! apply sPHENIX nominal beam parameter with 2mrad crossing as defined in sPH-TRG-2020-001
      Input::ApplysPHENIXBeamParameter(INPUTMANAGER::HepMCPileupInputManager);
    }
  // register all input generators with Fun4All
  InputRegister();

  // set up production relatedstuff
  Enable::PRODUCTION = true;

  //======================
  // Write the DST
  //======================
  /*
    Enable::DSTOUT = false;
    Enable::DSTOUT_COMPRESS = false;
    DstOut::OutputDir = outdir;
    DstOut::OutputFile = outputFile0;
  */

  //Option to convert DST to human command readable TTree for quick poke around the outputs
  //Enable::DSTREADER = true;

  // turn the display on (default off)
  Enable::DISPLAY = false;

  //======================
  // What to run
  //======================

  // QA, main switch
  Enable::QA = false;

  //Global options (enabled for all enables subsystems - if implemented)
  //Enable::ABSORBER = true;
  //Enable::OVERLAPCHECK = true;
  //Enable::VERBOSITY = 1;

  //Enable::MBD = true;
  //Enable::MBD_SUPPORT = true; // save hist in MBD/BBC support structure
  //Enable::MBDFAKE = true;  // Smeared vtx and t0, use if you don't want real MBD/BBC in simulation

  //Enable::PIPE = true;
  //Enable::PIPE_ABSORBER = true;
  //Enable::INTT = false;
  //Enable::INTT_ABSORBER = true; // enables layerwise support structure readout
  //Enable::INTT_SUPPORT = true; // enable global support structure readout
  Enable::INTT_CELL = Enable::INTT && true;
  Enable::INTT_CLUSTER = Enable::INTT_CELL && true;
  Enable::INTT_QA = Enable::INTT_CLUSTER && Enable::QA && true;

  Enable::TPC = false;
  Enable::TPC_ABSORBER = true;
  Enable::TPC_CELL = Enable::TPC && true;
  Enable::TPC_CLUSTER = Enable::TPC_CELL && true;
  Enable::TPC_QA = Enable::TPC_CLUSTER && Enable::QA && true;

  Enable::MICROMEGAS = false;
  Enable::MICROMEGAS_CELL = Enable::MICROMEGAS && true;
  Enable::MICROMEGAS_CLUSTER = Enable::MICROMEGAS_CELL && true;
  Enable::MICROMEGAS_QA = Enable::MICROMEGAS_CLUSTER && Enable::QA && true;

  Enable::TRACKING_TRACK = false;
  Enable::TRACKING_EVAL = Enable::TRACKING_TRACK && true;
  Enable::TRACKING_QA = Enable::TRACKING_TRACK && Enable::QA && true;

  //  cemc electronics + thin layer of W-epoxy to get albedo from cemc
  //  into the tracking, cannot run together with CEMC
  //  Enable::CEMCALBEDO = true;

  Enable::CEMC = true;
  // Enable::CEMC_ABSORBER = true;
  Enable::CEMC_CELL = Enable::CEMC && true;
  Enable::CEMC_TOWER = Enable::CEMC_CELL && true;
  Enable::CEMC_CLUSTER = Enable::CEMC_TOWER && true;
  //Enable::CEMC_EVAL = false;//Enable::CEMC_CLUSTER && true;
  //Enable::CEMC_QA = false;//Enable::CEMC_CLUSTER && Enable::QA && true;

  Enable::HCALIN =false;
  Enable::HCALIN_ABSORBER = true;
  Enable::HCALIN_CELL = Enable::HCALIN && true;
  Enable::HCALIN_TOWER = Enable::HCALIN_CELL && true;
  Enable::HCALIN_CLUSTER = Enable::HCALIN_TOWER && true;
  Enable::HCALIN_EVAL = Enable::HCALIN_CLUSTER && true;
  Enable::HCALIN_QA = Enable::HCALIN_CLUSTER && Enable::QA && true;

	
  //Enable::MAGNET = false;
  //Enable::MAGNET_ABSORBER = false;


  Enable::HCALOUT = false;
  Enable::HCALOUT_ABSORBER = true;
  Enable::HCALOUT_CELL = Enable::HCALOUT && true;
  Enable::HCALOUT_TOWER = Enable::HCALOUT_CELL && true;
  Enable::HCALOUT_CLUSTER = Enable::HCALOUT_TOWER && true;
  Enable::HCALOUT_EVAL = Enable::HCALOUT_CLUSTER && true;
  Enable::HCALOUT_QA = Enable::HCALOUT_CLUSTER && Enable::QA && true;

  Enable::EPD = false;

  Enable::BEAMLINE = true;
  //  Enable::BEAMLINE_ABSORBER = true;  // makes the beam line magnets sensitive volumes
  //  Enable::BEAMLINE_BLACKHOLE = true; // turns the beamline magnets into black holes
  Enable::ZDC = false;
  //  Enable::ZDC_ABSORBER = true;
  //  Enable::ZDC_SUPPORT = true;
  Enable::ZDC_TOWER = Enable::ZDC && true;
  Enable::ZDC_EVAL = Enable::ZDC_TOWER && true;

  //! forward flux return plug door. Out of acceptance and off by default.
  //Enable::PLUGDOOR = true;
  //Enable::PLUGDOOR_ABSORBER = true;

  //Enable::GLOBAL_RECO = true;
  //Enable::GLOBAL_FASTSIM = true;

  //Enable::KFPARTICLE = true;
  //Enable::KFPARTICLE_VERBOSITY = 1;
  //Enable::KFPARTICLE_TRUTH_MATCH = true;
  //Enable::KFPARTICLE_SAVE_NTUPLE = true;

  Enable::CALOTRIGGER = Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER && false;

  Enable::JETS = false;
  Enable::JETS_EVAL = Enable::JETS && true;
  Enable::JETS_QA = Enable::JETS && Enable::QA && true;

  // HI Jet Reco for p+Au / Au+Au collisions (default is false for
  // single particle / p+p-only simulations, or for p+Au / Au+Au
  // simulations which don't particularly care about jets)
  Enable::HIJETS = false && Enable::JETS && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER;

  // 3-D topoCluster reconstruction, potentially in all calorimeter layers
  Enable::TOPOCLUSTER = false && Enable::CEMC_TOWER && Enable::HCALIN_TOWER && Enable::HCALOUT_TOWER;
  // particle flow jet reconstruction - needs topoClusters!
  //Enable::PARTICLEFLOW = true && Enable::TOPOCLUSTER;
  // centrality reconstruction
  //Enable::CENTRALITY = true;

  // new settings using Enable namespace in GlobalVariables.C
  //Enable::BLACKHOLE = true;
  //Enable::BLACKHOLE_SAVEHITS = false; // turn off saving of bh hits
  //BlackHoleGeometry::visible = true;

  // Initialize the selected subsystems
  G4Init();

  //---------------------
  // GEANT4 Detector description
  //---------------------

		
  if (!Input::READHITS)
    {
      G4Setup();
    }
		

  //------------------
  // Detector Division
  //------------------

  if (Enable::MBD || Enable::MBDFAKE) Mbd_Reco();

  if (Enable::MVTX_CELL) Mvtx_Cells();
  if (Enable::INTT_CELL) Intt_Cells();
  if (Enable::TPC_CELL) TPC_Cells();
  if (Enable::MICROMEGAS_CELL) Micromegas_Cells();

  if (Enable::CEMC_CELL) CEMC_Cells();

  if (Enable::HCALIN_CELL) HCALInner_Cells();

  if (Enable::HCALOUT_CELL) HCALOuter_Cells();

  //-----------------------------
  // CEMC towering and clustering
  //-----------------------------

  if (Enable::CEMC_TOWER) CEMC_Towers();
  if (Enable::CEMC_CLUSTER) CEMC_Clusters();

  //-----------------------------
  // HCAL towering and clustering
  //-----------------------------

  if (Enable::HCALIN_TOWER) HCALInner_Towers();
  if (Enable::HCALIN_CLUSTER) HCALInner_Clusters();

  if (Enable::HCALOUT_TOWER) HCALOuter_Towers();
  if (Enable::HCALOUT_CLUSTER) HCALOuter_Clusters();

  // if enabled, do topoClustering early, upstream of any possible jet reconstruction
  if (Enable::TOPOCLUSTER) TopoClusterReco();

  //--------------
  // SVTX tracking
  //--------------
  if(Enable::TRACKING_TRACK)
    {
      TrackingInit();
    }
  if (Enable::MVTX_CLUSTER) Mvtx_Clustering();
  if (Enable::INTT_CLUSTER) Intt_Clustering();
  if (Enable::TPC_CLUSTER) TPC_Clustering();
  if (Enable::MICROMEGAS_CLUSTER) Micromegas_Clustering();

  if (Enable::TRACKING_TRACK)
    {
      Tracking_Reco();
    }
  //-----------------
  // Global Vertexing
  //-----------------
  /*
    if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
    {
    cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << endl;
    gSystem->Exit(1);
    }
    if (Enable::GLOBAL_RECO)
    {
    Global_Reco();
    }
    else if (Enable::GLOBAL_FASTSIM)
    {
    Global_FastSim();
    }
  */	

  //-----------------
  // Centrality Determination
  //-----------------

  if (Enable::CENTRALITY)
    {
      Centrality();
    }

  //-----------------
  // Calo Trigger Simulation
  //-----------------

  if (Enable::CALOTRIGGER)
    {
      CaloTrigger_Sim();
    }

  //---------
  // Jet reco
  //---------

  if (Enable::JETS) Jet_Reco();
  if (Enable::HIJETS) HIJetReco();

  if (Enable::PARTICLEFLOW) ParticleFlow();

  /*
  //----------------------
  // Simulation evaluation
  //----------------------
  string outputroot = outputFile0;
  string remove_this = ".root";
  size_t pos = outputroot.find(remove_this);
  if (pos != string::npos)
  {
  outputroot.erase(pos, remove_this.length());
  }

  if (Enable::TRACKING_EVAL) Tracking_Eval(outputroot + "_g4svtx_eval.root");

  if (Enable::CEMC_EVAL) CEMC_Eval(outputroot + "_g4cemc_eval.root");
  
  if (Enable::HCALIN_EVAL) HCALInner_Eval(outputroot + "_g4hcalin_eval.root");

  if (Enable::HCALOUT_EVAL) HCALOuter_Eval(outputroot + "_g4hcalout_eval.root");

  if (Enable::JETS_EVAL) Jet_Eval(outputroot + "_g4jet_eval.root");

  if (Enable::DSTREADER) G4DSTreader(outputroot + "_DSTReader.root");

  if (Enable::USER) UserAnalysisInit();
  */
  
  /*
  //======================
  // Run KFParticle on evt
  //======================
  if (Enable::KFPARTICLE && Input::UPSILON) KFParticle_Upsilon_Reco();
  if (Enable::KFPARTICLE && Input::DZERO) KFParticle_D0_Reco();
  */
	
  //----------------------
  // Standard QAs
  //----------------------

  if (Enable::CEMC_QA) CEMC_QA();
  if (Enable::HCALIN_QA) HCALInner_QA();
  if (Enable::HCALOUT_QA) HCALOuter_QA();

  if (Enable::JETS_QA) Jet_QA();

  if (Enable::MVTX_QA) Mvtx_QA();
  if (Enable::INTT_QA) Intt_QA();
  if (Enable::TPC_QA) TPC_QA();
  if (Enable::MICROMEGAS_QA) Micromegas_QA();
  if (Enable::TRACKING_QA) Tracking_QA();

  if (Enable::TRACKING_QA && Enable::CEMC_QA && Enable::HCALIN_QA && Enable::HCALOUT_QA) QA_G4CaloTracking();

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  if (Enable::PRODUCTION)
    {
      Production_CreateOutputDir();
    }

  if (Enable::DSTOUT)
    {
      string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;
      Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
      if (Enable::DSTOUT_COMPRESS)
	{
	  ShowerCompress();
	  DstCompress(out);
	}
      se->registerOutputManager(out);
    }
  //-----------------
  // Event processing
  //-----------------
  if (Enable::DISPLAY)
    {
      DisplayOn();

      gROOT->ProcessLine("Fun4AllServer *se = Fun4AllServer::instance();");
      gROOT->ProcessLine("PHG4Reco *g4 = (PHG4Reco *) se->getSubsysReco(\"PHG4RECO\");");

      cout << "-------------------------------------------------" << endl;
      cout << "You are in event display mode. Run one event with" << endl;
      cout << "se->run(1)" << endl;
      cout << "Run Geant4 command with following examples" << endl;
      gROOT->ProcessLine("displaycmd()");

      return 0;
    }

  // if we decide to use CaloCalib Class for cluster/tower level processing
	
  /*
    CaloCalibEmc_eta *eval_pi0 = new CaloCalibEmc_eta("dummy0", outputFile0);
    // this call is needed for embedding
    //eval_pi2->set_centrality_nclusters_cut(350);  // which uses more central events
    // than we will for data to enhance Bkg
    // to match the enhanced signal from embed
    se->registerSubsystem(eval_pi0);
    cout << "successful registration of pi0 " << endl;
  */


  /*
    Fun4AllInputManager *truthCalo = new Fun4AllDstInputManager("DSTCaloTruth");
    std::cout << "Adding truth file  " << std::endl;
    truthCalo->AddFile(inputFile2);
    //truthCalo -> AddListFile(truthFile,1);
    se -> registerInputManager(truthCalo);
  */

  
  pi0ClusterAna *eval = new pi0ClusterAna("pi0ClusterAna", outputFile1);
  eval->SaveAllMode(true); //"false" default means we save only final state sec. par.
  se->registerSubsystem(eval);
  cout << "successful analysis of truth particle " << endl;


  // if we run the particle generator and use 0 it'll run forever
  // for embedding it runs forever if the repeat flag is set
  if (nEvents == 0 && !Input::HEPMC && !Input::READHITS && INPUTEMBED::REPEAT)
    {
      cout << "using 0 for number of events is a bad idea when using particle generators" << endl;
      cout << "it will run forever, so I just return without running anything" << endl;
      return 0;
    }

  // if we use a negative number of events we go back to the command line here
  if (nEvents == 0)
    {
      return 0;
    }

  se->skip(skip);
  se->run(nEvents);

  //-----
  // QA output
  //-----

  //if (Enable::QA) QA_Output(outputroot + "_qa.root");

  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;
  /*
    if (Enable::PRODUCTION)
    {
    Production_MoveOutput();
    }
  */
  gSystem->Exit(0);
  return 0;
}
#endif
