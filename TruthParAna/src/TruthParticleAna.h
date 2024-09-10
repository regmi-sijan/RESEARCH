// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef TRUTHPARTICLEANA_H
#define	TRUTHPARTICLEANA_H

#include <fun4all/SubsysReco.h>
#include <TLorentzVector.h>
#include <string>
#include <vector>
#include <THn.h>
#include <THnSparse.h>
#include <iostream>

class PHCompositeNode;
class PHG4Particle;
class CaloEvalStack;
class Fun4AllHistoManager;
class TFile;
class TH1F;
class TH2F;
class TH3F;
class TF1;
class TH1;
class TNtuple;
class TTree;
class TString;
class TCanvas;

class TruthParticleAna : public SubsysReco
{
 public:

  TruthParticleAna(const std::string &name, const std::string &outName);

  ~TruthParticleAna() override;

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

	void Loop(int nevts, TString _filename, TTree * intree = 0);
	void Loop_back_up(int nevts, TString _filename, TTree * intree = 0);
	void Loop_ETAtoPIONS(int nevts, TString _filename, TTree * intree = 0);
	void Loop_background_position_swappping(int nevts, TString _filename, TTree * intree = 0);

	void Loop_background_position_swappping_v2(int nevts, TString _filename, TTree * intree = 0);
	
	void Loop_background_event_mixing(int nevts, TString _filename, TString ref_file,TTree * intree = 0);
	void Loop_final_processing(int nevts, TString _filename, TString _outputfile, TTree * intree = 0);
	
	void Loop_background_position_swapppingF(int nevts, TString _filename, TTree * intree = 0);
	void Loop_background_position_swapppingB(int nevts, TString _filename, TTree * intree = 0);
	void Loop_background_90deg_rotation(int nevts, TString _filename, TTree * intree = 0);
	void check_primary_secondary(int nevts, TString _filename, TTree * intree = 0);

 private:

	THnF *nD_InvMass_f; 
	THnF *nD_InvMass_b; 

  float getEta(PHG4Particle *particle);
  float FourVecEta(TLorentzVector& four_vector);

  TTree *clusters_Towers;
  TTree *truth_photon;
  TTree *truth_pi0;

  float alphaCut = -1.;
	int n_event;
	int ndim; 
  
	TTree *_eventTree = nullptr;
	
	int _eventNumber = -1;
	int _nFourVector = -1; // store the number of clusters recorded in each event (our focus will only be in EMCal)
  float _fv_energy[100000] = {0}; // truth energy of four vector
  float _fv_px[100000] = {0}; // 3 of four momentum (px, py and pz) of the cluster
  float _fv_py[100000] = {0};
  float _fv_pz[100000] = {0};
  float _fv_Eta[100000] = {0}; // pseudorapidity of the four vector
  float _pid_primary[100000] = {0}; // associated primary particle
  float _pid_secondary[100000] = {0}; // associated secondary particle
	float _primary_id[100000] = {0}; // primary id of the particle
	float _parent_id[100000] = {0}; // parent id of the particle
	float _embedding[100000] = {0}; // parent id of the particle

	TTree *_data = nullptr;
	
	int _dataNum = -1;
	int _nPair = -1; // store the number of pairs for foreground or background
  int _fb_flag[100000] = {0}; // flag for foreground or background (1 for foreground and -1 for background
  // multiple (different) information for pairs
	float _InvMass[100000] = {0};
  float _pT[100000] = {0};
  float _Eta[100000] = {0};
  float _DelPhi[100000] = {0};
  float _DelEta[100000] = {0};

  
  TFile *out;
	TFile *f_temp;
  //Fun4AllHistoManager *hm = nullptr;
	
  std::string Outfile;

	TH3F *Pi0_pt_eta_mass = nullptr;
	TH3F *Eta_pt_eta_mass = nullptr;
	TH1F *Pi0_pt = nullptr;
	TH1F *Eta_pt = nullptr;
	TH1F *delphi_f = nullptr;
	TH1F *delphi_b = nullptr;
	TH1F *deleta_f = nullptr;
	TH1F *deleta_b = nullptr;
	
	TH1F *pairInvMassTotal = nullptr;
	TH1F *pairInvMassTotalBkgd = nullptr;

	TH3F *pairInvMassPtEta = nullptr;
	TH3F *pairInvMassPtEtaBkgd = nullptr;

	TH3F *pairpTDelPhiDelEta = nullptr;
	TH3F *pairpTDelPhiDelEtaBkgd  = nullptr;

	TH3F *pairInvMass_SecP1_SecP2 = nullptr;
	TH3F *pairInvMass_PriP1_PriP2 = nullptr;
	TH2F *check_PriP1_PriP2 = nullptr;

	TH1F *Pair_PT_f = nullptr; 
	TH1F *Pair_PT_b = nullptr; 

	TH1F *delR_f = nullptr;
	TH1F *delR_b = nullptr;

	TH3F *pairInvMassPt1Pt2 = nullptr;
	TH3F *pairInvMassPX1PX2 = nullptr;
	TH3F *pairInvMassPY1PY2 = nullptr;
	TH3F *pairInvMassPZ1PZ2 = nullptr;
	TH3F *pairInvMassE1E2 = nullptr;
	TH3F *pairInvMassEta1Eta2 = nullptr;

	TH3F *pairInvMassPtdelR = NULL;
	TH3F *pairInvMassPtdelRBkgd = NULL;
	
	TH2F *DelPhi_DelEta_f	= nullptr;
	TH2F *DelPhi_DelEta_b	= nullptr;
	
	TH3F *pairInvMasspTDelEta = nullptr;
	TH3F *pairInvMasspTDelEtaBkgd = nullptr;

	TH2F *DelPhi_PairpT_f = nullptr; 
	TH2F *DelPhi_PairpT_b = nullptr;
	
	TH2F *DelEta_PairpT_f = nullptr;
	TH2F *DelEta_PairpT_b = nullptr;
	
	TH3F *Pri1_Pri2_pairpT = nullptr;
	
	TH2F *DelR_pairpT_f = nullptr;
	TH2F *DelR_pairpT_b = nullptr;
};

#endif // TRUTHPARTICLEANA_H
