// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PI0CLUSTERANA_H
#define PI0CLUSTERANA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <vector>

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

class pi0ClusterAna : public SubsysReco
{
 public:

  pi0ClusterAna(const std::string &name, const std::string &outName);

  ~pi0ClusterAna() override;

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


 private:

  float getEta(PHG4Particle *particle);
    
  TTree *clusters_Towers;
  TTree *truth_photon;
  TTree *truth_pi0;
  
  //CaloEvalStack *caloevalstack;
  
  //stuff for towers and clusters
  std::vector<float> m_eta_center;
  std::vector<float> m_phi_center;
  std::vector<float> m_tower_energy;
  std::vector<float> m_cluster_eta;
  std::vector<float> m_cluster_phi;
  std::vector<float> m_cluster_e;
  std::vector<float> m_cluster_chi2;
  std::vector<float> m_cluster_prob;
  std::vector<float> m_cluster_nTowers;
  
  //stuff for truth photons
  std::vector<float> m_asym;
  std::vector<float> m_deltaR;
  std::vector<float> m_lead_E;
  std::vector<float> m_sublead_E;
  std::vector<float> m_lead_phi;
  std::vector<float> m_lead_eta;
  std::vector<float> m_sublead_phi;
  std::vector<float> m_sublead_eta;


  float alphaCut = -1.;
	int n_event;
	
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

  TFile *out;
	TFile *f_temp;
  //Fun4AllHistoManager *hm = nullptr;
	
  std::string Outfile;

	TH3F *Pi0_pt_eta_mass = nullptr;
	TH3F *Eta_pt_eta_mass = nullptr;
	TH1F *Pi0_pt = nullptr;
	TH1F *Eta_pt = nullptr;
	TH1F *pairInvMassTotal = nullptr;

};

#endif // PI0CLUSTERANA_H
