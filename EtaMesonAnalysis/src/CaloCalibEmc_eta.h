// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CALOCALIBEMC_ETA_H
#define CALOCALIBEMC_ETA_H

#include <fun4all/SubsysReco.h>
#include <globalvertex/GlobalVertex.h>
#include <deque>

#include <string>

class PHCompositeNode;
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

class CaloCalibEmc_eta : public SubsysReco
{
 public:
  CaloCalibEmc_eta(const std::string &name = "CaloCalibEmc_eta",\
	 const std::string &fnm = "outJF");

  virtual ~CaloCalibEmc_eta() {}

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

	void Loop(std::vector<int> desired_triggers, std::vector<double> cuts, 
					int nevts, TString _filename, TTree * intree);
	
	void Loop_Event_Mixing(std::vector<double> evt_mix_parms, 
					std::vector<int> desired_triggers, std::vector<double> CutParms, 
					int nevts, TString _filename, TTree * intree);

	void Loop_Event_Mixing_V2(std::vector<double> evt_mix_parms, 
					std::vector<int> desired_triggers, std::vector<double> CutParms, 
					int nevts, TString _filename, TTree * intree);

	void Loop_Event_Mixing_Old(std::vector<double> evt_mix_parms,
											 std::vector<int> desired_triggers, std::vector<double> CutParms,
											 int nevts, TString _filename, TTree * intree);

	void Loop_Position_Swapping(std::vector<double> pos_swap_parms, 
					std::vector<int> desired_triggers, std::vector<double> CutParms, 
					int nevts, TString _filename, TTree * intree);
  
	//void set_centrality_nclusters_cut(int n){m_cent_nclus_cut=n;}
  
	void set_GlobalVertexType(GlobalVertex::VTXTYPE type) 
		{
			m_use_vertextype = true;
			m_vertex_type = type;
		}

  void set_requireVertex(bool state)
		{
			reqVertex = state;
			return;
		}

 private:
  bool reqMinBias = true;
  bool reqVertex = false;
  
	int m_ievent = 0;
  std::string m_Filename;
  TFile *cal_output = nullptr;
  std::string _caloname = "CEMC";

  bool m_use_vertextype {false};
  GlobalVertex::VTXTYPE m_vertex_type = GlobalVertex::UNDEFINED;

  bool m_UseTowerInfo; // decide whether to use "TowerInfo" (expt_data) or "RawTower" (sim_data)

  TH1F *pairInvMassTotal = NULL;
  TH1F *pairInvMassTotalBkgd = NULL;

	TH3F *pairInvMassPtEta = NULL;
	TH3F *pairInvMassPtEtaBkgd = NULL;

	TH3F *pairpTDelPhiDelEta = NULL;
	TH3F *pairpTDelPhiDelEtaBkgd = NULL;
  
	TH3F *pairInvMassPtdelR = NULL;
  TH3F *pairInvMassPtdelRBkgd = NULL;
	
	TH2F *DelR_pairpT_f = NULL;
	TH2F *DelR_pairpT_b = NULL;
  
	TH2F *h_etaphi_clus = NULL;
	
  TTree *_eventTree = nullptr; // initializing Ttree to store data
  // TTree variables
  int _eventNumber = -1;
  int _nClusters = -1;
  int _triggerVector[64] = {0};
	float _nCentrality = -1.0;
	float _vertex[3] = {0};
  float _clusterIDs[10000] = {0};
  float _clusterEnergies[10000] = {0};
  float _clusterPts[10000] = {0};
  float _clusterEtas[10000] = {0};
  float _clusterPhis[10000] = {0};
	float _clusterChi2[10000] = {0};
  
	int maxTowerEta = -1;
  int maxTowerPhi = -1;

  int _maxTowerEtas[10000] = {0};
  int _maxTowerPhis[10000] = {0};

  float alphaCut = -1.;
	int m_cent_nclus_cut;
  
  TFile * f_temp;
};

#endif  //   CALOCALIBEMC_ETA_H
