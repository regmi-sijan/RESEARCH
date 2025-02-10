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
  
  void SaveAllMode(bool state) // to decide if we want to save information of all secondary particle (default is false)
  {
    AllMode=state;
    return;
  }
 
protected:
  bool AllMode{false}; // "false" = save only final secondary particles
 
private:

  float getEta(PHG4Particle *particle);
    
  int n_event;
	
  TTree *_eventTree = nullptr;
	
  int _eventNumber = -1;
  int _nFourVector = -1; // store the number of clusters (truth info) recorded in each event (our focus will only be in EMCal)
  float _fv_energy[10000] = {0}; // truth energy of four vector
  float _fv_px[10000] = {0}; // 3 momentum (px, py and pz) of the cluster (truth particle)
  float _fv_py[10000] = {0};
  float _fv_pz[10000] = {0};
  float _fv_Eta[10000] = {0}; // pseudorapidity of the four vector
  float _pid_particle[10000]; // pid of associated particle
  float _primary_id[10000]; // primary id of the particle
  float _parent_id[10000]; // parent id of the particle
  float _embedding[10000]; // parent id of the particle
  float _track_id[10000]; // track id of particle
  

  TFile *out;
	
  std::string Outfile;

};

#endif // PI0CLUSTERANA_H
