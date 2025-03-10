#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TSystem.h"
#include "TVector3.h"

#include <iostream>
#include <cmath>

// Convert pseudorapidity (eta) to polar angle (theta)
double Eta2Theta(double etaVal) { 
  return 2.0 * atan(exp(-etaVal)); 
}

void SuperFastMCTheta() {
  
  const char* outputName = "pi0_decay.root";
  int numPions = 1e7;
  float maxPt = 20.0;
  int NumBins = 500;
  
  // Set the particle type: 0 for pi0, 1 for eta.
  // Change this flag to simulate the decay of the desired meson.
  int particleType = 0; // 0 -> pi0, 1 -> eta
  
  TFile* out = new TFile(outputName, "UPDATE");

  // Functions for generating flat distributions in pt, rapidity, and phi.
  TF1* flatPt = new TF1("flatpt", "1", 0.02, maxPt);
  float MYPI = 3.14159;
  TF1* flat_rapidity = new TF1("flatRap", "1", -1, 1);
  TF1* flat_phi = new TF1("flatPhi", "1", 0, 2 * MYPI);
  TF1* real_theta = new TF1("flatTheta", "sin(x)", 0, MYPI);
  
  // Woods-Saxon function defined but not directly used.
  TF1 woodsSax("woodsSax", "1.0 - 1.0/(1.0 + exp((0.5*x - 0.25)/0.02))", 0.0, maxPt);
  woodsSax.SetNpx(20000);

  // Create histograms.
  TH1F* picorr = new TH1F("picorr", "corrected pi0 spect, according to fit values", NumBins, 0.0, maxPt);
  TH1F* gamcorr = new TH1F("gamcorr", "correct gamma spectrum according to pi0 fit values", NumBins, 0.0, maxPt);
  TH2F* invMass = new TH2F("invMassPt", "", 160, 0.0, 0.8, 4, 2.0, 6.0);
  invMass->Sumw2();
  
  TH1F* pipt = new TH1F(*picorr);
  pipt->Reset();
  
  TH3F* ptp_pt1_dang = new TH3F("ptp_pt1_dang", "", 50, 1.0, 20.0, 50, 1.0, 20.0, 600, 0, 0.600);
  
  TH1F* ratiopt = new TH1F(*pipt);
  ratiopt->SetName("ratiopt");
  std::cout << "processing " << ratiopt->GetName() << std::endl;
    
  TH1F* gampt = new TH1F(*pipt);
  TH1F* gampipt = new TH1F(*pipt);
  TH1F* ratiogampipt = new TH1F(*pipt);

  // Create additional histograms for eta if needed.
  TH1F* hetapt = new TH1F(*pipt);
  TH1F* gametapt = new TH1F(*pipt);
  TH1F* ratiogametapt = new TH1F(*pipt);

  // Set histogram names and titles.
  pipt->SetName("pipt");
  gampt->SetName("gampt");
  gampipt->SetName("gampipt");
  ratiogampipt->SetName("ratiogampipt");

  pipt->SetTitle("Pi0 Pt Dist");
  gampt->SetTitle("All Gamma Pt Dist");
  ratiopt->SetTitle("All Gamma / Pi0 Ratio (Pt)");
  gampipt->SetTitle("Pi0 Gam. Pt Dist");
  ratiogampipt->SetTitle("Gamma_Pi0 / Pi0 Ratio (Pt)");

  TH1F* thetaDist = new TH1F("thetaDist", "distribution of opening angles of the two gamma", 3140, 0, 3.14);
	  		
  // Define the transverse momentum spectral functions.
  TF1* pi0pt = new TF1("pi0pt", "[0]*(x^0.5)/([1]+x)^[2]", 0.001, 10);
  pi0pt->SetParameters(1.0e11, 0.5, 8.0);
  
  TF1* etapt = new TF1("etapt", "[0]*x/([1]+x)^[2]", 0.001, 10);
  etapt->SetParameters(1.0e8, 1.1, 5.0);

  // Define the resolution and smearing functions.
  TF1* f1 = new TF1("f1reso", "sqrt([0]*[0]/x+[1]*[1])", 0.5, 20.);
  f1->SetParameter(0, 0.165);
  f1->SetParameter(1, 0.26);
	  
  TF1* fgaus = new TF1("mygaus", "gaus", -1.5, 1.5);
  fgaus->SetParameter(0, 1e6);
  fgaus->SetParameter(1, 0);

  TH2F* h2esm = new TH2F("h2esm", "", 80, 0, 20, 200, 0, 20);
  TH2F* h2resoe = new TH2F("h2resoe", "", 80, 0, 20, 200, -1, 1);

  // Main Monte Carlo loop.
  for (int i = 0; i < numPions; i++) {
    // Select the mass based on the chosen particle type.
    float massPi = (particleType == 0) ? 0.13497 : 0.547;
    
    float thetaOrig = Eta2Theta(flat_rapidity->GetRandom());
    float ptPi = flatPt->GetRandom();
    float pPi = ptPi / sin(thetaOrig);
    float ePi = sqrt(pPi * pPi + massPi * massPi);

    if (ePi < 100) {
      // Build the pion (or eta) 3-vector and 4-vector.
      TVector3 pion3v(2.0, 0.0, 0.0);
      pion3v.SetTheta(thetaOrig);
      pion3v.SetPhi(flat_phi->GetRandom());
      pion3v.SetMag(pPi);
		
      TLorentzVector pcylorig;
      pcylorig.SetVect(pion3v);
      pcylorig.SetE(ePi);
		
      // Simulate the two-photon decay in the meson's rest frame.
      TLorentzVector photon1;
      TLorentzVector photon2;
      TLorentzVector* g1 = &photon1;
      TLorentzVector* g2 = &photon2;

      TVector3 photon3v(2.0, 0.0, 0.0);
      photon3v.SetMag(massPi / 2);
      photon3v.SetTheta(real_theta->GetRandom());
      photon3v.SetPhi(flat_phi->GetRandom());
      g1->SetVect(photon3v);
      g1->SetE(massPi / 2);
		
      g2->SetVect(-photon3v);
      g2->SetE(massPi / 2);
		
      // Boost the photons from the meson's rest frame to the lab frame.
      g1->Boost(pcylorig.BoostVector());
      g2->Boost(pcylorig.BoostVector());
		
      TLorentzVector pho1, pho2, pi0lv;

      float tt_clus_energy = g1->E();
      float reso1 = f1->Eval(tt_clus_energy);
      fgaus->SetParameter(2, reso1);
      float smear1 = fgaus->GetRandom();
      tt_clus_energy += tt_clus_energy * smear1;
      h2esm->Fill(g1->E(), tt_clus_energy);
      h2resoe->Fill(g1->E(), tt_clus_energy / g1->E() - 1.0);
      float tt_clus_pt = g1->Pt();
      tt_clus_pt += tt_clus_pt * smear1;
      float tt_clus_phi = g1->Phi();
      float tt_clus_eta = g1->Eta();

      float tt_clus_energy2 = g2->E();
      float reso2 = f1->Eval(tt_clus_energy2);
      fgaus->SetParameter(2, reso2);
      float smear2 = fgaus->GetRandom();
      tt_clus_energy2 += tt_clus_energy2 * smear2;
      h2esm->Fill(g2->E(), tt_clus_energy2);
      h2resoe->Fill(g2->E(), tt_clus_energy2 / g2->E() - 1.0);
      float tt_clus_pt2 = g2->Pt();
      tt_clus_pt2 += tt_clus_pt2 * smear2;
      float tt_clus_phi2 = g2->Phi();
      float tt_clus_eta2 = g2->Eta();

      pho1.SetPtEtaPhiE(tt_clus_pt, tt_clus_eta, tt_clus_phi, tt_clus_energy);
      pho2.SetPtEtaPhiE(tt_clus_pt2, tt_clus_eta2, tt_clus_phi2, tt_clus_energy2);
		
      // Reconstruct the invariant mass.
      pi0lv = pho1 + pho2;
      float pairInvMass = pi0lv.M();

      float weight = 1.0;
      // Use the appropriate spectral function based on the particle type.
      if (particleType == 0) {
        weight = pi0pt->Eval(ptPi);
        pipt->Fill(ptPi, weight);
      } else {
        weight = etapt->Eval(ptPi);
        hetapt->Fill(ptPi, weight);
      }
		
      float alph = fabs(pho1.E() - pho2.E()) / (pho1.E() + pho2.E());
      if (((pho1.E() > 1.0 && pho2.E() > 0.60) ||
           (pho2.E() > 1.0 && pho1.E() > 0.60)) &&
          fabs(pi0lv.Pt()) > 2.00 && alph < 0.6) {
        invMass->Fill(pairInvMass, pi0lv.Pt(), weight);
      }

      if (pho2.E() > 0.60 && pho1.E() > 0.6 && alph < 1.0) {
        ptp_pt1_dang->Fill(pi0lv.Pt(), pho1.Pt(), pho1.DeltaR(pho2));
      }
		    
      gampt->Fill(g1->Pt(), weight);
      gampt->Fill(g2->Pt(), weight);
      gamcorr->Fill(g1->Pt(), weight);
      gamcorr->Fill(g2->Pt(), weight);
		
      thetaDist->Fill(g1->Vect().Angle(g2->Vect()));
      if (particleType == 0) {
        gampipt->Fill(g1->Pt(), weight);
        gampipt->Fill(g2->Pt(), weight);
      } else {
        gametapt->Fill(g1->Pt(), weight);
        gametapt->Fill(g2->Pt(), weight);
      }
    }
  }
	  
  // Write the histograms and functions to file.
  invMass->Write();
  ptp_pt1_dang->Write();
  pi0pt->Write();
  etapt->Write();

  for (int k = 1; k <= NumBins; k++) {
    float binCent = picorr->GetBinCenter(k);
    float ptval = pi0pt->Eval(binCent);
    picorr->SetBinContent(k, ptval);
  }
  
  picorr->Scale(1e-5);
  gamcorr->Scale(1e-15);

  ratiogampipt->Divide(gampipt, pipt);
  ratiopt->Divide(gampt, pipt, 1.02, 1.0);
	  
  out->Write(0, TObject::kOverwrite);
  std::cout << "Done" << std::endl;
  out->Close();
	gSystem->Exit(0);
}

