// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.cc
// 
/*

 Description: an example analyzer using SUSY Ntuples

 Implementation:
 The macro is driven by ana.C in this directory.
 Fills MET and diEMPt histograms for gg and ff events with >= 1 good jets. gg and ff
 skimming can be done simultaneously by setting CopyEvents to true in ana.C.

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.cc,v 1.21 2013/05/11 14:14:45 yiiyama Exp $
//

#include <TH1F.h>
#include <TFile.h>

#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <map>
#include <utility>

// REMOVE THE LINES BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h" // to access the JEC scales
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" // to access the uncertainties

#include "SusyEventAnalyzer.h"

template<typename T>
bool
PtGreater(const T* p1, const T* p2) {
  return p1->momentum.Pt() > p2->momentum.Pt();
}

template<typename T1, typename T2>
bool
isSameObject(const T1& p1, const T2& p2, double dRCut = 0.1)
{
  float dEta = p1.momentum.Eta() - p2.momentum.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.momentum.Phi() - p2.momentum.Phi());
  return dEta*dEta + dPhi*dPhi < dRCut * dRCut;
}

void
photonEffectiveAreas(double _eta, double* _effA)
{
  double& effACH(_effA[0]);
  double& effANH(_effA[1]);
  double& effAPh(_effA[2]);

  // Source: CutBasedPhotonID2012 twiki
  if(_eta < 1.){
    effACH = 0.012;
    effANH = 0.03;
    effAPh = 0.148;
  }
  else if(_eta < 1.479){
    effACH = 0.010;
    effANH = 0.057;
    effAPh = 0.13;
  }
  else if(_eta < 2.){
    effACH = 0.014;
    effANH = 0.039;
    effAPh = 0.112;
  }
  else if(_eta < 2.2){
    effACH = 0.012;
    effANH = 0.015;
    effAPh = 0.216;
  }
  else if(_eta < 2.3){
    effACH = 0.016;
    effANH = 0.024;
    effAPh = 0.262;
  }
  else if(_eta < 2.4){
    effACH = 0.02;
    effANH = 0.039;
    effAPh = 0.26;
  }
  else{
    effACH = 0.012;
    effANH = 0.072;
    effAPh = 0.266;
  }
}

void
electronEffectiveAreas(double _eta, double &_effA)
{
  // Source: EgammaEARhoCorrection twiki
  if(_eta < 1.)
    _effA = 0.13;
  else if(_eta < 1.479)
    _effA = 0.14;
  else if(_eta < 2.)
    _effA = 0.07;
  else if(_eta < 2.2)
    _effA = 0.09;
  else if(_eta < 2.3)
    _effA = 0.11;
  else if(_eta < 2.4)
    _effA = 0.11;
  else
    _effA = 0.14;
}

////////// MAIN ANALYSIS FUNCTION //////////
void
SusyEventAnalyzer::Run()
{
  int nCnt[20];
  for(int i(0); i<20; i++) nCnt[i] = 0;

  ////////// TEXT OUTPUT //////////

  std::ofstream outFile;
  if(logFileName != "cout"){
    outFile.open(logFileName);
    if(!outFile.is_open()){
      std::cerr << "Log output " << logFileName << " could not be opened." << std::endl;
      return;
    }
  }
  std::ostream& out(logFileName == "cout" ? std::cout : outFile);

  ////////// SWITCH FOR USING USER-PROVIDED JEC //////////
  bool useCustomJEC(false);

  // REMOVE THE LINES BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
  FactorizedJetCorrector* jetCorrection(0);
  JetCorrectionUncertainty* jecUncertainty(0);

  TFile* ggFile(0);
  TFile* ffFile(0);
  TFile* fout(0);

  try{

    ////////// INITIALIZE SKIMMED CLONE //////////
    TTree* ggTree(0);
    TTree* ffTree(0);
    if(copyEvents){

      if(printLevel > 0) out << "Open file for skim output" << std::endl;

      ggFile = TFile::Open("susyEvents_" + outputName + "_gg.root", "RECREATE");
      ffFile = TFile::Open("susyEvents_" + outputName + "_ff.root", "RECREATE");
      if(!ggFile || ggFile->IsZombie() || !ffFile || ffFile->IsZombie()){
        std::cerr << "Cannot open output file susyEvents_" << outputName << ".root" << std::endl;
        throw std::runtime_error("IOError");
      }
      else{
        ggFile->cd();
        ggTree = new TTree("susyTree", "SUSY Event");
        ggTree->SetAutoSave(10000000);

        ffFile->cd();
        ffTree = new TTree("susyTree", "SUSY Event");
        ffTree->SetAutoSave(10000000);

        /* Register the output trees to the Event object - the branches will be booked internally */
        event.addOutput(*ggTree);
        event.addOutput(*ffTree);
      }

    }

    ////////// INITIALIZE HISTOGRAM OUTPUT //////////
    if(printLevel > 0) out << "Open file for histograms" << std::endl;

    fout = TFile::Open("hist_" + outputName + ".root", "RECREATE");
    if(!fout || fout->IsZombie()){
      std::cerr << "Cannot open output file hist_" << outputName << ".root" << std::endl;
      throw std::runtime_error("IOError");
    }

    if(printLevel > 0) out << "Define histograms" << std::endl;

    TH1F* h_met_gg(new TH1F("h_met_gg","#gamma#gamma #slash{E}_{T};#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
    TH1F* h_met_ff(new TH1F("h_met_ff","ff #slash{E}_{T};#slash{E}_{T} (GeV);Events / GeV", 200, 0., 200.));
    TH1F* h_diEMPt_gg(new TH1F("h_diEMPt_gg", "#gamma#gamma diEM P_{T};diEM P_{T} (GeV);Events / GeV", 200, 0., 200.));
    TH1F* h_diEMPt_ff(new TH1F("h_diEMPt_ff", "ff diEM P_{T};diEM P_{T} (GeV);Events / GeV", 200, 0., 200.));
    TH1F* h_PFDiEMPt_gg(new TH1F("h_PFDiEMPt_gg", "#gamma#gamma PF diEM P_{T};diEM P_{T} (GeV);Events / GeV", 200, 0., 200.));
    TH1F* h_PFDiEMPt_ff(new TH1F("h_PFDiEMPt_ff", "ff PF diEM P_{T};diEM P_{T} (GeV);Events / GeV", 200, 0., 200.));

    ////////// INITIALIZE JEC //////////
    // REMOVE THE LINES BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
    if(useCustomJEC){
      if(printLevel > 0) out << "Initialize jet energy corrections" << std::endl;

      std::string jecSourcePrefix("../jec/FT_53_V21_AN3_");
      jetCorrection = new FactorizedJetCorrector("L1FastJet:L2Relative:L3Absolute:L2Relative",
                                       jecSourcePrefix + "L1FastJet_AK5PF.txt:" +
                                       jecSourcePrefix + "L2Relative_AK5PF.txt:" +
                                       jecSourcePrefix + "L3Absolute_AK5PF.txt:" +
                                       jecSourcePrefix + "L2RelativeL3AbsoluteResidual_AK5PF.txt");
      jecUncertainty = new JetCorrectionUncertainty(jecSourcePrefix + "Uncertainty_AK5PF.txt");
    }

    /////////////////////////////////////
    ////////// MAIN EVENT LOOP //////////
    /////////////////////////////////////

    if(printLevel > 0) out << "Start event loop" << std::endl;

    long iEntry(0);
    while(iEntry != processNEvents && event.getEntry(iEntry++) != 0){

      if(printLevel > 1 || iEntry % printInterval == 0)
        out << "Begin processing event " << iEntry - 1 << ". Current: run=" << event.runNumber << ", event=" << event.eventNumber << std::endl;

      ////////// DUMP EVENT CONTENT //////////
      if(printLevel > 2) event.Print(out);

      /* total number of events */
      nCnt[0]++;

      ////////// GOOD LUMI FILTER //////////
      if(printLevel > 1) out << "Apply good lumi list." << std::endl;
      if(!IsGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;

      ////////// MET FILTER //////////
      if(event.isRealData){
        if(printLevel > 1) out << "Apply MET filter." << std::endl;
        if(!event.passMetFilters()) continue;
      }

      ////////// TRIGGER //////////
      if(printLevel > 1) out << "Apply HLT cut." << std::endl;
      if(!PassTriggers()) continue;

      ////////// REQUIRE AT LEAST ONE GOOD VERTEX //////////
      if(printLevel > 1) out << "Require at least one good vertex." << std::endl;
      unsigned nV(event.vertices.size());
      unsigned iV(0);
      for(; iV != nV; ++iV){
        susy::Vertex& vertex(event.vertices[iV]);
        if(vertex.ndof >= 4 && std::abs(vertex.position.Z()) < 24. && vertex.position.Perp() < 2.) break;
      }
      if(iV == nV) continue;

      susy::Vertex const& primVtx(event.vertices[iV]);

      if(printLevel > 1) out << "Event passes preselection." << std::endl;

      /* number of events passing preselection */
      nCnt[1]++;

      ////////// SELECT GOOD AND FAKE PHOTONS //////////
      if(printLevel > 1) out << "Find good barrel photons in the event" << std::endl;

      std::vector<susy::Photon const*> goodPhotons;
      std::vector<susy::Photon const*> fakePhotons;

      susy::PhotonCollection const& photons(event.photons["photons"]);
      unsigned nPhoton(photons.size());
      for(unsigned iP(0); iP != nPhoton; ++iP){
        susy::Photon const& photon(photons[iP]);

        double pt(photon.momentum.Pt());
        if(pt < 25.) continue;

        double absEta(std::abs(photon.momentum.Eta()));
        if(absEta > susy::etaGapBegin) continue; /* use barrel photons only */

        if(photon.nPixelSeeds > 0) continue;

        if(photon.hadTowOverEm > 0.05) continue;

        double effA[3];
        photonEffectiveAreas(absEta, effA);

        if(photon.neutralHadronIso - event.rho25 * effA[1] - 0.04 * pt > 3.5) continue;

        if(photon.photonIso - event.rho25 * effA[2] - 0.005 * pt > 1.3) continue;

        double chIso(photon.chargedHadronIso - event.rho25 * effA[0]);

        if(photon.sigmaIetaIeta < 0.012 && chIso < 2.6) goodPhotons.push_back(&photon);
        else if(photon.sigmaIetaIeta < 0.014 && chIso < 15.) fakePhotons.push_back(&photon);
      }

      ////////// SORT SELECTED PHOTONS //////////
      std::sort(goodPhotons.begin(), goodPhotons.end(), PtGreater<susy::Photon>);
      std::sort(fakePhotons.begin(), fakePhotons.end(), PtGreater<susy::Photon>);

      /* leading photon has to have Pt greater than 40 GeV */
      if(goodPhotons.size() != 0 && goodPhotons[0]->momentum.Pt() < 40.) goodPhotons.clear();
      if(fakePhotons.size() != 0 && fakePhotons[0]->momentum.Pt() < 40.) fakePhotons.clear();

      if(goodPhotons.size() < 2 && fakePhotons.size() < 2) continue;
      if(goodPhotons.size() >= 2 && fakePhotons.size() >= 2){
        out << "Run " << event.runNumber << " Lumi " << event.luminosityBlockNumber << " Event " << event.eventNumber;
        out << " has >= 2 good photons AND >= 2 fake photons. Check for pathologies!" << std::endl;
      }

      /* number of events with two good photons */
      if(goodPhotons.size() >= 2) nCnt[2]++;
      /* number of events with two fake photons */
      if(fakePhotons.size() >= 2) nCnt[3]++;


      ////////// VETO EVENTS WITH LOOSE MUONS //////////
      if(printLevel > 1) out << "Find loose muons in the event" << std::endl;

      susy::MuonCollection const& muons(event.muons["muons"]);
      unsigned nMuon(muons.size());
      unsigned iMuon(0);
      for(; iMuon != nMuon; ++iMuon){
        susy::Muon const& muon(muons[iMuon]);

        if(muon.momentum.Pt() < 10.) continue;

        if(muon.isPFMuon() && (muon.isTrackerMuon() || muon.isGlobalMuon())) break;
      }
      if(iMuon != nMuon) continue;

      /* number of events with two good photons and no muon */
      if(goodPhotons.size() >= 2) nCnt[4]++;
      /* number of events with two fake photons and no muon */
      if(fakePhotons.size() >= 2) nCnt[5]++;


      ////////// VETO EVENTS WITH LOOSE ELECTRONS //////////
      if(printLevel > 1) out << "Find loose electrons in the event" << std::endl;

      susy::ElectronCollection const& electrons(event.electrons["gsfElectrons"]);
      unsigned nEle(electrons.size());
      unsigned iEle(0);
      for(; iEle != nEle; ++iEle){
        susy::Electron const& electron(electrons[iEle]);

        double pt(electron.momentum.Pt());
        if(pt < 10.) continue;

        double absEta(std::abs(electron.momentum.Eta()));
        bool isBarrel(absEta < susy::etaGapBegin);
        if((!isBarrel && absEta < susy::etaGapEnd) || absEta > susy::etaMax) continue;

        unsigned iP(0);
        for(; iP != goodPhotons.size(); ++iP)
          if(isSameObject(electron, *goodPhotons[iP])) break;
        if(iP != goodPhotons.size()) continue;

        if(std::abs(electron.deltaEtaSuperClusterTrackAtVtx) > isBarrel ? 0.007 : 0.01) continue;
        if(std::abs(electron.deltaPhiSuperClusterTrackAtVtx) > isBarrel ? 0.8 : 0.7) continue;
        if(electron.sigmaIetaIeta > isBarrel ? 0.01 : 0.03) continue;
        if(isBarrel && electron.hcalOverEcalBc > 0.15) continue;

        if(electron.gsfTrack == 0) continue;
        if(std::abs(electron.gsfTrack->dxy(primVtx.position)) > 0.04) continue;
        if(std::abs(electron.gsfTrack->dz(primVtx.position)) > 0.2) continue;

        double effA;
        electronEffectiveAreas(absEta, effA);
        if((electron.chargedHadronIso + electron.neutralHadronIso + electron.photonIso - event.rho25 * effA) / pt > 0.15) continue;

        /* this is a loose electron (working point "Veto") */
        break;
      }
      if(iEle != nEle) continue;

      /* number of events with two good photons and no lepton */
      if(goodPhotons.size() >= 2) nCnt[6]++;
      /* number of events with two fake photons and no lepton */
      if(fakePhotons.size() >= 2) nCnt[7]++;

      ////////// SELECT JETS //////////
      if(printLevel > 1) out << "Find good PFJets in the event" << std::endl;

      susy::PFJetCollection const& pfJets(event.pfJets["ak5"]);
      unsigned nJets(pfJets.size());
      unsigned nGoodJets(0);
      for(unsigned iJ(0); iJ != nJets; ++iJ){
        susy::PFJet const& jet(pfJets[iJ]);

        double eta(jet.momentum.Eta());
        if(std::abs(eta) > 2.6) continue;

        if(jet.nConstituents < 2) continue;
        if(jet.chargedMultiplicity == 0) continue;

        double energy(jet.momentum.E());

        if(jet.neutralHadronEnergy / energy > 0.99) continue;
        if(jet.neutralEmEnergy / energy > 0.99) continue;
        if(jet.chargedHadronEnergy / energy < 1.e-6) continue;
        if(jet.chargedEmEnergy / energy > 0.99) continue;

        unsigned iP(0);
        for(; iP != goodPhotons.size(); ++iP)
          if(isSameObject(jet, *goodPhotons[iP])) break;
        if(iP != goodPhotons.size()) continue;

        double pt(jet.momentum.Pt());

        float jecScale;
        float jecScaleUncertainty;
        if(useCustomJEC){
          jetCorrection->setJetEta(eta);
          jetCorrection->setJetPt(pt);
          jetCorrection->setJetA(jet.jetArea);
          jetCorrection->setRho(event.rho);

          jecScale = jetCorrection->getCorrection();

          jecUncertainty->setJetEta(eta);
          jecUncertainty->setJetPt(pt);

          try{
            jecScaleUncertainty = jecUncertainty->getUncertainty(true);
          }
          catch(std::exception& e){
            std::cerr << "Cannot get uncertainty for jet Pt = " << pt << " Eta = " << eta << ". Setting to -1." << std::endl;
            jecScaleUncertainty = -1.;
          }
        }
        else{
          // "Residual" correction for data is already included in the ntuplizer (see line 2279 of SusyNtuplizer.cc)
          jecScale = jet.jecScaleFactors.find("L1FastL2L3")->second;

          jecScaleUncertainty = jet.jecUncertainty;
        }

        if(jecScaleUncertainty > 0.2) continue;

        TLorentzVector corrP4(jecScale * jet.momentum);

        if(corrP4.Pt() < 30.) continue;

        ++nGoodJets;
      }

      if(nGoodJets > 0){
        /* number of events with two good photons, no lepton, and >=1 good jets */
        if(goodPhotons.size() >= 2) nCnt[8]++;
        /* number of events with two fake photons, no lepton, and >=1 good jets */
        if(fakePhotons.size() >= 2) nCnt[9]++;
      }

      ////////// MET //////////
      TVector2 const& metV(event.metMap["pfType01CorrectedMet"].mEt);

      if(nGoodJets > 0 && metV.Mod() > 50.){
        /* number of events with two good photons, no lepton, >=1 good jets, and MET > 50 GeV */
        if(goodPhotons.size() >= 2) nCnt[10]++;
        /* number of events with two fake photons, no lepton, >=1 good jets, and MET > 50 GeV */
        if(fakePhotons.size() >= 2) nCnt[11]++;
      }

      ////////// PF PARTICLES (DEMONSTRATION OF THE WORKAROUND FOR THE BUG IN TAG cms538v0 / cms538v1) ///////////
      std::map<std::pair<int, float>, susy::PFParticle const*> uniquePF;

      unsigned nPF(event.pfParticles.size());
      for(unsigned iPF(0); iPF != nPF; ++iPF){
        susy::PFParticle const& particle(event.pfParticles[iPF]);
        uniquePF[std::pair<int, float>(particle.pdgId, particle.momentum.Pt())] = &particle;
      }

      std::vector<susy::PFParticle const*> pfParticles;
      for(std::map<std::pair<int, float>, susy::PFParticle const*>::iterator pfItr(uniquePF.begin()); pfItr != uniquePF.end(); ++pfItr)
        pfParticles.push_back(pfItr->second);

      nPF = pfParticles.size();

      ////////// FILL HISTOGRAMS //////////

      if(nGoodJets > 0){
        if(goodPhotons.size() >= 2){
          h_met_gg->Fill(metV.Mod());
          TLorentzVector diEMP(goodPhotons[0]->momentum);
          diEMP += goodPhotons[1]->momentum;
          h_diEMPt_gg->Fill(diEMP.Pt());

          TLorentzVector PFDiEMP;
          for(unsigned iPF(0); iPF != nPF; ++iPF)
            if(pfParticles[iPF]->momentum.DeltaR(goodPhotons[0]->momentum) < 0.3 || pfParticles[iPF]->momentum.DeltaR(goodPhotons[1]->momentum) < 0.3)
              PFDiEMP += pfParticles[iPF]->momentum;
          h_PFDiEMPt_gg->Fill(PFDiEMP.Pt());
        }
        if(fakePhotons.size() >= 2){
          h_met_ff->Fill(metV.Mod());
          TLorentzVector diEMP(fakePhotons[0]->momentum);
          diEMP += fakePhotons[1]->momentum;
          h_diEMPt_ff->Fill(diEMP.Pt());

          TLorentzVector PFDiEMP;
          for(unsigned iPF(0); iPF != nPF; ++iPF)
            if(pfParticles[iPF]->momentum.DeltaR(fakePhotons[0]->momentum) < 0.3 || pfParticles[iPF]->momentum.DeltaR(fakePhotons[1]->momentum) < 0.3)
              PFDiEMP += pfParticles[iPF]->momentum;
          h_PFDiEMPt_ff->Fill(PFDiEMP.Pt());
        }
      }

      ////////// FILL SKIMS //////////
      if(copyEvents){
        if(goodPhotons.size() >= 2) ggTree->Fill();
        if(fakePhotons.size() >= 2) ffTree->Fill();
      }

    }

    ////////// END OF EVENT LOOP //////////

    out << " ----------------- Job Summary ----------------- " << std::endl;
    out << " Total events                                                 : " << nCnt[0] << std::endl;
    if(nCnt[0] > 0){
      out << " passed preselection                                          : " << nCnt[1] << " (" << nCnt[1]/float(nCnt[0]) << ")" << std::endl;
      out << " goodPhotons >= 2                                             : " << nCnt[2] << " (" << nCnt[2]/float(nCnt[1]) << ")" << std::endl;
      out << " fakePhotons >= 2                                             : " << nCnt[3] << " (" << nCnt[3]/float(nCnt[1]) << ")" << std::endl;
      out << " goodPhotons >= 2 && 0 muon                                   : " << nCnt[4] << " (" << nCnt[4]/float(nCnt[1]) << ")" << std::endl;
      out << " fakePhotons >= 2 && 0 muon                                   : " << nCnt[5] << " (" << nCnt[5]/float(nCnt[1]) << ")" << std::endl;
      out << " goodPhotons >= 2 && 0 lepton                                 : " << nCnt[6] << " (" << nCnt[6]/float(nCnt[1]) << ")" << std::endl;
      out << " fakePhotons >= 2 && 0 lepton                                 : " << nCnt[7] << " (" << nCnt[7]/float(nCnt[1]) << ")" << std::endl;
      out << " goodPhotons >= 2 && 0 lepton && >= 1 good jet                : " << nCnt[8] << " (" << nCnt[8]/float(nCnt[1]) << ")" << std::endl;
      out << " goodPhotons >= 2 && 0 lepton && >= 1 good jet                : " << nCnt[9] << " (" << nCnt[9]/float(nCnt[1]) << ")" << std::endl;
      out << " goodPhotons >= 2 && 0 lepton && >= 1 good jet && MET > 50 GeV: " << nCnt[10] << " (" << nCnt[10]/float(nCnt[1]) << ")" << std::endl;
      out << " goodPhotons >= 2 && 0 lepton && >= 1 good jet && MET > 50 GeV: " << nCnt[11] << " (" << nCnt[11]/float(nCnt[1]) << ")" << std::endl;
    }

    if(printLevel > 0) out << "Save outputs" << std::endl;

    ////////// END OF EVENT LOOP //////////

    fout->cd();
    fout->Write();

    if(copyEvents){
      event.releaseTree(*ggTree);
      event.releaseTree(*ffTree);

      ggFile->cd();
      ggFile->Write();

      ffFile->cd();
      ffFile->Write();
    }
  }
  catch(std::exception& e){
    std::cerr << e.what() << std::endl;

    event.releaseTrees();
  }

  delete jetCorrection;
  delete jecUncertainty;

  delete ggFile;
  delete ffFile;
  delete fout;
}

