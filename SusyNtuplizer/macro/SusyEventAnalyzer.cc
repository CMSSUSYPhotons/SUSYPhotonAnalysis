// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.cc
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.cc,v 1.4 2011/04/19 20:15:20 dwjang Exp $
//

#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>

#include "SusyEventAnalyzer.h"
#include "SusyEventPrinter.h"

template<typename T> bool EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}


void SusyEventAnalyzer::InitializePerEvent() {

}


bool SusyEventAnalyzer::isSameObject(TLorentzVector& p1, TLorentzVector& p2) {

  float dEta = p1.Eta() - p2.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
  float dR = std::sqrt(dEta*dEta + dPhi*dPhi);
  if(dR < 0.5) return true;
  return false;
}


float SusyEventAnalyzer::d0correction(TVector3& beamSpot, susy::Track& track) const {

  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}


bool SusyEventAnalyzer::PassTrigger(TString path) {
  bool pass = false;
  for(susy::TriggerMap::iterator it = event->hltMap.begin(); it != event->hltMap.end(); it++) {
    if(it->first.Contains(path) && (int(it->second.second)) ) {
      pass = true;
      break;
    }
  }
  return pass;
}


bool SusyEventAnalyzer::PassTriggers() {
  bool pass = false;
  for(std::vector<TString>::iterator it = hltNames.begin(); it != hltNames.end(); it++) {
    if(PassTrigger(*it)) {
      pass = true;
      break;
    }
  }
  return pass;
}



void SusyEventAnalyzer::Loop() {

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();

  std::cout << "total events in files  : " << nentries << std::endl;

  if(processNEvents <= 0 || processNEvents > nentries) processNEvents = nentries;

  std::cout << "events to be processed : " << processNEvents << std::endl; 


  if(printLevel > 0) std::cout << "Initialize event counters." << std::endl;
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  int nFiltered = 0;
  TTree* filterTree = 0;

  if(enableFilter) {
    TFile* filterFile = new TFile(filtered_file_name,"RECREATE");
    filterTree = (TTree*) fChain->GetTree()->CloneTree(0);
    filterTree->SetAutoSave();
  }


  // open hist file and define histograms

  TFile* fout = new TFile("hist_"+ds+".root","RECREATE");

  fout->cd();

  TH1F* h_vtxZ = new TH1F("vtxZ","Z position of the primary vertex;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_bsZ = new TH1F("bsZ","Z position of the beam spot;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_met = new TH1F("met","missing transverse energy;#slash{E}_{T} (GeV);Events",200,0.0,1000.0);
  TH1F* h_sumEt = new TH1F("sumEt","Scalar sum of all calorimeter energy;#sigmaE_{T} (GeV);Events",200,0.0,2000.0);



  // to check duplicated events
  std::map<int, std::set<int> > allEvents;

  // start event looping

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry < processNEvents; jentry++) {

    if(printLevel > 0) std::cout << "Get the tree contents." << std::endl;

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      std::cout << int(jentry) << " events processed with run="
		<< event->runNumber << ", event=" << event->eventNumber << std::endl;
    }


    if(printLevel > 0) std::cout << "Initialize any global variables to be reset per event." << std::endl;

    InitializePerEvent();


    if(printLevel > 0) std::cout << "Apply good run list." << std::endl;

    
    // uncomment this to use the Json file to flag good data (or bad depending on your outlook)    
    // if(!isInJson(event->runNumber,event->luminosityBlockNumber)) continue;

    Print(*event);


    if(printLevel > 0) std::cout << "Check duplicated events for data only." << std::endl;

    bool duplicateEvent = ! (allEvents[event->runNumber].insert(event->eventNumber)).second;
    if(event->isRealData && duplicateEvent) continue;
 

    if(printLevel > 0) std::cout << "Setup object vectors." << std::endl;

    // loose + tight cuts come from syncronization exercise of Susy RA4 group

    // loose objects have all standard cuts except for isolation
    std::vector<susy::Muon*>     loose_muons;
    std::vector<susy::Electron*> loose_electrons;
    std::vector<susy::Photon*>   loose_photons;
    std::vector<susy::CaloJet*>  loose_jets;

    // tight objects hava isolation cuts applied on top of loose objects
    std::vector<susy::Muon*>     tight_muons;
    std::vector<susy::Electron*> tight_electrons;
    std::vector<susy::Photon*>   tight_photons;
    std::vector<susy::CaloJet*>  tight_jets;


    if(printLevel > 0) std::cout << "Find primary vertex in the event." << std::endl;

    TVector3* primVtx = 0;
    if(event->vertices.size() > 0) primVtx = &(event->vertices[0]);

    if(primVtx) h_vtxZ->Fill(primVtx->Z());
    h_bsZ->Fill(event->beamSpot.Z());


    if(printLevel > 0) std::cout << "Find loose and tight muons in the event." << std::endl;

    for(std::vector<susy::Muon>::iterator it = event->muons.begin();
	it != event->muons.end(); it++) {

      float pt = it->momentum.Pt();
      if(pt < 10) continue;

      if(std::abs(it->momentum.Eta()) > 2.1) continue;

      if(it->combinedTrackIndex < 0) continue;

      susy::Track& combinedTrack = event->tracks[it->combinedTrackIndex];

      float normChi2 = combinedTrack.chi2 / combinedTrack.ndof;
      if(normChi2 > 10) continue;

      float d0 = d0correction(event->beamSpot,combinedTrack);
      if(std::abs(d0) > 0.2) continue;

      if(it->nValidHits < 11) continue;

      loose_muons.push_back(&*it);

      std::map<TString,UChar_t>::iterator idPair = it->idPairs.find("GlobalMuonPromptTight");
      if(idPair == it->idPairs.end()) continue;
      if(int(idPair->second) == 0) continue;

      float relIso = (it->ecalIsoR03 + it->hcalIsoR03 + it->trackIsoR03)/pt;
      if(relIso > 0.1) continue;

      tight_muons.push_back(&*it);

    } // for muon

    if(printLevel > 0) std::cout << "Find loose and tight electrons in the event." << std::endl;

    std::map<TString, std::vector<susy::Electron> >::iterator eleMap = event->electrons.find("gsfElectrons");
    if(eleMap == event->electrons.end()) {
      if(event->electrons.size() > 0) std::cout << "electron collection is not available!!!" << std::endl;
    }
    else {
      for(std::vector<susy::Electron>::iterator it = eleMap->second.begin();
	  it != eleMap->second.end(); it++) {

	float pt = it->momentum.Pt();
	if(pt < 20) continue;

	float absEta = std::abs(it->momentum.Eta());
	if(absEta > 2.5) continue;
	if(absEta > 1.47 && absEta < 1.567) continue;

	if(it->gsfTrackIndex < 0) continue;

	susy::Track& gsfTrack = event->tracks[it->gsfTrackIndex];
	float d0 = d0correction(event->beamSpot, gsfTrack);
	//	if(std::abs(d0) > 0.2) continue;

	bool same = false;
	for(std::vector<susy::Muon*>::iterator m_it = tight_muons.begin();
	    m_it != tight_muons.end(); m_it++){
	  if(isSameObject(it->momentum,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}

	if(same) continue;

	loose_electrons.push_back(&*it);

	std::map<TString,Float_t>::iterator idPair = it->idPairs.find("eidRobustTight");
	if(idPair == it->idPairs.end()) continue;
	if(idPair->second < 0.5) continue;

	float relIso = (it->dr04EcalRecHitSumEt + it->dr04HcalTowerSumEt())/pt;
	if(relIso > 0.1) continue;

	tight_electrons.push_back(&*it);

      }// for electron
    }// else 


    if(printLevel > 0) std::cout << "Find loose and tight photons in the event." << std::endl;

    std::map<TString, std::vector<susy::Photon> >::iterator phoMap = event->photons.find("photons");
    if(phoMap == event->photons.end()) {
      if(event->photons.size() > 0) std::cout << "photon collection is not available!" << std::endl;
    }
    else {
      for(std::vector<susy::Photon>::iterator it = phoMap->second.begin();
	  it != phoMap->second.end(); it++) {

	// fiducial cuts
	if(std::abs(it->caloPosition.Eta()) > 2.5) continue;
	if(it->isEBEtaGap() || it->isEBPhiGap() || it->isEERingGap() || it->isEEDeeGap() || it->isEBEEGap()) continue;

	// Et cuts, for the time being, 20 GeV
	if(it->momentum.Et() < 20.0) continue;

	bool same = false;
	for(std::vector<susy::Muon*>::iterator m_it = tight_muons.begin();
	    m_it != tight_muons.end(); m_it++){
	  if(isSameObject(it->momentum,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	for(std::vector<susy::Electron*>::iterator m_it = tight_electrons.begin();
	    m_it != tight_electrons.end(); m_it++){
	  if(isSameObject(it->momentum,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	// loose & tight ID variables
	bool loosePhoton = false;
	bool tightPhoton = false;

	std::map<TString,UChar_t>::iterator id_it = it->idPairs.find("PhotonCutBasedIDLoose");
	if(id_it != it->idPairs.end() && int(id_it->second) != 0) loosePhoton = true;

	id_it = it->idPairs.find("PhotonCutBasedIDTight");
	if(id_it != it->idPairs.end() && int(id_it->second) != 0) tightPhoton = true;

	if(loosePhoton) {
	  loose_photons.push_back(&*it);
	}

	if(tightPhoton) {
	  tight_photons.push_back(&*it);
	}

      }// for photon
    }// else


    if(printLevel > 0) std::cout << "Find loose and tight jets in the event." << std::endl;
      
    std::map<TString,susy::CaloJetCollection>::iterator jetColl_it = event->caloJets.find("ak5");
    if(jetColl_it == event->caloJets.end()){
      if(event->caloJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
    }
    else {

      susy::CaloJetCollection& jetColl = jetColl_it->second;

      for(std::vector<susy::CaloJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	float pt = it->momentum.Pt();
	if(pt < 30) continue;

	if(it->emEnergyFraction > 0.9) continue;

	if(std::abs(it->momentum.Eta()) > 2.4) continue;

	bool same = false;
	for(std::vector<susy::Muon*>::iterator m_it = tight_muons.begin();
	    m_it != tight_muons.end(); m_it++){
	  if(isSameObject(it->momentum,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	for(std::vector<susy::Electron*>::iterator m_it = tight_electrons.begin();
	    m_it != tight_electrons.end(); m_it++){
	  if(isSameObject(it->momentum,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(it->momentum,(*m_it)->momentum)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;


	loose_jets.push_back(&*it);

	tight_jets.push_back(&*it);

      }// for jet
    }// else

    if(printLevel > 0) std::cout << "Apply trigger selection in the event." << std::endl;

    bool passHLT = (useTrigger ? PassTriggers() : true);

    if(printLevel > 0) std::cout << "Select which met will be used in the event." << std::endl;

    std::map<TString, susy::MET>::iterator met_it = event->metMap.find("tcMet");
    if(met_it == event->metMap.end()) {
      std::cout << "MET map is not available!!!" << std::endl;
      continue;
    }
    susy::MET* met = &(met_it->second);

    if(printLevel > 0) {
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "              event summary" << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "loose_muons       : " << loose_muons.size() << std::endl;
      std::cout << "loose_electrons   : " << loose_electrons.size() << std::endl;
      std::cout << "loose_photons     : " << loose_photons.size() << std::endl;
      std::cout << "loose_jets        : " << loose_jets.size() << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "tight_muons       : " << tight_muons.size() << std::endl;
      std::cout << "tight_electrons   : " << tight_electrons.size() << std::endl;
      std::cout << "tight_photons     : " << tight_photons.size() << std::endl;
      std::cout << "tight_jets        : " << tight_jets.size() << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "met               : " << met->met() << std::endl;
    } 


    if(printLevel > 0) std::cout << "Apply event level cuts from now on..." << std::endl;


    // filter conditions

    if(enableFilter) {
      bool filterThis = (met->met() > 20.0) && (loose_photons.size() > 0);
      if(filterThis) {
	nFiltered++;
	filterTree->Fill();
      }
    }// if(enableFilter)


    nCnt[0]++; // total number of events

    if(!passHLT) continue;

    nCnt[1]++;

    if(tight_photons.size() == 0) continue;

    nCnt[2]++;

    h_met->Fill(met->met());
    h_sumEt->Fill(met->sumEt);


    if(met->met() < 20.0) continue;

    nCnt[3]++;

  } // for jentry


  // end of event loop and print summary

  std::cout << " ----------------- Job Summary ----------------- " << std::endl;
  std::cout << " Total events            : " << nCnt[0] << " (" << nCnt[0]/float(nCnt[0]) << ")" << std::endl;
  std::cout << " HLT passed              : " << nCnt[1] << " (" << nCnt[1]/float(nCnt[0]) << ")" << std::endl;
  std::cout << " photons > 0             : " << nCnt[2] << " (" << nCnt[2]/float(nCnt[0]) << ")" << std::endl;
  std::cout << " met > 20 GeV            : " << nCnt[3] << " (" << nCnt[3]/float(nCnt[0]) << ")" << std::endl;
  if(enableFilter) {
    std::cout << " --------------- Filtered events --------------- " << std::endl;
    std::cout << " filtered events         : " << nFiltered << " (" << nFiltered/float(nCnt[0]) << ")" << std::endl;
  }
  std::cout << " ----------------------------------------------- " << std::endl;

  // close the output file

  fout->cd();
  fout->Write();
  fout->Close();

  if(enableFilter) {
    filterTree->GetCurrentFile()->cd();
    filterTree->GetCurrentFile()->Write();
    filterTree->GetCurrentFile()->Close();
  }

}

