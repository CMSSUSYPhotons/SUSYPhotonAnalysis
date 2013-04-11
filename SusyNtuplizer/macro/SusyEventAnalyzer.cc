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
// $Id: SusyEventAnalyzer.cc,v 1.15 2012/08/31 11:33:53 bfrancis Exp $
//

#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TPRegexp.h>
#include <TFile.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>
#include <iostream>
#include <fstream>

#include "SusyEventAnalyzer.h"
#include "SusyEventPrinter.h"

template<typename T>
bool
EtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Et() > p2->momentum.Et());
}

template<typename T>
bool
PtGreater(const T* p1, const T* p2) {
  return (p1->momentum.Pt() > p2->momentum.Pt());
}

template<typename T1, typename T2>
bool
isSameObject(const T1& p1, const T2& p2)
{
  float dEta = p1.momentum.Eta() - p2.momentum.Eta();
  float dPhi = TVector2::Phi_mpi_pi(p1.momentum.Phi() - p2.momentum.Phi());
  float dR2 = dEta*dEta + dPhi*dPhi;
  if(dR2 < 0.25) return true;
  return false;
}

float
d0correction(TVector3& beamSpot, susy::Track& track)
{
  float d0 = track.d0() - beamSpot.X()*std::sin(track.phi()) + beamSpot.Y()*std::cos(track.phi());
  return d0;
}


SusyEventAnalyzer::SusyEventAnalyzer(TTree& tree) :
  fTree(&tree),
  event(0),
  outputName("analysis"),
  printLevel(0),
  printInterval(1000),
  processNEvents(-1),
  hltNames(),
  copyEvents(false),
  goodLumiList(),
  currentLumi(0, 0),
  currentLumiIsGood(true)
{
  event = new susy::Event;
  event->bindTree(tree);
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
  delete event;
}

void
SusyEventAnalyzer::IncludeAJson(TString const& _fileName)
{
  if(_fileName == "") return;

  std::ifstream inputFile(_fileName);
  if(!inputFile.is_open()){
    std::cerr << "Cannot open JSON file " << _fileName << std::endl;
    return;
  }

  std::string line;
  TString jsonText;
  while(true){
    std::getline(inputFile, line);
    if(!inputFile.good()) break;
    jsonText += line;
  }
  inputFile.close();

  TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
  TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

  TArrayI positions(2);
  positions[1] = 0;
  while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
    TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
    TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

    unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
    std::set<unsigned>& lumis(goodLumiList[run]);

    TArrayI lumiPos(2);
    lumiPos[1] = 0;
    while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
      TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
      int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
      int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
      for(int lumi(begin); lumi <= end; ++lumi)
        lumis.insert(lumi);
    }
  }
}

bool
SusyEventAnalyzer::IsGoodLumi(UInt_t run, UInt_t lumi) const
{
  if(goodLumiList.size() == 0) return true;
  if(run == currentLumi.first && lumi == currentLumi.second) return currentLumiIsGood;
  currentLumi.first = run;
  currentLumi.second = lumi;
  currentLumiIsGood = false;

  std::map<unsigned, std::set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
  if(rItr != goodLumiList.end()){
    std::set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
    if(lItr != rItr->second.end()) currentLumiIsGood = true;
  }

  return currentLumiIsGood;
}

bool
SusyEventAnalyzer::PassTriggers() const
{
  unsigned nT(hltNames.size());
  if(nT == 0) return true;

  for(unsigned iT(0); iT != nT; ++iT)
    if(event->hltMap.pass(hltNames[iT])) return true;

  return false;
}

bool
SusyEventAnalyzer::InitializePerEvent()
{
  return true;
}

void
SusyEventAnalyzer::Run()
{
  const int NCNT = 20;
  int nCnt[NCNT];
  for(int i=0; i<NCNT; i++) nCnt[i] = 0;

  if(printLevel > 0) std::cout << "Open file for histograms" << std::endl;
  TFile* fout = TFile::Open("hist_" + outputName + ".root", "RECREATE");
  if(!fout || fout->IsZombie()){
    std::cerr << "Cannot open output file hist_" << outputName << ".root" << std::endl;
    return;
  }

  TTree* copyTree = 0;

  if(copyEvents){
    if(printLevel > 0) std::cout << "Open file for skim output" << std::endl;
    TFile* filterFile = TFile::Open("susyEvents_" + outputName + ".root", "RECREATE");
    if(!filterFile || filterFile->IsZombie())
      std::cerr << "Cannot open output file susyEvents_" << outputName << ".root" << std::endl;
    else{
      event->bindTree(*copyTree, false);
      copyTree->SetAutoSave(10000000);
    }
  }

  if(printLevel > 0) std::cout << "Define histograms" << std::endl;
    
  fout->cd();

  TH1F* h_vtxZ = new TH1F("vtxZ","Z position of the primary vertex;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_bsZ = new TH1F("bsZ","Z position of the beam spot;Z (cm);Events",100,-50.0,50.0);
  TH1F* h_met = new TH1F("met","missing transverse energy;#slash{E}_{T} (GeV);Events",200,0.0,1000.0);
  TH1F* h_sumEt = new TH1F("sumEt","Scalar sum of all calorimeter energy;#sigmaE_{T} (GeV);Events",200,0.0,2000.0);

  if(printLevel > 0) std::cout << "Start event loop" << std::endl;

  long iEntry(0);
  while(iEntry != processNEvents && fTree->GetEntry(iEntry++) != 0){

    if(printLevel > 0 || iEntry % printInterval == 0)
      std::cout << iEntry - 1 << " events processed with run=" << event->runNumber << ", event=" << event->eventNumber << std::endl;

    if(printLevel > 0) std::cout << "Initialize any global variables to be reset per event->" << std::endl;

    // remove events that fail initialization
    if(!InitializePerEvent()) continue;

    if(printLevel > 1) Print(*event);

    if(printLevel > 0) std::cout << "Apply good run list." << std::endl;

    // remove events not in good lumi list
    if(!IsGoodLumi(event->runNumber, event->luminosityBlockNumber)) continue;

    if(printLevel > 0) std::cout << "Apply MET filter." << std::endl;

    // remove events filtered by optional met filters
    if(event->isRealData && !event->passMetFilters()) continue;

    nCnt[0]++; // total number of events

    if(printLevel > 0) std::cout << "Apply HLT cut." << std::endl;

    if(!PassTriggers()) continue;

    if(printLevel > 0) std::cout << "Event passes presele." << std::endl;

    nCnt[1]++; // total number of events

    if(printLevel > 0) std::cout << "Set object references in Event" << std::endl;

    event->fillRefs();

    // classify photon objects

    // loose objects have all standard cuts except for isolation
    std::vector<susy::Photon*>   loose_photons;

    // tight objects hava isolation cuts applied on top of loose objects
    std::vector<susy::Photon*>   tight_photons;

    // same as tight except for nPixelSeeds > 0
    std::vector<susy::Photon*>   ele_photons;

    // same as tight except for reversing either trackIso or sigmaIetaIeta
    std::vector<susy::Photon*>   fake_photons;

    std::vector<susy::CaloJet*>  caloJets;
    std::vector<susy::PFJet*>    pfJets;

    if(printLevel > 0) std::cout << "Find primary vertex in the event->" << std::endl;

    TVector3* primVtx = 0;
    if(event->vertices.size() > 0) primVtx = &(event->vertices[0].position);

    if(primVtx) h_vtxZ->Fill(primVtx->Z());
    h_bsZ->Fill(event->beamSpot.Z());

    if(printLevel > 0) std::cout << "Find loose and tight photons in the event->" << std::endl;

    std::map<TString, susy::PhotonCollection>::iterator phoMap = event->photons.find("photons");

    if(phoMap != event->photons.end()) {

      susy::PhotonCollection& phoColl = phoMap->second;

      for(susy::PhotonCollection::iterator it = phoColl.begin();
	  it != phoColl.end(); it++) {

	// fiducial cuts. Look for only barrel now
	if(!it->isEB()) continue;

	// Et cuts, 25 GeV for trailing photons. Will apply tighter for the leading one.
	if(it->momentum.Et() < 25.0) continue;

        // optional Spike cleaning
        if(it->r9 > 1.0) continue;

        // H/E (in trigger, 0.15 for EB, 0.10 for EE)
        bool heCut = (it->hadronicOverEm < 0.05);
        
        // sigma_ietaieta (in trigger 0.014 for EB, 0.034 for EE)
        bool sIetaCut = (it->sigmaIetaIeta < 0.013);

        // Ecal Isolation
        bool ecalIsoCut = (it->ecalRecHitSumEtConeDR04 < 4.2 + 0.006 * it->momentum.Et());

        // Hcal Isolation
        bool hcalIsoCut = (it->hcalTowerSumEtConeDR04() < 2.2 + 0.0025 * it->momentum.Et());

        // Track Isolation
        bool trackIsoCut = (it->trkSumPtHollowConeDR04 < 2.0 + 0.001 * it->momentum.Et());

        bool pixelCut = (it->nPixelSeeds == 0);

        // loose & tight ID variables
        bool looseCut = heCut && ecalIsoCut && hcalIsoCut;
        bool tightCut = looseCut && pixelCut && sIetaCut && trackIsoCut;
        bool eleClass  = looseCut && !pixelCut && sIetaCut && trackIsoCut;
        bool fakeClass = looseCut && pixelCut && !(sIetaCut && trackIsoCut);

        if(looseCut) {
          loose_photons.push_back(&*it);
        }
        if(tightCut) {
          tight_photons.push_back(&*it);
        }
        if(eleClass) {
          ele_photons.push_back(&*it);
        }
        if(fakeClass) {
          fake_photons.push_back(&*it);
        }

      }// for photon
    }// else

    // sort photons by Et
    std::sort(loose_photons.begin(),loose_photons.end(),EtGreater<susy::Photon>);
    std::sort(tight_photons.begin(),tight_photons.end(),EtGreater<susy::Photon>);
    std::sort(ele_photons.begin(),ele_photons.end(),EtGreater<susy::Photon>);
    std::sort(fake_photons.begin(),fake_photons.end(),EtGreater<susy::Photon>);


    if(printLevel > 0) std::cout << "Find caloJets in the event->" << std::endl;
      
    std::map<TString,susy::CaloJetCollection>::iterator caloJets_it = event->caloJets.find("ak5");

    if(caloJets_it != event->caloJets.end()){

      susy::CaloJetCollection& jetColl = caloJets_it->second;

      for(std::vector<susy::CaloJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
	if (s_it == it->jecScaleFactors.end()) {
	  std::cout << "JEC is not available for this jet!!!" << std::endl;
	  continue;
	}
	float scale = s_it->second;

        if(printLevel > 2) std::cout << "CaloJet stored (" << scale << ")" << std::endl;

	TLorentzVector corrP4 = scale * it->momentum;

	if(std::abs(corrP4.Eta()) > 3.0) continue;

	bool same = false;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(*it, **m_it)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	//	if(pt < 20) continue;

	caloJets.push_back(&*it);

      }// for jet
    }// else

    std::sort(caloJets.begin(),caloJets.end(),EtGreater<susy::CaloJet>);


    if(printLevel > 0) std::cout << "Find pfJets in the event->" << std::endl;
      
    std::map<TString,susy::PFJetCollection>::iterator pfJets_it = event->pfJets.find("ak5");
    if(pfJets_it == event->pfJets.end()){
      if(event->pfJets.size() > 0) std::cout << "JetCollection is not available!!!" << std::endl;
    }
    else {

      susy::PFJetCollection& jetColl = pfJets_it->second;

      for(std::vector<susy::PFJet>::iterator it = jetColl.begin();
	  it != jetColl.end(); it++) {

	std::map<TString,Float_t>::iterator s_it = it->jecScaleFactors.find("L2L3");
	if (s_it == it->jecScaleFactors.end()) {
	  std::cout << "JEC is not available for this jet!!!" << std::endl;
	  continue;
	}
	float scale = s_it->second;

        if(printLevel > 2) std::cout << "PFJet stored (" << scale << ")" << std::endl;

	TLorentzVector corrP4 = scale * it->momentum;

	if(std::abs(corrP4.Eta()) > 3.0) continue;

	bool same = false;

	for(std::vector<susy::Photon*>::iterator m_it = tight_photons.begin();
	    m_it != tight_photons.end(); m_it++){
	  if(isSameObject(*it, **m_it)){
	    same = true;
	    break;
	  }
	}
	if(same) continue;

	//	if(pt < 20) continue;

	pfJets.push_back(&*it);

      }// for jet
    }// else

    std::sort(pfJets.begin(),pfJets.end(),EtGreater<susy::PFJet>);


    if(printLevel > 0) std::cout << "Select which met will be used in the event->" << std::endl;

    std::map<TString, susy::MET>::iterator met_it = event->metMap.find("pfType01SysShiftCorrectedMet");
    if(met_it == event->metMap.end()) {
      std::cout << "MET is not available!!!" << std::endl;
      continue;
    }
    susy::MET& met = met_it->second;

    if(printLevel > 0) {
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "              event summary" << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "loose_photons     : " << loose_photons.size() << std::endl;
      std::cout << "tight_photons     : " << tight_photons.size() << std::endl;
      std::cout << "ele_photons       : " << ele_photons.size() << std::endl;
      std::cout << "fake_photons      : " << fake_photons.size() << std::endl;
      std::cout << "caloJets          : " << caloJets.size() << std::endl;
      std::cout << "pfJets            : " << pfJets.size() << std::endl;
      std::cout << "------------------------------------------" << std::endl;
      std::cout << "met               : " << met.met() << std::endl;
    } 



    if(printLevel > 0) std::cout << "Apply event level cuts from now on..." << std::endl;

    // Event level cuts

    if(loose_photons.size() == 0) continue;

    // End event level cuts

    if(copyEvents) copyTree->Fill();

    nCnt[2]++;

    h_met->Fill(met.met());
    h_sumEt->Fill(met.sumEt);

    // two photons
    if(tight_photons.size() >= 2) {
      nCnt[3]++;
    }

    // one photon + one electron
    if(tight_photons.size() >= 1 && ele_photons.size() >= 1) {
      nCnt[4]++;
    }

    // two electrons
    if(ele_photons.size() >= 2) {
      nCnt[5]++;
    }

    // one photon + one fake
    if(tight_photons.size() >= 1 && fake_photons.size() >= 1) {
      nCnt[6]++;
    }

    // two fakes
    if(fake_photons.size() >= 2) {
      nCnt[7]++;
    }

    if(met.met() < 50.0) continue;

    nCnt[8]++;

  } // for iEntry


  // end of event loop and print summary

  std::cout << " ----------------- Job Summary ----------------- " << std::endl;
  std::cout << " Total events            : " << nCnt[0] << std::endl;
  std::cout << " HLT passed              : " << nCnt[1] << " (" << nCnt[1]/float(nCnt[0]) << ") wrt total events" << std::endl;
  std::cout << " loose_photons > 0       : " << nCnt[2] << " (" << nCnt[2]/float(nCnt[1]) << ") wrt HLT" << std::endl;
  std::cout << " gg events               : " << nCnt[3] << " (" << nCnt[3]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ge events               : " << nCnt[4] << " (" << nCnt[4]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ee events               : " << nCnt[5] << " (" << nCnt[5]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " gf events               : " << nCnt[6] << " (" << nCnt[6]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " ff events               : " << nCnt[7] << " (" << nCnt[7]/float(nCnt[1]) << ")" << std::endl;
  std::cout << " met > 50 GeV            : " << nCnt[8] << " (" << nCnt[8]/float(nCnt[1]) << ")" << std::endl;

  if(copyEvents){
    std::cout << " --------------- Filtered events --------------- " << std::endl;
    std::cout << " filtered events         : " << copyTree->GetEntries() << " (" << copyTree->GetEntries()/float(nCnt[0]) << ")" << std::endl;
  }
  std::cout << " ----------------------------------------------- " << std::endl;

  // close the output file

  fout->cd();
  fout->Write();
  delete fout;

  if(copyEvents){
    TFile* copyFile(copyTree->GetCurrentFile());
    copyFile->cd();
    copyFile->Write();
    event->releaseTree(*copyTree);
    delete copyFile;
  }

}

