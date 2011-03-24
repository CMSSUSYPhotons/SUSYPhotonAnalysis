// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.h
// 
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.h,v 1.1 2011/03/24 23:46:27 dwjang Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <iostream>

#include "../src/SusyEvent.h"


class SusyEventAnalyzer {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  susy::Event     *event;

  // List of branches
  TBranch        *b_Event;

  SusyEventAnalyzer(TTree *tree=0);
  virtual ~SusyEventAnalyzer();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();                          // event loop for main analysis

  // utility functions
  bool isGoodRun(int run);
  bool isSameObject(TLorentzVector& p1, TLorentzVector& p2);
  float d0correction(TVector3& beamSpot, susy::Track& track) const;

  // parameter configuration functions
  void Initialize();         // global variables needed to be initialized just once
  void InitializePerEvent(); // global variables needed to be initialized per event
  void SetDataset(TString& v) {          ds = v; }
  void SetPrintInterval(int v) {         printInterval = v; }
  void SetPrintLevel(int v) {            printLevel = v; }
  void SetProcessNEvents(int v) {        processNEvents = v; }
  void SetUseTrigger(bool v) {           useTrigger = v; }
  void SetHltName(TString v) {           hltName = v; }
  void SetFilter(bool v) {               enableFilter = v; }
  void SetFilteredFileName(TString v) {  filtered_file_name = v; }

 private:

  TString ds;               // dataset name to be used for output histfile name

  // printLevel
  // 0 : default - no printout
  // 1 : print functional step in every event
  // 2 : print values in collections
  int printLevel;           // print frequency

  int printInterval;        // print level for event content: defined in Event.h
  int processNEvents;       // number of events to be processed
  bool useTrigger;          // flag for using trigger bit selection.
  TString hltName;          // HLT trigger path name
  bool enableFilter;        // filter events of interest
  TString filtered_file_name; // filtered output file name

};

#endif

#ifdef SusyEventAnalyzer_cxx
SusyEventAnalyzer::SusyEventAnalyzer(TTree *tree)
{
  if (tree == 0) {
    std::cout << "Error!!! There is no file containing a tree." << std::endl;
  }
  Init(tree);
  Initialize();
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t SusyEventAnalyzer::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t SusyEventAnalyzer::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
  }
  return centry;
}

void SusyEventAnalyzer::Init(TTree *tree)
{
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  //   fChain->SetMakeClass(1);

  event = new susy::Event;

  fChain->SetBranchAddress("susyEvent", &event, &b_Event);
}

void SusyEventAnalyzer::Initialize() {

  ds = "test";
  printLevel = 0;
  printInterval = 1000;
  processNEvents = -1;
  useTrigger = false;
  enableFilter = false;
  filtered_file_name = "filtered.root";

}

#endif // #ifdef SusyEventAnalyzer_cxx
