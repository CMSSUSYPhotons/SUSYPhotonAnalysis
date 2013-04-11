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
// $Id: SusyEventAnalyzer.h,v 1.6 2012/08/31 11:33:53 bfrancis Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TChain.h>
#include <TString.h>

#include <map>
#include <set>

#include "../src/SusyEvent.h"

class SusyEventAnalyzer {
 public :
  TTree          *fTree;   // pointer to the analyzed TTree or TChain

  // Declaration of leaf types
  susy::Event*    event;

  SusyEventAnalyzer(TTree&);
  virtual ~SusyEventAnalyzer();

  virtual bool InitializePerEvent();
  virtual void Run();                          // event loop for main analysis

  bool IsGoodLumi(UInt_t, UInt_t lumi) const;      // JSON based good run list cut...
  bool PassTriggers() const; // return true if any of names in hltNames are fired

  // parameter configuration functions
  void IncludeAJson(TString const&);  // Call to pull in a json file 
  void SetOutput(TString const& v) { outputName = v; }
  void SetPrintInterval(int v) { printInterval = v; }
  void SetPrintLevel(int v) { printLevel = v; }
  void SetProcessNEvents(int v) { processNEvents = v; }
  void AddHltName(TString const& v) { hltNames.push_back(v + "_v*"); }
  void CopyEvents(bool v) { copyEvents = v; }

 private:

  TString outputName;               // dataset name to be used for output histfile name

  // printLevel
  // 0 : default - no printout
  // 1 : print functional step in every event
  // 2 : print values in collections
  int printLevel;

  unsigned printInterval;           // print frequency
  int processNEvents;       // number of events to be processed
  std::vector<TString> hltNames;          // HLT trigger path names
  bool copyEvents;        // filter events of interest

  std::map<unsigned, std::set<unsigned> > goodLumiList;
  mutable std::pair<unsigned, unsigned> currentLumi;
  mutable bool currentLumiIsGood;
};

#endif

