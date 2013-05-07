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
// $Id: SusyEventAnalyzer.h,v 1.8 2013/05/05 12:01:32 yiiyama Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TChain.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TArrayI.h>

#include <map>
#include <set>
#include <fstream>

#include "../src/SusyEvent.h"

class SusyEventAnalyzer {
public:
  SusyEventAnalyzer(TTree&);
  virtual ~SusyEventAnalyzer();

  /* Main analyzer function to be defined */
  virtual void Run();

  /* parameter configuration functions */
  void IncludeAJson(TString const&);
  void SetOutput(TString const& v) { outputName = v; }
  void SetLogFile(TString const& v) { logFileName = v;}
  void SetPrintInterval(int v) { printInterval = v; }
  void SetPrintLevel(int v) { printLevel = v; }
  void SetProcessNEvents(int v) { processNEvents = v; }
  void AddHltName(TString const& v) { hltNames.push_back(v + "_v*"); }
  void CopyEvents(bool v) { copyEvents = v; }

protected:
  bool IsGoodLumi(UInt_t, UInt_t) const;
  bool PassTriggers() const;

  /* container of all event data */
  susy::Event event;
  /* input tree */
  TTree *fTree;
  /* suffix of the output file */
  TString outputName;
  /* log file name */
  TString logFileName;
  /* verbosity - 0 => no printout, 1 => print function control flow, 2 => print event processing flow, 3 => print event dump */
  int printLevel;
  /* print frequency */
  unsigned printInterval;
  /* maximum number of events */
  int processNEvents;
  /* HLT path names */
  std::vector<TString> hltNames;
  /* switch for saving skims */
  bool copyEvents;
  /* good lumi list */
  std::map<unsigned, std::set<unsigned> > goodLumiList;
  mutable std::pair<unsigned, unsigned> currentLumi;
  mutable bool currentLumiIsGood;
};

SusyEventAnalyzer::SusyEventAnalyzer(TTree& tree) :
  event(),
  fTree(&tree),
  outputName("analysis"),
  logFileName(outputName + ".log"),
  printLevel(0),
  printInterval(1000),
  processNEvents(-1),
  hltNames(),
  copyEvents(false),
  goodLumiList(),
  currentLumi(0, 0),
  currentLumiIsGood(true)
{
  event.setInput(tree);
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
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
    if(event.hltMap.pass(hltNames[iT])) return true;

  return false;
}

#endif

