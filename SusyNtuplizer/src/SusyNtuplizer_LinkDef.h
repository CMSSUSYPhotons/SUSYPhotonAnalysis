// -*- C++ -*-
//
// Package:    SusyNtuplizer
/*

 Description: Dictionary generator for SusyEvent ojects

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyNtuplizer_LinkDef.h,v 1.5 2011/05/13 22:37:59 dwjang Exp $
//

#include "SusyEvent.h"

#ifdef __CINT__

#pragma link off all   globals;
#pragma link off all   classes;
#pragma link off all   functions;

#pragma link C++ class  std::pair<Int_t, UChar_t>+;
#pragma link C++ class  std::map<TString, UChar_t>+;
#pragma link C++ class  std::map<TString, Float_t>+;
#pragma link C++ class  std::map<TString, std::pair<Int_t, UChar_t> >+;
#pragma link C++ class  std::pair<TString, std::pair<Int_t, UChar_t> >+;
#pragma link C++ class  std::map<TString, susy::MET>+;
#pragma link C++ class  std::map<TString, TVector3>+;
#pragma link C++ class  std::map<TString, TLorentzVector>+;
#pragma link C++ class  std::pair<TString, UChar_t>+;
#pragma link C++ class  std::pair<TString, Float_t>+;
#pragma link C++ class  std::pair<TString, susy::MET>+;
#pragma link C++ class  std::pair<TString, TVector3>+;
#pragma link C++ class  std::pair<TString, TLorentzVector>+;
#pragma link C++ class  susy::Particle+;
#pragma link C++ class  susy::CorrMETData+;
#pragma link C++ class  susy::MET+;
#pragma link C++ class  susy::Cluster+;
#pragma link C++ class  susy::SuperCluster+;
#pragma link C++ class  susy::Track+;
#pragma link C++ class  susy::Muon+;
#pragma link C++ class  susy::Photon+;
#pragma link C++ class  susy::Electron+;
#pragma link C++ class  susy::CaloJet+;
#pragma link C++ class  susy::PFJet+;
#pragma link C++ class  susy::JPTJet+;
#pragma link C++ class  susy::PFParticle+;
#pragma link C++ class  std::map<TString, std::vector<susy::Electron> >+;
#pragma link C++ class  std::map<TString, std::vector<susy::Photon> >+;
#pragma link C++ class  std::map<TString, std::vector<susy::CaloJet> >+;
#pragma link C++ class  std::map<TString, std::vector<susy::PFJet> >+;
#pragma link C++ class  std::map<TString, std::vector<susy::JPTJet> >+;
#pragma link C++ class  std::map<TString, std::vector<susy::PFParticle> >+;
#pragma link C++ class  susy::Event+;

#endif
