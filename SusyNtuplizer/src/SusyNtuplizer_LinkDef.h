// -*- C++ -*-
//
// Package:    SusyNtuplizer
/*

 Description: Dictionary generator for SusyEvent ojects

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyNtuplizer_LinkDef.h,v 1.9 2012/11/28 12:01:04 yiiyama Exp $
//

#include "../src/SusyEvent.h"

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedefs;

//// LINKS NEEDED TO BE ABLE TO WRITE SUSYTREE

#pragma link C++ class std::pair<TString, UChar_t>+;
#pragma link C++ class std::map<TString, UChar_t>+;
#pragma link C++ class std::pair<TString, Float_t>+;
#pragma link C++ class std::map<TString, Float_t>+;
#pragma link C++ class std::pair<TString, TVector3>+;
#pragma link C++ class std::map<TString, TVector3>+;
#pragma link C++ class std::pair<TString, TLorentzVector>+;
#pragma link C++ class std::map<TString, TLorentzVector>+;

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace susy;

#pragma link C++ class susy::PUSummaryInfo+;
#pragma link C++ class susy::PUSummaryInfoCollection+;
#pragma link C++ class susy::Particle+;
#pragma link C++ class susy::ParticleCollection+;
#pragma link C++ class susy::MET+;
#pragma link C++ class susy::Vertex+;
#pragma link C++ class susy::VertexCollection+;
#pragma link C++ class susy::Cluster+;
#pragma link C++ class susy::ClusterCollection+;
#pragma link C++ class susy::SuperCluster+;
#pragma link C++ class susy::SuperClusterCollection+;
#pragma link C++ class susy::Track+;
#pragma link C++ class susy::TrackCollection+;
#pragma link C++ class susy::Muon+;
#pragma link C++ class susy::MuonCollection+;
#pragma link C++ class susy::Electron+;
#pragma link C++ class susy::ElectronCollection+;
#pragma link C++ class susy::Photon+;
#pragma link C++ class susy::PhotonCollection+;
#pragma link C++ class susy::CaloJet+;
#pragma link C++ class susy::CaloJetCollection+;
#pragma link C++ class susy::PFJet+;
#pragma link C++ class susy::PFJetCollection+;
#pragma link C++ class susy::JPTJet+;
#pragma link C++ class susy::JPTJetCollection+;
#pragma link C++ class susy::PFParticle+;
#pragma link C++ class susy::PFParticleCollection+;

//// LINKS NEEDED TO BE ABLE TO READ SUSYTREE FROM INTERACTIVE SESSION

#pragma link C++ enum susy::MetFilters;
#pragma link C++ enum susy::BTagDiscriminators;
#pragma link C++ enum susy::PUJetIdAlgorithms;

#pragma link C++ class std::vector<const susy::Cluster*>+;
#pragma link C++ class std::vector<const susy::Track*>+;
#pragma link C++ class std::vector<const susy::PFParticle*>+;

#pragma link C++ class std::map<TString, susy::MET>+;
#pragma link C++ class std::map<TString, susy::MuonCollection>+;
#pragma link C++ class std::map<TString, susy::ElectronCollection>+;
#pragma link C++ class std::map<TString, susy::PhotonCollection>+;
#pragma link C++ class std::map<TString, susy::CaloJetCollection>+;
#pragma link C++ class std::map<TString, susy::PFJetCollection>+;
#pragma link C++ class std::map<TString, susy::JPTJetCollection>+;
#pragma link C++ class std::map<TString, susy::PFParticleCollection>+;
#pragma link C++ class std::pair<UInt_t, Bool_t>+;
#pragma link C++ class susy::TriggerMap+;
#pragma link C++ class susy::Event+;

#endif
