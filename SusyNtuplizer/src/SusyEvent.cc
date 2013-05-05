// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEvent.cc
//
/*

Description: Objects definitions used for SusyNtuples

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEvent.cc,v 1.31 2013/04/12 09:53:27 yiiyama Exp $
//

#include "SusyEvent.h"

#include "TTree.h"
#include "TChain.h"

#include <limits>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <cstdio>

// Use for initialization of variables that can have a value 0 and are not always filled
float const BIGVALUE(std::numeric_limits<float>::max());

// Print utility functions
std::ostream& indent(std::ostream& os)
{ return os << "\t"; }
template<class T> std::ostream& operator<<(std::ostream& os, std::vector<T> const& vect)
{ for(unsigned i(0); i != vect.size(); ++i) os << vect[i] << ", "; return os; }
template<class T> std::ostream& operator<<(std::ostream& os, std::map<TString, T> const& map)
{ for(typename std::map<TString, T>::const_iterator itr(map.begin()); itr != map.end(); ++itr) os << "(\"" << itr->first << "\" => " << itr->second << "), "; return os; }
std::ostream& operator<<(std::ostream& os, TVector2 const& v)
{ os << "(x, y) = (" << v.X() << ", " << v.Y() << ") (rho, phi) = (" << v.Mod() << ", " << v.Phi() << ")"; return os; }
std::ostream& operator<<(std::ostream& os, TVector3 const& v)
{ os << "(x, y, z) = (" << v.X() << ", " << v.Y() << ", " << v.Z() << ") (rho, theta, phi) = (" << v.Mag() << ", " << v.Theta() << ", " << v.Phi() << ")"; return os; }
std::ostream& operator<<(std::ostream& os, TLorentzVector const& v)
{ os << "(x, y, z, t) = (" << v.X() << ", " << v.Y() << ", " << v.Z() << ", " << v.T() << ") (P, eta, phi, E) = (" << v.P() << ", " << v.Eta() << ", " << v.Phi() << ", " << v.E() << ")"; return os; }
template<class T> TString bin(T bits)
{ unsigned const len(sizeof(T) * 8); TString str; for(unsigned i(0); i != len; ++i) str.Append('0' + ((bits >> (len - i - 1)) & 1)); return str; }

void
susy::PUSummaryInfo::Init()
{
  BX = 0;
  numInteractions = 0;
  trueNumInteractions = 0;
  zPositions.clear();
  sumPTLowPT.clear();
  sumPTHighPT.clear();
  numTracksLowPT.clear();
  numTracksHighPT.clear();
  instLumi.clear();
  dataMixerRun.clear();
  dataMixerEvt.clear();
  dataMixerLumiSection.clear();
}

void
susy::PUSummaryInfo::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "BX: " << int(BX) << std::endl;
  indent(os) << "numInteractions: " << int(numInteractions) << std::endl;
  indent(os) << "trueNumInteractions: " << trueNumInteractions << std::endl;
  indent(os) << "zPositions: " << zPositions << std::endl;
  indent(os) << "sumPTLowPT: " << sumPTLowPT << std::endl;
  indent(os) << "sumPTHighPT: " << sumPTHighPT << std::endl;
  indent(os) << "numTracksLowPT: " << numTracksLowPT << std::endl;
  indent(os) << "numTracksHighPT: " << numTracksHighPT << std::endl;
  indent(os) << "instLumi: " << instLumi << std::endl;
  indent(os) << "dataMixerRun: " << dataMixerRun << std::endl;
  indent(os) << "dataMixerEvt: " << dataMixerEvt << std::endl;
  indent(os) << "dataMixerLumiSection: " << dataMixerLumiSection << std::endl;
}

void
susy::Particle::Init()
{
  status         = 0;
  charge         = 0;
  motherIndex    = -1;
  pdgId          = 0;
  vertex        *= 0;
  momentum      *= 0;

  mother         = 0;
}

void
susy::Particle::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "status: " << int(status) << std::endl;
  indent(os) << "charge: " << int(charge) << std::endl;
  indent(os) << "motherIndex: " << motherIndex << std::endl;
  indent(os) << "pdgId: " << pdgId << std::endl;
  indent(os) << "vertex: " << vertex << std::endl;
  indent(os) << "momentum: " << momentum << std::endl;
}

void
susy::Particle::fillRefs(Event const* _evt)
{
  mother = motherIndex != -1 ? &_evt->genParticles[motherIndex] : 0;
}

void
susy::PFParticle::Init()
{
  pdgId                       = 0;
  charge                      = 0;
  isPU                        = kFALSE;
  ecalEnergy                  = 0;
  hcalEnergy                  = 0;

  vertex                     *= 0;
  momentum                   *= 0;
}

void
susy::PFParticle::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "pdgId: " << pdgId << std::endl;
  indent(os) << "charge: " << int(charge) << std::endl;
  indent(os) << "isPU: " << isPU << std::endl;
  indent(os) << "ecalEnergy: " << ecalEnergy << std::endl;
  indent(os) << "hcalEnergy: " << hcalEnergy << std::endl;
  indent(os) << "vertex: " << vertex << std::endl;
  indent(os) << "momentum: " << momentum << std::endl;
}

void
susy::MET::Init()
{
  sumEt        = 0;
  significance = 0;
  mEt         *= 0;
}

void
susy::MET::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "sumEt: " << sumEt << std::endl;
  indent(os) << "significance: " << significance << std::endl;
  indent(os) << "mEt: " << mEt << std::endl;
}

void
susy::Vertex::Init()
{
  tracksSize = 0;
  sumPt2     = 0;
  chi2       = 0;
  ndof       = 0;
  position  *= 0;
}

void
susy::Vertex::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "tracksSize: " << tracksSize << std::endl;
  indent(os) << "sumPt2: " << sumPt2 << std::endl;
  indent(os) << "chi2: " << chi2 << std::endl;
  indent(os) << "ndof: " << ndof << std::endl;
  indent(os) << "position: " << position << std::endl;
}

void
susy::Cluster::Init()
{
  nCrystals = 0;
  energy    = 0;
  position *= 0;
}

void
susy::Cluster::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "nCrystals: " << int(nCrystals) << std::endl;
  indent(os) << "energy: " << energy << std::endl;
  indent(os) << "position: " << position << std::endl;
}

void
susy::SuperCluster::Init()
{
  seedClusterIndex = -1;
  energy           = 0;
  preshowerEnergy  = 0;
  phiWidth         = 0;
  etaWidth         = 0;
  position        *= 0;
  basicClusterIndices.clear();

  seedCluster = 0;
  basicClusters.clear();
}

void
susy::SuperCluster::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "seedClusterIndex: " << seedClusterIndex << std::endl;
  indent(os) << "energy: " << energy << std::endl;
  indent(os) << "preshowerEnergy: " << preshowerEnergy << std::endl;
  indent(os) << "phiWidth: " << phiWidth << std::endl;
  indent(os) << "etaWidth: " << etaWidth << std::endl;
  indent(os) << "position: " << position << std::endl;
  indent(os) << "basicClusterIndices: " << basicClusterIndices << std::endl;
}

void
susy::SuperCluster::fillRefs(Event const* _evt)
{
  seedCluster = seedClusterIndex != -1 ? &_evt->clusters[seedClusterIndex] : 0;

  basicClusters.assign(basicClusterIndices.size(), 0);
  for(unsigned iC(0); iC != basicClusterIndices.size(); ++iC)
    basicClusters[iC] = &_evt->clusters[basicClusterIndices[iC]];
}

void
susy::Track::Init()
{
  algorithm                        = 0;
  quality                          = 0;
  numberOfValidHits                = 0;
  numberOfValidTrackerHits         = 0;
  numberOfValidMuonHits            = 0;
  numberOfValidPixelHits           = 0;
  numberOfValidStripHits           = 0;
  vertexIndex                      = -1;
  chi2                             = 0;
  ndof                             = 0;
  charge                           = 0;
  for (int i=0; i<5; i++) error[i] = 0;
  ptError                          = 0;
  vertex                          *= 0;
  momentum                        *= 0;

  assignedVertex = 0;
}

void
susy::Track::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "algorithm: " << int(algorithm) << std::endl;
  indent(os) << "quality: " << int(quality) << std::endl;
  indent(os) << "numberOfValidHits: " << int(numberOfValidHits) << std::endl;
  indent(os) << "numberOfValidTrackerHits: " << int(numberOfValidTrackerHits) << std::endl;
  indent(os) << "numberOfValidMuonHits: " << int(numberOfValidMuonHits) << std::endl;
  indent(os) << "numberOfValidPixelHits: " << int(numberOfValidPixelHits) << std::endl;
  indent(os) << "numberOfValidStripHits: " << int(numberOfValidStripHits) << std::endl;
  indent(os) << "vertexIndex: " << vertexIndex << std::endl;
  indent(os) << "chi2: " << chi2 << std::endl;
  indent(os) << "ndof: " << ndof << std::endl;
  indent(os) << "charge: " << charge << std::endl;
  indent(os) << "error: ";
  for(unsigned i(0); i != 5; ++i) os << i << ". " << error[i] << " ";
  os << std::endl;
  indent(os) << "ptError: " << ptError << std::endl;
  indent(os) << "vertex: " << vertex << std::endl;
  indent(os) << "momentum: " << momentum << std::endl;
}

void
susy::Track::fillRefs(Event const* _evt)
{
  assignedVertex = vertexIndex != -1 ? &_evt->vertices[vertexIndex] : 0;
}

void
susy::Photon::Init()
{
  fidBit                          = 0;
  nPixelSeeds                     = 0;
  passelectronveto                = 0;
  hadronicOverEm                  = 0;
  hadTowOverEm                    = 0;
  hadronicDepth1OverEm            = 0;
  hadronicDepth2OverEm            = 0;
  e1x2                            = BIGVALUE;
  e1x5                            = 0;
  e2x5                            = 0;
  e3x3                            = 0;
  e5x5                            = 0;
  maxEnergyXtal                   = 0;
  sigmaEtaEta                     = 0;
  sigmaIetaIeta                   = 0;
  sigmaIphiIphi                   = BIGVALUE;
  r9                              = 0;

  ecalRecHitSumEtConeDR04         = 0;
  hcalDepth1TowerSumEtConeDR04    = 0;
  hcalDepth2TowerSumEtConeDR04    = 0;
  hcalIsoConeDR04_2012            = 0;
  trkSumPtSolidConeDR04           = 0;
  trkSumPtHollowConeDR04          = 0;
  nTrkSolidConeDR04               = 0;
  nTrkHollowConeDR04              = 0;

  ecalRecHitSumEtConeDR03         = 0;
  hcalDepth1TowerSumEtConeDR03    = 0;
  hcalDepth2TowerSumEtConeDR03    = 0;
  hcalIsoConeDR03_2012            = 0;
  trkSumPtSolidConeDR03           = 0;
  trkSumPtHollowConeDR03          = 0;
  nTrkSolidConeDR03               = 0;
  nTrkHollowConeDR03              = 0;

  chargedHadronIso                = 0;
  neutralHadronIso                = 0;
  photonIso                       = 0;

  worstOtherVtxChargedHadronIso   = 0.;
  worstOtherVtxChargedHadronIsoVtxIdx = -1;

  chargedHadronIsoDeposit         = BIGVALUE;
  neutralHadronIsoDeposit         = BIGVALUE;
  photonIsoDeposit                = BIGVALUE;

  seedTime                        = BIGVALUE;

  mipChi2                         = 0.;
  mipTotEnergy                    = 0.;
  mipSlope                        = 0.;
  mipIntercept                    = 0.;
  mipNhitCone                     = 0.;
  mipIsHalo                       = kFALSE;

  convInfo                        = kFALSE;
  convDist                        = 0.;
  convDcot                        = 0.;
  convVtxChi2                     = -1.;
  convVtxNdof                     = 0;
  convVertex                     *= 0;
  convDxy                         = 0.;
  convDz                          = 0.;
  convLxy                         = 0.;
  convLz                          = 0.;
  convZofPVFromTracks             = 0.;
  convTrackChargeProd             = 0;
  convTrack1nHit                  = 0;
  convTrack2nHit                  = 0;
  convTrack1chi2                  = 0.;
  convTrack2chi2                  = 0.;
  convTrack1pT                    = 0.;
  convTrack2pT                    = 0.;
  convTrack1InnerZ                = 0.;
  convTrack2InnerZ                = 0.;
  convTrack1InnerX                = 0.;
  convTrack2InnerX                = 0.;
  convTrack1InnerY                = 0.;
  convTrack2InnerY                = 0.;
  convTrack1Signedd0              = 0.;
  convTrack2Signedd0              = 0.;

  superClusterIndex               = -1;
  superClusterPreshowerEnergy     = 0;
  superClusterPhiWidth            = 0;
  superClusterEtaWidth            = 0;
  caloPosition                   *= 0;

  MVAregEnergyAndErr.first        = 0;
  MVAregEnergyAndErr.second       = 0;
  MVAcorrMomentum                *= 0;

  momentum                       *= 0;

  superCluster = 0;
  worstOtherVtxChargedHadronIsoVtx = 0;
}

void
susy::Photon::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "fidBit: " << bin(fidBit) << std::endl;
  indent(os) << "nPixelSeeds: " << nPixelSeeds << std::endl;
  indent(os) << "passelectronveto: " << passelectronveto << std::endl;
  indent(os) << "hadronicOverEm: " << hadronicOverEm << std::endl;
  indent(os) << "hadTowOverEm: " << hadTowOverEm << std::endl;
  indent(os) << "hadronicDepth1OverEm: " << hadronicDepth1OverEm << std::endl;
  indent(os) << "hadronicDepth2OverEm: " << hadronicDepth2OverEm << std::endl;
  indent(os) << "e1x2: " << e1x2 << std::endl;
  indent(os) << "e1x5: " << e1x5 << std::endl;
  indent(os) << "e2x5: " << e2x5 << std::endl;
  indent(os) << "e3x3: " << e3x3 << std::endl;
  indent(os) << "e5x5: " << e5x5 << std::endl;
  indent(os) << "maxEnergyXtal: " << maxEnergyXtal << std::endl;
  indent(os) << "sigmaEtaEta: " << sigmaEtaEta << std::endl;
  indent(os) << "sigmaIetaIeta: " << sigmaIetaIeta << std::endl;
  indent(os) << "sigmaIphiIphi: " << sigmaIphiIphi << std::endl;
  indent(os) << "r9: " << r9 << std::endl;
  
  indent(os) << "ecalRecHitSumEtConeDR04: " << ecalRecHitSumEtConeDR04 << std::endl;
  indent(os) << "hcalDepth1TowerSumEtConeDR04: " << hcalDepth1TowerSumEtConeDR04 << std::endl;
  indent(os) << "hcalDepth2TowerSumEtConeDR04: " << hcalDepth2TowerSumEtConeDR04 << std::endl;
  indent(os) << "hcalIsoConeDR04_2012: " << hcalIsoConeDR04_2012 << std::endl;
  indent(os) << "trkSumPtSolidConeDR04: " << trkSumPtSolidConeDR04 << std::endl;
  indent(os) << "trkSumPtHollowConeDR04: " << trkSumPtHollowConeDR04 << std::endl;
  indent(os) << "nTrkSolidConeDR04: " << int(nTrkSolidConeDR04) << std::endl;
  indent(os) << "nTrkHollowConeDR04: " << int(nTrkHollowConeDR04) << std::endl;
  
  indent(os) << "ecalRecHitSumEtConeDR03: " << ecalRecHitSumEtConeDR03 << std::endl;
  indent(os) << "hcalDepth1TowerSumEtConeDR03: " << hcalDepth1TowerSumEtConeDR03 << std::endl;
  indent(os) << "hcalDepth2TowerSumEtConeDR03: " << hcalDepth2TowerSumEtConeDR03 << std::endl;
  indent(os) << "hcalIsoConeDR03_2012: " << hcalIsoConeDR03_2012 << std::endl;
  indent(os) << "trkSumPtSolidConeDR03: " << trkSumPtSolidConeDR03 << std::endl;
  indent(os) << "trkSumPtHollowConeDR03: " << trkSumPtHollowConeDR03 << std::endl;
  indent(os) << "nTrkSolidConeDR03: " << int(nTrkSolidConeDR03) << std::endl;
  indent(os) << "nTrkHollowConeDR03: " << int(nTrkHollowConeDR03) << std::endl;
  
  indent(os) << "chargedHadronIso: " << chargedHadronIso << std::endl;
  indent(os) << "neutralHadronIso: " << neutralHadronIso << std::endl;
  indent(os) << "photonIso: " << photonIso << std::endl;
  
  indent(os) << "worstOtherVtxChargedHadronIso: " << worstOtherVtxChargedHadronIso << std::endl;
  indent(os) << "worstOtherVtxChargedHadronIsoVtxIdx: " << worstOtherVtxChargedHadronIsoVtxIdx << std::endl;
  
  indent(os) << "chargedHadronIsoDeposit: " << chargedHadronIsoDeposit << std::endl;
  indent(os) << "neutralHadronIsoDeposit: " << neutralHadronIsoDeposit << std::endl;
  indent(os) << "photonIsoDeposit: " << photonIsoDeposit << std::endl;
  
  indent(os) << "seedTime: " << seedTime << std::endl;
  
  indent(os) << "mipChi2: " << mipChi2 << std::endl;
  indent(os) << "mipTotEnergy: " << mipTotEnergy << std::endl;
  indent(os) << "mipSlope: " << mipSlope << std::endl;
  indent(os) << "mipIntercept: " << mipIntercept << std::endl;
  indent(os) << "mipNhitCone: " << mipNhitCone << std::endl;
  indent(os) << "mipIsHalo: " << mipIsHalo << std::endl;
  
  indent(os) << "convInfo: " << convInfo << std::endl;
  indent(os) << "convDist: " << convDist << std::endl;
  indent(os) << "convDcot: " << convDcot << std::endl;
  indent(os) << "convVtxChi2: " << convVtxChi2 << std::endl;
  indent(os) << "convVtxNdof: " << convVtxNdof << std::endl;
  indent(os) << "convVertex: " << convVertex << std::endl;
  indent(os) << "convDxy: " << convDxy << std::endl;
  indent(os) << "convDz: " << convDz << std::endl;
  indent(os) << "convLxy: " << convLxy << std::endl;
  indent(os) << "convLz: " << convLz << std::endl;
  indent(os) << "convZofPVFromTracks: " << convZofPVFromTracks << std::endl;
  indent(os) << "convTrackChargeProd: " << convTrackChargeProd << std::endl;
  indent(os) << "convTrack1nHit: " << convTrack1nHit << std::endl;
  indent(os) << "convTrack2nHit: " << convTrack2nHit << std::endl;
  indent(os) << "convTrack1chi2: " << convTrack1chi2 << std::endl;
  indent(os) << "convTrack2chi2: " << convTrack2chi2 << std::endl;
  indent(os) << "convTrack1pT: " << convTrack1pT << std::endl;
  indent(os) << "convTrack2pT: " << convTrack2pT << std::endl;
  indent(os) << "convTrack1InnerZ: " << convTrack1InnerZ << std::endl;
  indent(os) << "convTrack2InnerZ: " << convTrack2InnerZ << std::endl;
  indent(os) << "convTrack1InnerX: " << convTrack1InnerX << std::endl;
  indent(os) << "convTrack2InnerX: " << convTrack2InnerX << std::endl;
  indent(os) << "convTrack1InnerY: " << convTrack1InnerY << std::endl;
  indent(os) << "convTrack2InnerY: " << convTrack2InnerY << std::endl;
  indent(os) << "convTrack1Signedd0: " << convTrack1Signedd0 << std::endl;
  indent(os) << "convTrack2Signedd0: " << convTrack2Signedd0 << std::endl;
  
  indent(os) << "superClusterIndex: " << superClusterIndex << std::endl;
  indent(os) << "superClusterPreshowerEnergy: " << superClusterPreshowerEnergy << std::endl;
  indent(os) << "superClusterPhiWidth: " << superClusterPhiWidth << std::endl;
  indent(os) << "superClusterEtaWidth: " << superClusterEtaWidth << std::endl;
  indent(os) << "caloPosition: " << caloPosition << std::endl;
  
  indent(os) << "MVAregEnergyAndErr: " << MVAregEnergyAndErr.first << ", " << MVAregEnergyAndErr.second << std::endl;
  indent(os) << "MVAcorrMomentum: " << MVAcorrMomentum << std::endl;
  
  indent(os) << "momentum: " << momentum << std::endl;
}

void
susy::Photon::fillRefs(Event const* _evt)
{
  superCluster = superClusterIndex != -1 ? &_evt->superClusters[superClusterIndex] : 0;

  worstOtherVtxChargedHadronIsoVtx = worstOtherVtxChargedHadronIsoVtxIdx != -1 ? &_evt->vertices[worstOtherVtxChargedHadronIsoVtxIdx] : 0;
}

void
susy::Electron::Init()
{
  fidBit                           = 0;
  scPixCharge                      = 0;
  boolPack                         = 0;
  convFlag                         = -1;

  eSuperClusterOverP               = 0;
  eSeedClusterOverP                = 0;
  eSeedClusterOverPout             = 0;
  eEleClusterOverPout              = 0;
  deltaEtaSuperClusterTrackAtVtx   = 0;
  deltaEtaSeedClusterTrackAtCalo   = 0;
  deltaEtaEleClusterTrackAtCalo    = 0;
  deltaPhiSuperClusterTrackAtVtx   = 0;
  deltaPhiSeedClusterTrackAtCalo   = 0;
  deltaPhiEleClusterTrackAtCalo    = 0;

  shFracInnerHits                  = 0;

  sigmaEtaEta                      = 0;
  sigmaIetaIeta                    = 0;
  sigmaIphiIphi                    = 0;
  e1x5                             = 0;
  e2x5Max                          = 0;
  e5x5                             = 0;
  r9                               = 0;
  hcalDepth1OverEcal               = 0;
  hcalDepth2OverEcal               = 0;
  hcalOverEcalBc                   = 0;

  dr03TkSumPt                      = 0;
  dr03EcalRecHitSumEt              = 0;
  dr03HcalDepth1TowerSumEt         = 0;
  dr03HcalDepth2TowerSumEt         = 0;
  dr03HcalDepth1TowerSumEtBc       = 0;
  dr03HcalDepth2TowerSumEtBc       = 0;

  dr04TkSumPt                      = 0;
  dr04EcalRecHitSumEt              = 0;
  dr04HcalDepth1TowerSumEt         = 0;
  dr04HcalDepth2TowerSumEt         = 0;
  dr04HcalDepth1TowerSumEtBc       = 0;
  dr04HcalDepth2TowerSumEtBc       = 0;

  convDist                         = 0.;
  convDcot                         = 0.;
  convRadius                       = 0.;

  chargedHadronIso                 = BIGVALUE;
  neutralHadronIso                 = BIGVALUE;
  photonIso                        = BIGVALUE;

  mvaStatus                        = 0;
  mva                              = 0;

  mvaTrig                          = BIGVALUE;
  mvaNonTrig                       = BIGVALUE;

  bremClass                        = 0;
  fbrem                            = 0;

  ecalEnergy                       = 0;
  ecalEnergyError                  = 0;
  trackMomentumError               = 0;

  gsfTrackIndex                    = -1;
  closestCtfTrackIndex             = -1;
  electronClusterIndex             = -1;
  superClusterIndex                = -1;

  nMissingHits                     = 999;
  passConversionVeto               = 0;

  trackPositionAtVtx               *= 0;
  trackPositionAtCalo              *= 0;
  trackMomentumAtVtx               *= 0;
  trackMomentumAtCalo              *= 0;
  trackMomentumOut                 *= 0;
  trackMomentumAtEleClus           *= 0;
  trackMomentumAtVtxWithConstraint *= 0;

  vertex     *= 0;
  momentum   *= 0;

  gsfTrack = 0;
  closestCtfTrack = 0;
  electronCluster = 0;
  superCluster = 0;
}

void
susy::Electron::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "fidBit: " << bin(fidBit) << std::endl;
  indent(os) << "scPixCharge: " << int(scPixCharge) << std::endl;
  indent(os) << "boolPack: " << bin(boolPack) << std::endl;
  indent(os) << "convFlag: " << convFlag << std::endl;
  
  indent(os) << "eSuperClusterOverP: " << eSuperClusterOverP << std::endl;
  indent(os) << "eSeedClusterOverP: " << eSeedClusterOverP << std::endl;
  indent(os) << "eSeedClusterOverPout: " << eSeedClusterOverPout << std::endl;
  indent(os) << "eEleClusterOverPout: " << eEleClusterOverPout << std::endl;
  indent(os) << "deltaEtaSuperClusterTrackAtVtx: " << deltaEtaSuperClusterTrackAtVtx << std::endl;
  indent(os) << "deltaEtaSeedClusterTrackAtCalo: " << deltaEtaSeedClusterTrackAtCalo << std::endl;
  indent(os) << "deltaEtaEleClusterTrackAtCalo: " << deltaEtaEleClusterTrackAtCalo << std::endl;
  indent(os) << "deltaPhiSuperClusterTrackAtVtx: " << deltaPhiSuperClusterTrackAtVtx << std::endl;
  indent(os) << "deltaPhiSeedClusterTrackAtCalo: " << deltaPhiSeedClusterTrackAtCalo << std::endl;
  indent(os) << "deltaPhiEleClusterTrackAtCalo: " << deltaPhiEleClusterTrackAtCalo << std::endl;
  
  indent(os) << "shFracInnerHits: " << shFracInnerHits << std::endl;
  
  indent(os) << "sigmaEtaEta: " << sigmaEtaEta << std::endl;
  indent(os) << "sigmaIetaIeta: " << sigmaIetaIeta << std::endl;
  indent(os) << "sigmaIphiIphi: " << sigmaIphiIphi << std::endl;
  indent(os) << "e1x5: " << e1x5 << std::endl;
  indent(os) << "e2x5Max: " << e2x5Max << std::endl;
  indent(os) << "e5x5: " << e5x5 << std::endl;
  indent(os) << "r9: " << r9 << std::endl;
  indent(os) << "hcalDepth1OverEcal: " << hcalDepth1OverEcal << std::endl;
  indent(os) << "hcalDepth2OverEcal: " << hcalDepth2OverEcal << std::endl;
  indent(os) << "hcalOverEcalBc: " << hcalOverEcalBc << std::endl;
  
  indent(os) << "dr03TkSumPt: " << dr03TkSumPt << std::endl;
  indent(os) << "dr03EcalRecHitSumEt: " << dr03EcalRecHitSumEt << std::endl;
  indent(os) << "dr03HcalDepth1TowerSumEt: " << dr03HcalDepth1TowerSumEt << std::endl;
  indent(os) << "dr03HcalDepth2TowerSumEt: " << dr03HcalDepth2TowerSumEt << std::endl;
  indent(os) << "dr03HcalDepth1TowerSumEtBc: " << dr03HcalDepth1TowerSumEtBc << std::endl;
  indent(os) << "dr03HcalDepth2TowerSumEtBc: " << dr03HcalDepth2TowerSumEtBc << std::endl;
  
  indent(os) << "dr04TkSumPt: " << dr04TkSumPt << std::endl;
  indent(os) << "dr04EcalRecHitSumEt: " << dr04EcalRecHitSumEt << std::endl;
  indent(os) << "dr04HcalDepth1TowerSumEt: " << dr04HcalDepth1TowerSumEt << std::endl;
  indent(os) << "dr04HcalDepth2TowerSumEt: " << dr04HcalDepth2TowerSumEt << std::endl;
  indent(os) << "dr04HcalDepth1TowerSumEtBc: " << dr04HcalDepth1TowerSumEtBc << std::endl;
  indent(os) << "dr04HcalDepth2TowerSumEtBc: " << dr04HcalDepth2TowerSumEtBc << std::endl;
  
  indent(os) << "convDist: " << convDist << std::endl;
  indent(os) << "convDcot: " << convDcot << std::endl;
  indent(os) << "convRadius: " << convRadius << std::endl;
  
  indent(os) << "chargedHadronIso: " << chargedHadronIso << std::endl;
  indent(os) << "neutralHadronIso: " << neutralHadronIso << std::endl;
  indent(os) << "photonIso: " << photonIso << std::endl;
  
  indent(os) << "mvaStatus: " << mvaStatus << std::endl;
  indent(os) << "mva: " << mva << std::endl;
  
  indent(os) << "mvaTrig: " << mvaTrig << std::endl;
  indent(os) << "mvaNonTrig: " << mvaNonTrig << std::endl;
  
  indent(os) << "bremClass: " << int(bremClass) << std::endl;
  indent(os) << "fbrem: " << fbrem << std::endl;
  
  indent(os) << "ecalEnergy: " << ecalEnergy << std::endl;
  indent(os) << "ecalEnergyError: " << ecalEnergyError << std::endl;
  indent(os) << "trackMomentumError: " << trackMomentumError << std::endl;
  
  indent(os) << "gsfTrackIndex: " << gsfTrackIndex << std::endl;
  indent(os) << "closestCtfTrackIndex: " << closestCtfTrackIndex << std::endl;
  indent(os) << "electronClusterIndex: " << electronClusterIndex << std::endl;
  indent(os) << "superClusterIndex: " << superClusterIndex << std::endl;
  
  indent(os) << "nMissingHits: " << nMissingHits << std::endl;
  indent(os) << "passConversionVeto: " << passConversionVeto << std::endl;
  
  indent(os) << "trackPositionAtVtx: " << trackPositionAtVtx << std::endl;
  indent(os) << "trackPositionAtCalo: " << trackPositionAtCalo << std::endl;
  indent(os) << "trackMomentumAtVtx: " << trackMomentumAtVtx << std::endl;
  indent(os) << "trackMomentumAtCalo: " << trackMomentumAtCalo << std::endl;
  indent(os) << "trackMomentumOut: " << trackMomentumOut << std::endl;
  indent(os) << "trackMomentumAtEleClus: " << trackMomentumAtEleClus << std::endl;
  indent(os) << "trackMomentumAtVtxWithConstraint: " << trackMomentumAtVtxWithConstraint << std::endl;
  
  indent(os) << "vertex: " << vertex << std::endl;
  indent(os) << "momentum: " << momentum << std::endl;
}

void
susy::Electron::fillRefs(Event const* _evt)
{
  gsfTrack = gsfTrackIndex != -1 ? &_evt->tracks[gsfTrackIndex] : 0;
  closestCtfTrack = closestCtfTrackIndex != -1 ? &_evt->tracks[closestCtfTrackIndex] : 0;
  electronCluster = electronClusterIndex != -1 ? &_evt->clusters[electronClusterIndex] : 0;
  superCluster = superClusterIndex != -1 ? &_evt->superClusters[superClusterIndex] : 0;
}

void
susy::Muon::Init()
{
  type                    = 0;
  bestTrackType           = 0;
  highPtBestTrackType     = 0;

  qualityFlags            = 0;

  nChambers               = 0;
  nMatches                = 0;
  stationMask             = 0;
  nMatchedStations        = 0;
  nValidHits              = 0;
  nValidTrackerHits       = 0;
  nValidMuonHits          = 0;
  nPixelLayersWithMeasurement = 0;
  nStripLayersWithMeasurement = 0;

  timeNDof                = 0;
  timeDirection           = 0;
  timeAtIp                = 0;
  timeAtIpError           = 0;
  caloCompatibility       = 0;
  segmentCompatibility    = 0;
  emEnergy                = 0;
  hadEnergy               = 0;
  trackIsoR03             = 0;
  ecalIsoR03              = 0;
  hcalIsoR03              = 0;
  trackIsoR05             = 0;
  ecalIsoR05              = 0;
  hcalIsoR05              = 0;

  sumChargedHadronPt03    = 0;
  sumChargedParticlePt03  = 0;
  sumNeutralHadronEt03    = 0;
  sumPhotonEt03           = 0;
  sumNeutralHadronEtHighThreshold03 = 0;
  sumPhotonEtHighThreshold03        = 0;
  sumPUPt03               = 0;

  sumChargedHadronPt04    = 0;
  sumChargedParticlePt04  = 0;
  sumNeutralHadronEt04    = 0;
  sumPhotonEt04           = 0;
  sumNeutralHadronEtHighThreshold04 = 0;
  sumPhotonEtHighThreshold04        = 0;
  sumPUPt04               = 0;

  trackIndex              = -1;
  standAloneTrackIndex    = -1;
  combinedTrackIndex      = -1;
  tpfmsTrackIndex         = -1;
  pickyTrackIndex         = -1;
  dytTrackIndex           = -1;

  momentum               *= 0;

  innerTrack = 0;
  outerTrack = 0;
  globalTrack = 0;
  tpfmsTrack = 0;
  pickyTrack = 0;
  dytTrack = 0;
  bestTrack = 0;
  highPtBestTrack = 0;
}

void
susy::Muon::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "type: " << int(type) << std::endl;
  indent(os) << "bestTrackType: " << int(bestTrackType) << std::endl;
  indent(os) << "highPtBestTrackType: " << int(highPtBestTrackType) << std::endl;
  
  indent(os) << "qualityFlags: " << bin(qualityFlags) << std::endl;
  
  indent(os) << "nChambers: " << int(nChambers) << std::endl;
  indent(os) << "nMatches: " << int(nMatches) << std::endl;
  indent(os) << "stationMask: " << int(stationMask) << std::endl;
  indent(os) << "nMatchedStations: " << int(nMatchedStations) << std::endl;
  indent(os) << "nValidHits: " << int(nValidHits) << std::endl;
  indent(os) << "nValidTrackerHits: " << int(nValidTrackerHits) << std::endl;
  indent(os) << "nValidMuonHits: " << int(nValidMuonHits) << std::endl;
  indent(os) << "nPixelLayersWithMeasurement: " << int(nPixelLayersWithMeasurement) << std::endl;
  indent(os) << "nStripLayersWithMeasurement: " << int(nStripLayersWithMeasurement) << std::endl;
  
  indent(os) << "timeNDof: " << int(timeNDof) << std::endl;
  indent(os) << "timeDirection: " << int(timeDirection) << std::endl;
  indent(os) << "timeAtIp: " << timeAtIp << std::endl;
  indent(os) << "timeAtIpError: " << timeAtIpError << std::endl;
  indent(os) << "caloCompatibility: " << caloCompatibility << std::endl;
  indent(os) << "segmentCompatibility: " << segmentCompatibility << std::endl;
  indent(os) << "emEnergy: " << emEnergy << std::endl;
  indent(os) << "hadEnergy: " << hadEnergy << std::endl;
  indent(os) << "trackIsoR03: " << trackIsoR03 << std::endl;
  indent(os) << "ecalIsoR03: " << ecalIsoR03 << std::endl;
  indent(os) << "hcalIsoR03: " << hcalIsoR03 << std::endl;
  indent(os) << "trackIsoR05: " << trackIsoR05 << std::endl;
  indent(os) << "ecalIsoR05: " << ecalIsoR05 << std::endl;
  indent(os) << "hcalIsoR05: " << hcalIsoR05 << std::endl;
  
  indent(os) << "sumChargedHadronPt03: " << sumChargedHadronPt03 << std::endl;
  indent(os) << "sumChargedParticlePt03: " << sumChargedParticlePt03 << std::endl;
  indent(os) << "sumNeutralHadronEt03: " << sumNeutralHadronEt03 << std::endl;
  indent(os) << "sumPhotonEt03: " << sumPhotonEt03 << std::endl;
  indent(os) << "sumNeutralHadronEtHighThreshold03: " << sumNeutralHadronEtHighThreshold03 << std::endl;
  indent(os) << "sumPhotonEtHighThreshold03: " << sumPhotonEtHighThreshold03 << std::endl;
  indent(os) << "sumPUPt03: " << sumPUPt03 << std::endl;
  
  indent(os) << "sumChargedHadronPt04: " << sumChargedHadronPt04 << std::endl;
  indent(os) << "sumChargedParticlePt04: " << sumChargedParticlePt04 << std::endl;
  indent(os) << "sumNeutralHadronEt04: " << sumNeutralHadronEt04 << std::endl;
  indent(os) << "sumPhotonEt04: " << sumPhotonEt04 << std::endl;
  indent(os) << "sumNeutralHadronEtHighThreshold04: " << sumNeutralHadronEtHighThreshold04 << std::endl;
  indent(os) << "sumPhotonEtHighThreshold04: " << sumPhotonEtHighThreshold04 << std::endl;
  indent(os) << "sumPUPt04: " << sumPUPt04 << std::endl;
  
  indent(os) << "trackIndex: " << trackIndex << std::endl;
  indent(os) << "standAloneTrackIndex: " << standAloneTrackIndex << std::endl;
  indent(os) << "combinedTrackIndex: " << combinedTrackIndex << std::endl;
  indent(os) << "tpfmsTrackIndex: " << tpfmsTrackIndex << std::endl;
  indent(os) << "pickyTrackIndex: " << pickyTrackIndex << std::endl;
  indent(os) << "dytTrackIndex: " << dytTrackIndex << std::endl;
  
  indent(os) << "momentum: " << momentum << std::endl;
}

void
susy::Muon::fillRefs(Event const* _evt)
{
  innerTrack = trackIndex != -1 ? &_evt->tracks[trackIndex] : 0;
  outerTrack = standAloneTrackIndex != -1 ? &_evt->tracks[standAloneTrackIndex] : 0;
  globalTrack = combinedTrackIndex != -1 ? &_evt->tracks[combinedTrackIndex] : 0;
  tpfmsTrack = tpfmsTrackIndex != -1 ? &_evt->tracks[tpfmsTrackIndex] : 0;
  pickyTrack = pickyTrackIndex != -1 ? &_evt->tracks[pickyTrackIndex] : 0;
  dytTrack = dytTrackIndex != -1 ? &_evt->tracks[dytTrackIndex] : 0;
  bestTrack = bestTrackIndex() != -1 ? &_evt->tracks[bestTrackIndex()] : 0;
  highPtBestTrack = highPtBestTrackIndex() != -1 ? &_evt->tracks[highPtBestTrackIndex()] : 0;
}

void
susy::CaloJet::Init()
{
  partonFlavour             = 0;
  jetCharge                 = 0;
  etaMean                   = 0;
  phiMean                   = 0;
  etaEtaMoment              = 0;
  etaPhiMoment              = 0;
  phiPhiMoment              = 0;
  maxDistance               = 0;
  jetArea                   = 0;
  pileup                    = 0;
  nPasses                   = 0;
  nConstituents             = 0;

  maxEInEmTowers            = 0;
  maxEInHadTowers           = 0;
  energyFractionHadronic    = 0;
  emEnergyFraction          = 0;
  hadEnergyInHB             = 0;
  hadEnergyInHO             = 0;
  hadEnergyInHE             = 0;
  hadEnergyInHF             = 0;
  emEnergyInEB              = 0;
  emEnergyInEE              = 0;
  emEnergyInHF              = 0;
  towersArea                = 0;
  n90                       = 0;
  n60                       = 0;

  fHPD                      = 0;
  fRBX                      = 0;
  n90Hits                   = 0;
  fSubDetector1             = 0;
  fSubDetector2             = 0;
  fSubDetector3             = 0;
  fSubDetector4             = 0;
  restrictedEMF             = 0;
  nHCALTowers               = 0;
  nECALTowers               = 0;
  approximatefHPD           = 0;
  approximatefRBX           = 0;
  hitsInN90                 = 0;
  numberOfHits2RPC          = 0;
  numberOfHits3RPC          = 0;
  numberOfHitsRPC           = 0;

  vertex                   *= 0;
  momentum                 *= 0;
  detectorP4               *= 0;

  jecScaleFactors.clear();
  jecUncertainty            = 0;
}

void
susy::CaloJet::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "partonFlavour: " << partonFlavour << std::endl;
  indent(os) << "jetCharge: " << jetCharge << std::endl;
  indent(os) << "etaMean: " << etaMean << std::endl;
  indent(os) << "phiMean: " << phiMean << std::endl;
  indent(os) << "etaEtaMoment: " << etaEtaMoment << std::endl;
  indent(os) << "etaPhiMoment: " << etaPhiMoment << std::endl;
  indent(os) << "phiPhiMoment: " << phiPhiMoment << std::endl;
  indent(os) << "maxDistance: " << maxDistance << std::endl;
  indent(os) << "jetArea: " << jetArea << std::endl;
  indent(os) << "pileup: " << pileup << std::endl;
  indent(os) << "nPasses: " << int(nPasses) << std::endl;
  indent(os) << "nConstituents: " << int(nConstituents) << std::endl;
  
  indent(os) << "maxEInEmTowers: " << maxEInEmTowers << std::endl;
  indent(os) << "maxEInHadTowers: " << maxEInHadTowers << std::endl;
  indent(os) << "energyFractionHadronic: " << energyFractionHadronic << std::endl;
  indent(os) << "emEnergyFraction: " << emEnergyFraction << std::endl;
  indent(os) << "hadEnergyInHB: " << hadEnergyInHB << std::endl;
  indent(os) << "hadEnergyInHO: " << hadEnergyInHO << std::endl;
  indent(os) << "hadEnergyInHE: " << hadEnergyInHE << std::endl;
  indent(os) << "hadEnergyInHF: " << hadEnergyInHF << std::endl;
  indent(os) << "emEnergyInEB: " << emEnergyInEB << std::endl;
  indent(os) << "emEnergyInEE: " << emEnergyInEE << std::endl;
  indent(os) << "emEnergyInHF: " << emEnergyInHF << std::endl;
  indent(os) << "towersArea: " << towersArea << std::endl;
  indent(os) << "n90: " << int(n90) << std::endl;
  indent(os) << "n60: " << int(n60) << std::endl;
  
  indent(os) << "fHPD: " << fHPD << std::endl;
  indent(os) << "fRBX: " << fRBX << std::endl;
  indent(os) << "n90Hits: " << n90Hits << std::endl;
  indent(os) << "fSubDetector1: " << fSubDetector1 << std::endl;
  indent(os) << "fSubDetector2: " << fSubDetector2 << std::endl;
  indent(os) << "fSubDetector3: " << fSubDetector3 << std::endl;
  indent(os) << "fSubDetector4: " << fSubDetector4 << std::endl;
  indent(os) << "restrictedEMF: " << restrictedEMF << std::endl;
  indent(os) << "nHCALTowers: " << int(nHCALTowers) << std::endl;
  indent(os) << "nECALTowers: " << int(nECALTowers) << std::endl;
  indent(os) << "approximatefHPD: " << approximatefHPD << std::endl;
  indent(os) << "approximatefRBX: " << approximatefRBX << std::endl;
  indent(os) << "hitsInN90: " << int(hitsInN90) << std::endl;
  indent(os) << "numberOfHits2RPC: " << int(numberOfHits2RPC) << std::endl;
  indent(os) << "numberOfHits3RPC: " << int(numberOfHits3RPC) << std::endl;
  indent(os) << "numberOfHitsRPC: " << int(numberOfHitsRPC) << std::endl;
  
  indent(os) << "vertex: " << vertex << std::endl;
  indent(os) << "momentum: " << momentum << std::endl;
  indent(os) << "detectorP4: " << detectorP4 << std::endl;

  indent(os) << "jecScaleFactors: " << jecScaleFactors << std::endl;
  indent(os) << "jecUncertainty: " << jecUncertainty << std::endl;
}

void
susy::CaloJet::fillRefs(Event const*)
{
}

void
susy::PFJet::Init()
{
  phyDefFlavour             = 0;
  algDefFlavour             = 0;
  jetCharge                 = 0;
  etaMean                   = 0;
  phiMean                   = 0;
  etaEtaMoment              = 0;
  etaPhiMoment              = 0;
  phiPhiMoment              = 0;
  maxDistance               = 0;
  jetArea                   = 0;
  pileup                    = 0;
  nPasses                   = 0;
  nConstituents             = 0;

  chargedHadronEnergy       = 0;
  neutralHadronEnergy       = 0;
  photonEnergy              = 0;
  electronEnergy            = 0;
  muonEnergy                = 0;
  HFHadronEnergy            = 0;
  HFEMEnergy                = 0;
  chargedEmEnergy           = 0;
  chargedMuEnergy           = 0;
  neutralEmEnergy           = 0;
  chargedHadronMultiplicity = 0;
  neutralHadronMultiplicity = 0;
  photonMultiplicity        = 0;
  electronMultiplicity      = 0;
  muonMultiplicity          = 0;
  HFHadronMultiplicity      = 0;
  HFEMMultiplicity          = 0;
  chargedMultiplicity       = 0;
  neutralMultiplicity       = 0;

  tracklist.clear();

  pfParticleList.clear();

  for(unsigned i(0); i != nPUJetIdAlgorithms; ++i){
    puJetIdDiscriminants[i] = 0.;
    puJetIdFlags[i] = -1;
  }

  for(unsigned i(0); i != nBTagDiscriminators; ++i)
    bTagDiscriminators[i] = 0.;

  for(unsigned i(0); i != nQGDiscriminators; ++i)
    qgDiscriminators[i] = -1.;

  momentum                 *= 0;

  jecScaleFactors.clear();
  jecUncertainty            = 0;

  tracks.clear();
  pfParticles.clear();
}

void
susy::PFJet::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "phyDefFlavour: " << phyDefFlavour << std::endl;
  indent(os) << "algDefFlavour: " << algDefFlavour << std::endl;
  indent(os) << "jetCharge: " << jetCharge << std::endl;
  indent(os) << "etaMean: " << etaMean << std::endl;
  indent(os) << "phiMean: " << phiMean << std::endl;
  indent(os) << "etaEtaMoment: " << etaEtaMoment << std::endl;
  indent(os) << "etaPhiMoment: " << etaPhiMoment << std::endl;
  indent(os) << "phiPhiMoment: " << phiPhiMoment << std::endl;
  indent(os) << "maxDistance: " << maxDistance << std::endl;
  indent(os) << "jetArea: " << jetArea << std::endl;
  indent(os) << "pileup: " << pileup << std::endl;
  indent(os) << "nPasses: " << int(nPasses) << std::endl;
  indent(os) << "nConstituents: " << int(nConstituents) << std::endl;
  
  indent(os) << "chargedHadronEnergy: " << chargedHadronEnergy << std::endl;
  indent(os) << "neutralHadronEnergy: " << neutralHadronEnergy << std::endl;
  indent(os) << "photonEnergy: " << photonEnergy << std::endl;
  indent(os) << "electronEnergy: " << electronEnergy << std::endl;
  indent(os) << "muonEnergy: " << muonEnergy << std::endl;
  indent(os) << "HFHadronEnergy: " << HFHadronEnergy << std::endl;
  indent(os) << "HFEMEnergy: " << HFEMEnergy << std::endl;
  indent(os) << "chargedEmEnergy: " << chargedEmEnergy << std::endl;
  indent(os) << "chargedMuEnergy: " << chargedMuEnergy << std::endl;
  indent(os) << "neutralEmEnergy: " << neutralEmEnergy << std::endl;
  indent(os) << "chargedHadronMultiplicity: " << int(chargedHadronMultiplicity) << std::endl;
  indent(os) << "neutralHadronMultiplicity: " << int(neutralHadronMultiplicity) << std::endl;
  indent(os) << "photonMultiplicity: " << int(photonMultiplicity) << std::endl;
  indent(os) << "electronMultiplicity: " << int(electronMultiplicity) << std::endl;
  indent(os) << "muonMultiplicity: " << int(muonMultiplicity) << std::endl;
  indent(os) << "HFHadronMultiplicity: " << int(HFHadronMultiplicity) << std::endl;
  indent(os) << "HFEMMultiplicity: " << int(HFEMMultiplicity) << std::endl;
  indent(os) << "chargedMultiplicity: " << int(chargedMultiplicity) << std::endl;
  indent(os) << "neutralMultiplicity: " << int(neutralMultiplicity) << std::endl;
  
  indent(os) << "tracklist: " << tracklist << std::endl;
  indent(os) << "pfParticleList: " << pfParticleList << std::endl;

  indent(os) << "puJetIdDiscriminants: ";
  for(unsigned i(0); i != nPUJetIdAlgorithms; ++i) os << i << ". " << puJetIdDiscriminants[i] << " ";
  os << std::endl;
  indent(os) << "puJetIdFlags: ";
  for(unsigned i(0); i != nPUJetIdAlgorithms; ++i) os << i << ". " << puJetIdFlags[i] << " ";
  os << std::endl;

  indent(os) << "bTagDiscriminators: ";
  for(unsigned i(0); i != nBTagDiscriminators; ++i) os << i << ". " << bTagDiscriminators[i] << " ";
  os << std::endl;

  indent(os) << "qgDiscriminators: ";
  for(unsigned i(0); i != nQGDiscriminators; ++i) os << i << ". " << qgDiscriminators[i] << " ";
  os << std::endl;
  
  indent(os) << "momentum: " << momentum << std::endl;
  
  indent(os) << "jecScaleFactors: " << jecScaleFactors << std::endl;
  indent(os) << "jecUncertainty: " << jecUncertainty << std::endl;
}

void
susy::PFJet::fillRefs(Event const* _evt)
{
  tracks.assign(tracklist.size(), 0);
  for(unsigned iT(0); iT != tracklist.size(); ++iT)
    tracks[iT] = &_evt->tracks[tracklist[iT]];

  pfParticles.assign(pfParticleList.size(), 0);
  for(unsigned iP(0); iP != pfParticleList.size(); ++iP)
    pfParticles[iP] = &_evt->pfParticles[pfParticleList[iP]];
}

void
susy::JPTJet::Init()
{
  partonFlavour             = 0;
  jetCharge                 = 0;
  etaMean                   = 0;
  phiMean                   = 0;
  etaEtaMoment              = 0;
  etaPhiMoment              = 0;
  phiPhiMoment              = 0;
  maxDistance               = 0;
  jetArea                   = 0;
  pileup                    = 0;
  nPasses                   = 0;
  nConstituents             = 0;

  chargedHadronEnergy       = 0;
  neutralHadronEnergy       = 0;
  chargedEmEnergy           = 0;
  neutralEmEnergy           = 0;
  chargedMultiplicity       = 0;
  muonMultiplicity          = 0;
  elecMultiplicity          = 0;
  getZSPCor                 = 0;

  momentum                 *= 0;

  jecScaleFactors.clear();
  jecUncertainty            = 0;
}

void
susy::JPTJet::Print(std::ostream& os/* = std::cout*/) const
{
  indent(os) << "partonFlavour: " << partonFlavour << std::endl;
  indent(os) << "jetCharge: " << jetCharge << std::endl;
  indent(os) << "etaMean: " << etaMean << std::endl;
  indent(os) << "phiMean: " << phiMean << std::endl;
  indent(os) << "etaEtaMoment: " << etaEtaMoment << std::endl;
  indent(os) << "etaPhiMoment: " << etaPhiMoment << std::endl;
  indent(os) << "phiPhiMoment: " << phiPhiMoment << std::endl;
  indent(os) << "maxDistance: " << maxDistance << std::endl;
  indent(os) << "jetArea: " << jetArea << std::endl;
  indent(os) << "pileup: " << pileup << std::endl;
  indent(os) << "nPasses: " << int(nPasses) << std::endl;
  indent(os) << "nConstituents: " << int(nConstituents) << std::endl;
  
  indent(os) << "chargedHadronEnergy: " << chargedHadronEnergy << std::endl;
  indent(os) << "neutralHadronEnergy: " << neutralHadronEnergy << std::endl;
  indent(os) << "chargedEmEnergy: " << chargedEmEnergy << std::endl;
  indent(os) << "neutralEmEnergy: " << neutralEmEnergy << std::endl;
  indent(os) << "chargedMultiplicity: " << int(chargedMultiplicity) << std::endl;
  indent(os) << "muonMultiplicity: " << int(muonMultiplicity) << std::endl;
  indent(os) << "elecMultiplicity: " << int(elecMultiplicity) << std::endl;
  indent(os) << "getZSPCor: " << getZSPCor << std::endl;
  
  indent(os) << "momentum: " << momentum << std::endl;

  indent(os) << "jecScaleFactors: " << jecScaleFactors << std::endl;
  indent(os) << "jecUncertainty: " << jecUncertainty << std::endl;
}

void
susy::JPTJet::fillRefs(Event const*)
{
}

susy::TriggerMap::TriggerMap(TString const& _type) :
  trigType_(_type),
  core_(),
  prescales_(0),
  decisions_(0),
  menuInEvent_(new TString),
  currentMenu_(""),
  treeNumber_(-1),
  inputTree_(0),
  outputTrees_()
{
}

susy::TriggerMap::~TriggerMap()
{
  releaseTrees(kTRUE);

  clear_();

  delete menuInEvent_;
}

void
susy::TriggerMap::Print(std::ostream& os/* = std::cout*/) const
{
  os << "Trigger Table: " << currentMenu_ << std::endl;
  for(CoreConstIter itr(core_.begin()); itr != core_.end(); ++itr)
    indent(os) << itr->first << "(prescale=" << *(itr->second.first) << ", fire=" << int(*(itr->second.second)) << ")" << std::endl;
}

Bool_t
susy::TriggerMap::pass(TString const& _path) const
{
  CoreConstIter itr(core_.lower_bound(_path));
  if(itr == core_.end()) return kFALSE;

  TString const& found(itr->first);

  if((_path.First('*') == _path.Length() - 1 && found.Index(_path.SubString(0, _path.Length() - 1)) == 0) || _path == found)
    return *itr->second.second != 0;
  else
    return kFALSE;
}

UInt_t
susy::TriggerMap::prescale(TString const& _path) const
{
  CoreConstIter itr(core_.lower_bound(_path));
  if(itr == core_.end()) return kFALSE;

  TString const& found(itr->first);

  if((_path.First('*') == _path.Length() - 1 && found.Index(_path.SubString(0, _path.Length() - 1)) == 0) || _path == found)
    return *itr->second.first;
  else
    return 0;
}

std::pair<UInt_t, Bool_t>
susy::TriggerMap::operator[](TString const& _path) const
{
  std::pair<UInt_t, Bool_t> value(0, kFALSE);

  CoreConstIter itr(core_.lower_bound(_path));
  if(itr == core_.end()) return value;

  TString const& found(itr->first);

  if((_path.First('*') == _path.Length() - 1 && found.Index(_path.SubString(0, _path.Length() - 1)) == 0) || _path == found){
    value.first = *itr->second.first;
    value.second = *itr->second.second != 0;
  }
  return value;
}

void
susy::TriggerMap::setInput(TTree& _input)
{
  if(inputTree_) releaseTree(*inputTree_, kTRUE);

  inputTree_ = &_input;

  inputTree_->SetBranchAddress(trigType_ + "Config", &menuInEvent_);

  currentMenu_ = "";
  treeNumber_ = -1;
}

void
susy::TriggerMap::addOutput(TTree& _output)
{
  outputTrees_.push_back(&_output);

  TBranch* configBranch(_output.GetBranch(trigType_ + "Config"));
  if(configBranch){
    // output is a halfway-filled susyTree
    _output.SetBranchAddress(trigType_ + "Config", &menuInEvent_);
  }
  else{
    configBranch = _output.Branch(trigType_ + "Config", "TString", &menuInEvent_);
    configBranch->SetFirstEntry(_output.GetEntries());
  }

  if(currentMenu_ == "") return;

  // look for current trigger table
  TBranch* bitBranch(0);
  TBranch* psBranch(0);
  if(configBranch){
    bitBranch = _output.GetBranch(trigType_ + "Prescales_" + currentMenu_);
    psBranch = _output.GetBranch(trigType_ + "Bits_" + currentMenu_);
  }

  if(bitBranch && psBranch){
    _output.SetBranchAddress(trigType_ + "Bits_" + currentMenu_, decisions_);
    _output.SetBranchAddress(trigType_ + "Prescales_" + currentMenu_, prescales_);
  }
  else{
    TString bitLeaves(formLeafList_());
    TString psLeaves(bitLeaves);
    psLeaves.ReplaceAll("/b", "/i");

    bitBranch = _output.Branch(trigType_ + "Bits_" + currentMenu_, decisions_, bitLeaves);
    psBranch = _output.Branch(trigType_ + "Prescales_" + currentMenu_, prescales_, psLeaves);

    bitBranch->SetFirstEntry(_output.GetEntries());
    psBranch->SetFirstEntry(_output.GetEntries());
  }
}

void
susy::TriggerMap::checkInput()
{
  // Check for menu or tree transition and set the branch addresses when necessary. Note that
  // when the input is a TChain, different constituent trees can have different trigger tables.
  // As such, SetBranchAddress is not issued on the inputTree itself (otherwise users will see
  // a lot of "SetBranchAddress: No branch XYZ found" messages) but only on the current tree.
  // The check for tree transition is required for this reason.

  if(!inputTree_) return;

  if(*menuInEvent_ == currentMenu_ && treeNumber_ == inputTree_->GetTreeNumber()) return;

  TBranch* bitBranch(inputTree_->GetBranch(trigType_ + "Bits_" + *menuInEvent_));
  TBranch* psBranch(inputTree_->GetBranch(trigType_ + "Prescales_" + *menuInEvent_));

  if(!bitBranch || !psBranch){
    // Input tree is not properly formatted
    releaseTree(*inputTree_, kTRUE);
    return;
  }

  if(*menuInEvent_ != currentMenu_){
    // New trigger table
    TObjArray* leaves(bitBranch->GetListOfLeaves());

    std::vector<std::string> paths;
    for(int iL(0); iL != leaves->GetEntries(); ++iL)
      paths.push_back(leaves->At(iL)->GetName());

    setMenu(*menuInEvent_, paths); // currentMenu is set here
  }

  treeNumber_ = inputTree_->GetTreeNumber();

  bitBranch->SetAddress(decisions_);
  psBranch->SetAddress(prescales_);

  inputTree_->GetEntry(inputTree_->GetReadEntry());
}

void
susy::TriggerMap::setMenu(TString const& _menuName, std::vector<std::string> const& _paths)
{
  // Configure the object for a new trigger table. Called either by the ntuplizer or internally
  // from checkInput() when the input configuration changes.

  if(_menuName == currentMenu_) return;

  releaseTrees(kFALSE);

  clear_();

  *menuInEvent_ = _menuName;
  currentMenu_ = _menuName;

  unsigned nP(_paths.size());
  unsigned iP(0);
  for(; iP != nP; ++iP)
    core_.insert(typename MapCore::value_type(_paths[iP], std::pair<UInt_t*, UChar_t*>(0, 0)));

  prescales_ = new UInt_t[nP];
  decisions_ = new UChar_t[nP];

  iP = 0;
  for(CoreIter itr(core_.begin()); itr != core_.end(); ++itr, ++iP){
    itr->second.first = prescales_ + iP;
    itr->second.second = decisions_ + iP;
  }

  std::fill_n(prescales_, nP, 0);
  std::fill_n(decisions_, nP, 0);

  if(outputTrees_.size() == 0) return;

  // If outputs exist, new branches must be created on each of them.

  TString bitLeaves(formLeafList_());
  TString psLeaves(bitLeaves);
  psLeaves.ReplaceAll("/b", "/i");

  for(unsigned iT(0); iT != outputTrees_.size(); ++iT){
    TTree* tree(outputTrees_[iT]);

    TBranch* bitBranch(tree->GetBranch(trigType_ + "Bits_" + currentMenu_));
    TBranch* psBranch(tree->GetBranch(trigType_ + "Prescales_" + currentMenu_));

    if(bitBranch)
      tree->SetBranchAddress(trigType_ + "Bits_" + currentMenu_, decisions_);
    else{
      bitBranch = tree->Branch(trigType_ + "Bits_" + currentMenu_, decisions_, bitLeaves);
      bitBranch->SetFirstEntry(tree->GetEntries());
    }
    if(psBranch)
      tree->SetBranchAddress(trigType_ + "Prescales_" + currentMenu_, prescales_);
    else{
      psBranch = tree->Branch(trigType_ + "Prescales_" + currentMenu_, prescales_, psLeaves);
      psBranch->SetFirstEntry(tree->GetEntries());
    }
  }
}

void
susy::TriggerMap::set(TString const& _path, UInt_t _prescale, Bool_t _decision)
{
  // Set the result of one trigger path

  CoreIter itr(core_.find(_path));
  if(itr == core_.end()) return;

  *itr->second.first = _prescale;
  *itr->second.second = _decision ? 1 : 0;
}

void
susy::TriggerMap::copy(TriggerMap const& _orig)
{
  unsigned nO(_orig.size());

  if(currentMenu_ != _orig.currentMenu_){
    std::vector<std::string> paths;
    for(CoreConstIter itr(_orig.core_.begin()); itr != _orig.core_.end(); ++itr)
      paths.push_back(itr->first.Data());
    setMenu(_orig.currentMenu_, paths);
    treeNumber_ = -1;
  }

  std::copy(_orig.prescales_, _orig.prescales_ + nO, prescales_);
  std::copy(_orig.decisions_, _orig.decisions_ + nO, decisions_);
}

void
susy::TriggerMap::releaseTree(TTree& _tree, Bool_t _fullRelease)
{
  TBranch* branch(0);
  if((branch = _tree.GetBranch(trigType_ + "Bits_" + currentMenu_)))
    _tree.ResetBranchAddress(branch);
  if((branch = _tree.GetBranch(trigType_ + "Prescales_" + currentMenu_)))
    _tree.ResetBranchAddress(branch);
  if(_fullRelease && (branch = _tree.GetBranch(trigType_ + "Config")))
    _tree.ResetBranchAddress(branch);

  if(!_fullRelease) return;

  if(inputTree_ == &_tree){
    inputTree_ = 0;
    treeNumber_ = -1;
  }
  else{
    for(unsigned iT(0); iT != outputTrees_.size(); ++iT){
      if(outputTrees_[iT] == &_tree){
	outputTrees_.erase(outputTrees_.begin() + iT);
	break;
      }
    }
  }
}

void
susy::TriggerMap::releaseTrees(Bool_t _fullRelease)
{
  if(inputTree_) releaseTree(*inputTree_, _fullRelease);

  if(_fullRelease){
    while(outputTrees_.size() != 0)
      releaseTree(*outputTrees_[0], kTRUE); // releaseTree removes the tree from the vector
  }
  else{
    for(unsigned iT(0); iT != outputTrees_.size(); ++iT)
      releaseTree(*outputTrees_[iT], kFALSE);
  }
}

void
susy::TriggerMap::clear_()
{
  core_.clear();

  delete [] prescales_;
  delete [] decisions_;
  prescales_ = 0;
  decisions_ = 0;

  currentMenu_ = "";
  treeNumber_ = -1;
}

TString
susy::TriggerMap::formLeafList_() const
{
  TString leaves;
  for(CoreConstIter itr(core_.begin()); itr != core_.end(); ++itr)
    leaves += itr->first + "/b:";

  leaves.Remove(leaves.Length() - 1);

  return leaves;
}

susy::Event::Event() :
  l1Map("l1"),
  hltMap("hlt"),
  inputTree_(0),
  outputTrees_()
{
  Init();

  metFilterMask = 
    (1 << kCSCBeamHalo) | (1 << kHcalNoise) | (1 << kEcalDeadCellTP) | (1 << kHcalLaserRECOUserStep) | (1 << kTrackingFailure) |
    (1 << kEEBadSC) | (1 << kLogErrorTooManyClusters) | (1 << kLogErrorTooManyTripletsPairs) | (1 << kLogErrorTooManySeeds);
}

susy::Event::~Event()
{
  releaseTrees();
}

void
susy::Event::Init()
{
  isRealData                  = kFALSE;
  runNumber                   = 0;
  eventNumber                 = 0;
  luminosityBlockNumber       = 0;
  bunchCrossing               = 0;
  avgInsRecLumi               = 0;
  intgRecLumi                 = 0;
  cosmicFlag                  = 0;
  rho                         = 0;
  rhoBarrel                   = 0;
  rho25                       = 0;
  metFilterBit                = 0xffffffff;

  beamSpot                   *= 0;

  vertices.clear();
  tracks.clear();
  superClusters.clear();
  clusters.clear();
  pfParticles.clear();
  pu.clear();
  genParticles.clear();

  for(std::map<TString, MET>::iterator itr(metMap.begin()); itr != metMap.end(); ++itr) itr->second.Init();
  for(std::map<TString, MuonCollection>::iterator itr(muons.begin()); itr != muons.end(); ++itr) itr->second.clear();
  for(std::map<TString, ElectronCollection>::iterator itr(electrons.begin()); itr != electrons.end(); ++itr) itr->second.clear();
  for(std::map<TString, PhotonCollection>::iterator itr(photons.begin()); itr != photons.end(); ++itr) itr->second.clear();
  for(std::map<TString, CaloJetCollection>::iterator itr(caloJets.begin()); itr != caloJets.end(); ++itr) itr->second.clear();
  for(std::map<TString, PFJetCollection>::iterator itr(pfJets.begin()); itr != pfJets.end(); ++itr) itr->second.clear();
  for(std::map<TString, JPTJetCollection>::iterator itr(jptJets.begin()); itr != jptJets.end(); ++itr) itr->second.clear();
  for(std::map<TString, Float_t>::iterator itr(gridParams.begin()); itr != gridParams.end(); ++itr) itr->second = 0.;
}

void
susy::Event::Print(std::ostream& os/* = std::cout*/) const
{
  os << "---------- run " << runNumber << " event " << eventNumber << " ----------" << std::endl;

  os << "isRealData: " << int(isRealData) << std::endl;
  os << "runNumber: " << runNumber << std::endl;
  os << "eventNumber: " << eventNumber << std::endl;
  os << "luminosityBlockNumber: " << luminosityBlockNumber << std::endl;
  os << "bunchCrossing: " << bunchCrossing << std::endl;
  os << "avgInsRecLumi: " << avgInsRecLumi << std::endl;
  os << "intgRecLumi: " << intgRecLumi << std::endl;
  os << "cosmicFlag: " << int(cosmicFlag) << std::endl;
  os << "rho: " << rho << std::endl;
  os << "rhoBarrel: " << rhoBarrel << std::endl;
  os << "rho25: " << rho25 << std::endl;
  os << "metFilterBit: " << bin(metFilterBit) << std::endl;
  os << "metFilterBit break down ===> " << std::endl;
  indent(os) << "CSCBeamHalo(" << passMetFilter(kCSCBeamHalo) << ") ";
  indent(os) << "HcalNoise(" << passMetFilter(kHcalNoise) << ") ";
  indent(os) << "EcalDeadCellTP(" << passMetFilter(kEcalDeadCellTP) << ") ";
  indent(os) << "EcalDeadCellBE(" << passMetFilter(kEcalDeadCellBE) << ") ";
  indent(os) << "HcalLaserOccupancy(" << passMetFilter(kHcalLaserOccupancy) << ") ";
  indent(os) << "TrackingFailure(" << passMetFilter(kTrackingFailure) << ") ";
  indent(os) << "EEBadSC(" << passMetFilter(kEEBadSC) << ") ";
  indent(os) << "HcalLaserEventList(" << passMetFilter(kHcalLaserEventList) << ") ";
  indent(os) << "EcalLaserCorr(" << passMetFilter(kEcalLaserCorr) << ") ";
  indent(os) << "ManyStripClus53X(" << passMetFilter(kManyStripClus53X) << ") ";
  indent(os) << "TooManyStripClus53X(" << passMetFilter(kTooManyStripClus53X) << ") ";
  indent(os) << "LogErrorTooManyClusters(" << passMetFilter(kLogErrorTooManyClusters) << ") ";
  indent(os) << "EERingOfFire(" << passMetFilter(kEERingOfFire) << ") ";
  indent(os) << "InconsistentMuon(" << passMetFilter(kInconsistentMuon) << ") ";
  indent(os) << "GreedyMuon(" << passMetFilter(kGreedyMuon) << ") ";
  indent(os) << "HcalLaserRECOUserStep(" << passMetFilter(kHcalLaserRECOUserStep) << ") ";
  os << std::endl << std::endl;
  
  os << "beamSpot: " << beamSpot << std::endl;

  os << "l1Map size(" << l1Map.size() << ") ======>" << std::endl;
  l1Map.Print(os);
  os << std::endl;

  os << "hltMap size(" << hltMap.size() << ") ======>" << std::endl;
  hltMap.Print(os);
  os << std::endl;

  os << "vertices size(" << vertices.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != vertices.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    vertices[i].Print(os);
  }
  os << std::endl;

  os << "tracks size(" << tracks.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != tracks.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    tracks[i].Print(os);
  }
  os << std::endl;

  os << "superClusters size(" << superClusters.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != superClusters.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    superClusters[i].Print(os);
  }
  os << std::endl;

  os << "clusters size(" << clusters.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != clusters.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    clusters[i].Print(os);
  }
  os << std::endl;

  os << "pfParticles size(" << pfParticles.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != pfParticles.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    pfParticles[i].Print(os);
  }
  os << std::endl;

  os << "metMap ======>" << std::endl;
  for(std::map<TString, MET>::const_iterator it = metMap.begin(); it != metMap.end(); it++){
    os << it->first << " ===>" << std::endl;
    it->second.Print(os);
  }
  os << std::endl;

  os << "muons =========>" << std::endl;
  for(std::map<TString, MuonCollection>::const_iterator it = muons.begin(); it != muons.end(); it++) {
    os << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(unsigned i(0); i != it->second.size(); ++i){
      os << "[" << i << "] ===>" << std::endl;
      it->second[i].Print(os);
    }
  }
  os << std::endl;

  os << "electrons ======>" << std::endl;
  for(std::map<TString, ElectronCollection>::const_iterator it = electrons.begin(); it != electrons.end(); it++) {
    os << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(unsigned i(0); i != it->second.size(); ++i){
      os << "[" << i << "] ===>" << std::endl;
      it->second[i].Print(os);
    }
  }
  os << std::endl;

  os << "photons ======>" << std::endl;
  for(std::map<TString, PhotonCollection>::const_iterator it = photons.begin(); it != photons.end(); it++) {
    os << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(unsigned i(0); i != it->second.size(); ++i){
      os << "[" << i << "] ===>" << std::endl;
      it->second[i].Print(os);
    }
  }
  os << std::endl;

  os << "caloJets ======>" << std::endl;
  for(std::map<TString, CaloJetCollection>::const_iterator it = caloJets.begin(); it != caloJets.end(); it++) {
    os << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(unsigned i(0); i != it->second.size(); ++i){
      os << "[" << i << "] ===>" << std::endl;
      it->second[i].Print(os);
    }
  }
  os << std::endl;

  os << "pfJets ======>" << std::endl;
  for(std::map<TString, PFJetCollection>::const_iterator it = pfJets.begin(); it != pfJets.end(); it++) {
    os << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(unsigned i(0); i != it->second.size(); ++i){
      os << "[" << i << "] ===>" << std::endl;
      it->second[i].Print(os);
    }
  }
  os << std::endl;

  os << "jptJets ======>" << std::endl;
  for(std::map<TString, JPTJetCollection>::const_iterator it = jptJets.begin(); it != jptJets.end(); it++) {
    os << it->first << " size(" << it->second.size() << ") ======>" << std::endl;
    for(unsigned i(0); i != it->second.size(); ++i){
      os << "[" << i << "] ===>" << std::endl;
      it->second[i].Print(os);
    }
  }
  os << std::endl;

  if(isRealData == 1) return;

  os << "pu size(" << pu.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != pu.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    pu[i].Print(os);
  }
  os << std::endl;

  os << "genParticles size(" << genParticles.size() << ") ======>" << std::endl;
  for(unsigned i(0); i != genParticles.size(); ++i){
    os << "[" << i << "] ===>" << std::endl;
    genParticles[i].Print(os);
  }
  os << std::endl;

  os << "gridParams ===>" << std::endl;
  indent(os) << gridParams << std::endl;

  os << std::endl << std::endl;
}

void
susy::Event::fillRefs()
{
  for(TrackCollection::iterator tItr(tracks.begin()); tItr != tracks.end(); ++tItr)
    tItr->fillRefs(this);
  for(SuperClusterCollection::iterator sItr(superClusters.begin()); sItr != superClusters.end(); ++sItr)
    sItr->fillRefs(this);
  for(ParticleCollection::iterator pItr(genParticles.begin()); pItr != genParticles.end(); ++pItr)
    pItr->fillRefs(this);

  for(std::map<TString, MuonCollection>::iterator cItr(muons.begin()); cItr != muons.end(); ++cItr)
    for(MuonCollection::iterator pItr(cItr->second.begin()); pItr != cItr->second.end(); ++pItr)
      pItr->fillRefs(this);
  for(std::map<TString, ElectronCollection>::iterator cItr(electrons.begin()); cItr != electrons.end(); ++cItr)
    for(ElectronCollection::iterator pItr(cItr->second.begin()); pItr != cItr->second.end(); ++pItr)
      pItr->fillRefs(this);
  for(std::map<TString, PhotonCollection>::iterator cItr(photons.begin()); cItr != photons.end(); ++cItr)
    for(PhotonCollection::iterator pItr(cItr->second.begin()); pItr != cItr->second.end(); ++pItr)
      pItr->fillRefs(this);
  for(std::map<TString, CaloJetCollection>::iterator cItr(caloJets.begin()); cItr != caloJets.end(); ++cItr)
    for(CaloJetCollection::iterator pItr(cItr->second.begin()); pItr != cItr->second.end(); ++pItr)
      pItr->fillRefs(this);
  for(std::map<TString, PFJetCollection>::iterator cItr(pfJets.begin()); cItr != pfJets.end(); ++cItr)
    for(PFJetCollection::iterator pItr(cItr->second.begin()); pItr != cItr->second.end(); ++pItr)
      pItr->fillRefs(this);
  for(std::map<TString, JPTJetCollection>::iterator cItr(jptJets.begin()); cItr != jptJets.end(); ++cItr)
    for(JPTJetCollection::iterator pItr(cItr->second.begin()); pItr != cItr->second.end(); ++pItr)
      pItr->fillRefs(this);
}

void
susy::Event::setInput(TTree& _tree)
{
  // To fix the branch status - seems odd but necessary due to a feature in TChain implementation.
  // TChain will enable a branch even if "*" is set to 0, when the branch address is set before a tree is loaded.
  _tree.LoadTree(0);

  _tree.SetBranchAddress("isRealData", &isRealData);
  _tree.SetBranchAddress("runNumber", &runNumber);
  _tree.SetBranchAddress("eventNumber", &eventNumber);
  _tree.SetBranchAddress("luminosityBlockNumber", &luminosityBlockNumber);
  _tree.SetBranchAddress("bunchCrossing", &bunchCrossing);
  _tree.SetBranchAddress("avgInsRecLumi", &avgInsRecLumi);
  _tree.SetBranchAddress("intgRecLumi", &intgRecLumi);
  _tree.SetBranchAddress("cosmicFlag", &cosmicFlag);
  _tree.SetBranchAddress("rho", &rho);
  _tree.SetBranchAddress("rhoBarrel", &rhoBarrel);
  _tree.SetBranchAddress("rho25", &rho25);
  _tree.SetBranchAddress("metFilterBit", &metFilterBit);

  if(_tree.GetBranchStatus("beamSpot.")) _tree.SetBranchAddress("beamSpot.", new TVector3*(&beamSpot));
  if(_tree.GetBranchStatus("vertices")) _tree.SetBranchAddress("vertices", new VertexCollection*(&vertices));
  if(_tree.GetBranchStatus("tracks")) _tree.SetBranchAddress("tracks", new TrackCollection*(&tracks));
  if(_tree.GetBranchStatus("superClusters")) _tree.SetBranchAddress("superClusters", new SuperClusterCollection*(&superClusters));
  if(_tree.GetBranchStatus("clusters")) _tree.SetBranchAddress("clusters", new ClusterCollection*(&clusters));
  if(_tree.GetBranchStatus("pfParticles")) _tree.SetBranchAddress("pfParticles", new PFParticleCollection*(&pfParticles));
  if(_tree.GetBranchStatus("pu")) _tree.SetBranchAddress("pu", new PUSummaryInfoCollection*(&pu));
  if(_tree.GetBranchStatus("genParticles")) _tree.SetBranchAddress("genParticles", new ParticleCollection*(&genParticles));

  metMap.clear();
  muons.clear();
  electrons.clear();
  photons.clear();
  caloJets.clear();
  pfJets.clear();
  jptJets.clear();
  gridParams.clear();

  TObjArray* branches(_tree.GetListOfBranches());

  for(int iB(0); iB != branches->GetEntries(); ++iB){
    TString bName(branches->At(iB)->GetName());
    if(!_tree.GetBranchStatus(bName)) continue;

    TString collectionName(bName(bName.First('_') + 1, bName.Length()));
    void** add(0);

    if(bName.Index("met_") == 0)
      add = reinterpret_cast<void**>(new MET*(&metMap[collectionName(0, collectionName.Length() - 1)]));
    else if(bName.Index("muons_") == 0)
      add = reinterpret_cast<void**>(new MuonCollection*(&muons[collectionName]));
    else if(bName.Index("electrons_") == 0)
      add = reinterpret_cast<void**>(new ElectronCollection*(&electrons[collectionName]));
    else if(bName.Index("photons_") == 0)
      add = reinterpret_cast<void**>(new PhotonCollection*(&photons[collectionName]));
    else if(bName.Index("caloJets_") == 0)
      add = reinterpret_cast<void**>(new CaloJetCollection*(&caloJets[collectionName]));
    else if(bName.Index("pfJets_") == 0)
      add = reinterpret_cast<void**>(new PFJetCollection*(&pfJets[collectionName]));
    else if(bName.Index("jptJets_") == 0)
      add = reinterpret_cast<void**>(new JPTJetCollection*(&jptJets[collectionName]));
    else if(bName.Index("gridParams_") == 0){
      _tree.SetBranchAddress(bName, &gridParams[collectionName]);
      continue;
    }
    else
      continue;

    _tree.SetBranchAddress(bName, add);
  }

  if(inputTree_) releaseTree(*inputTree_);
  inputTree_ = &_tree;
  hltMap.setInput(_tree);
  l1Map.setInput(_tree);
}

void
susy::Event::addOutput(TTree& _tree)
{
  _tree.Branch("isRealData", &isRealData, "isRealData/b");
  _tree.Branch("runNumber", &runNumber, "runNumber/i");
  _tree.Branch("eventNumber", &eventNumber, "eventNumber/i");
  _tree.Branch("luminosityBlockNumber", &luminosityBlockNumber, "luminosityBlockNumber/i");
  _tree.Branch("bunchCrossing", &bunchCrossing, "bunchCrossing/s");
  _tree.Branch("cosmicFlag", &cosmicFlag, "cosmicFlag/b");
  _tree.Branch("avgInsRecLumi", &avgInsRecLumi, "avgInsRecLumi/F");
  _tree.Branch("intgRecLumi", &intgRecLumi, "intgRecLumi/F");
  _tree.Branch("rho", &rho, "rho/F");
  _tree.Branch("rhoBarrel", &rhoBarrel, "rhoBarrel/F");
  _tree.Branch("rho25", &rho25, "rho25/F");
  _tree.Branch("metFilterBit", &metFilterBit, "metFilterBit/I");

  _tree.Branch("beamSpot.", "TVector3", new TVector3*(&beamSpot));
  _tree.Branch("vertices", "std::vector<susy::Vertex>", new VertexCollection*(&vertices));
  _tree.Branch("tracks", "std::vector<susy::Track>", new TrackCollection*(&tracks));
  _tree.Branch("superClusters", "std::vector<susy::SuperCluster>", new SuperClusterCollection*(&superClusters));
  _tree.Branch("clusters", "std::vector<susy::Cluster>", new ClusterCollection*(&clusters));
  _tree.Branch("pfParticles", "std::vector<susy::PFParticle>", new PFParticleCollection*(&pfParticles));
  _tree.Branch("pu", "std::vector<susy::PUSummaryInfo>", new PUSummaryInfoCollection*(&pu));
  _tree.Branch("genParticles", "std::vector<susy::Particle>", new ParticleCollection*(&genParticles));

  for(std::map<TString, MET>::iterator itr(metMap.begin()); itr != metMap.end(); ++itr)
    _tree.Branch("met_" + itr->first + ".", "susy::MET", new MET*(&itr->second));
  for(std::map<TString, MuonCollection>::iterator itr(muons.begin()); itr != muons.end(); ++itr)
    _tree.Branch("muons_" + itr->first, "std::vector<susy::Muon>", new MuonCollection*(&itr->second));
  for(std::map<TString, ElectronCollection>::iterator itr(electrons.begin()); itr != electrons.end(); ++itr)
    _tree.Branch("electrons_" + itr->first, "std::vector<susy::Electron>", new ElectronCollection*(&itr->second));
  for(std::map<TString, PhotonCollection>::iterator itr(photons.begin()); itr != photons.end(); ++itr)
    _tree.Branch("photons_" + itr->first, "std::vector<susy::Photon>", new PhotonCollection*(&itr->second));
  for(std::map<TString, CaloJetCollection>::iterator itr(caloJets.begin()); itr != caloJets.end(); ++itr)
    _tree.Branch("caloJets_" + itr->first, "std::vector<susy::CaloJet>", new CaloJetCollection*(&itr->second));
  for(std::map<TString, PFJetCollection>::iterator itr(pfJets.begin()); itr != pfJets.end(); ++itr)
    _tree.Branch("pfJets_" + itr->first, "std::vector<susy::PFJet>", new PFJetCollection*(&itr->second));
  for(std::map<TString, JPTJetCollection>::iterator itr(jptJets.begin()); itr != jptJets.end(); ++itr)
    _tree.Branch("jptJets_" + itr->first, "std::vector<susy::JPTJet>", new JPTJetCollection*(&itr->second));
  for(std::map<TString, Float_t>::iterator itr(gridParams.begin()); itr != gridParams.end(); ++itr)
    _tree.Branch("gridParams_" + itr->first, &itr->second, itr->first + "/F");

  outputTrees_.push_back(&_tree);
  hltMap.addOutput(_tree);
  l1Map.addOutput(_tree);
}

Int_t
susy::Event::getEntry(Long64_t _iEntry)
{
  if(!inputTree_) return 0;

  Long64_t result(inputTree_->GetEntry(_iEntry));

  if(result > 0){
    l1Map.checkInput();
    hltMap.checkInput();
  }

  fillRefs();

  return result;
}

void
susy::Event::releaseTree(TTree& _tree)
{
  int iTree(-2);
  if(&_tree == inputTree_) iTree = -1;
  else{
    for(unsigned iT(0); iT != outputTrees_.size(); ++iT){
      if(outputTrees_[iT] == &_tree){
	iTree = iT;
	break;
      }
    }
  }
  if(iTree == -2) return;

  TVector3** beamSpotP(_tree.GetBranchStatus("beamSpot.") ? reinterpret_cast<TVector3**>(_tree.FindBranch("beamSpot.")->GetAddress()) : 0);
  VertexCollection** verticesP(_tree.GetBranchStatus("vertices") ? reinterpret_cast<VertexCollection**>(_tree.FindBranch("vertices")->GetAddress()) : 0);
  TrackCollection** tracksP(_tree.GetBranchStatus("tracks") ? reinterpret_cast<TrackCollection**>(_tree.FindBranch("tracks")->GetAddress()) : 0);
  SuperClusterCollection** superClustersP(_tree.GetBranchStatus("superClusters") ? reinterpret_cast<SuperClusterCollection**>(_tree.FindBranch("superClusters")->GetAddress()) : 0);
  ClusterCollection** clustersP(_tree.GetBranchStatus("clusters") ? reinterpret_cast<ClusterCollection**>(_tree.FindBranch("clusters")->GetAddress()) : 0);
  PFParticleCollection** pfParticlesP(_tree.GetBranchStatus("pfParticles") ? reinterpret_cast<PFParticleCollection**>(_tree.FindBranch("pfParticles")->GetAddress()) : 0);
  PUSummaryInfoCollection** puP(_tree.GetBranchStatus("pu") ? reinterpret_cast<PUSummaryInfoCollection**>(_tree.FindBranch("pu")->GetAddress()) : 0);
  ParticleCollection** genParticlesP(_tree.GetBranchStatus("genParticles") ? reinterpret_cast<ParticleCollection**>(_tree.FindBranch("genParticles")->GetAddress()) : 0);

  TObjArray* branches(_tree.GetListOfBranches());
  std::vector<MET**> metPs;
  std::vector<MuonCollection**> muonPs;
  std::vector<ElectronCollection**> electronPs;
  std::vector<PhotonCollection**> photonPs;
  std::vector<CaloJetCollection**> caloJetPs;
  std::vector<PFJetCollection**> pfJetPs;
  std::vector<JPTJetCollection**> jptJetPs;

  for(int iB(0); iB != branches->GetEntries(); ++iB){
    TBranch* branch(static_cast<TBranch*>(branches->At(iB)));
    TString bName(branch->GetName());
    if(!_tree.GetBranchStatus(bName)) continue;

    TString collectionName(bName(bName.First('_') + 1, bName.Length()));

    if(bName.Index("met_") == 0 && metMap.find(collectionName(0, collectionName.Length() - 1)) != metMap.end())
      metPs.push_back(reinterpret_cast<MET**>(branch->GetAddress()));
    else if(bName.Index("muons_") == 0 && muons.find(collectionName) != muons.end())
      muonPs.push_back(reinterpret_cast<MuonCollection**>(branch->GetAddress()));
    else if(bName.Index("electrons_") == 0 && electrons.find(collectionName) != electrons.end())
      electronPs.push_back(reinterpret_cast<ElectronCollection**>(branch->GetAddress()));
    else if(bName.Index("photons_") == 0 && photons.find(collectionName) != photons.end())
      photonPs.push_back(reinterpret_cast<PhotonCollection**>(branch->GetAddress()));
    else if(bName.Index("caloJets_") == 0 && caloJets.find(collectionName) != caloJets.end())
      caloJetPs.push_back(reinterpret_cast<CaloJetCollection**>(branch->GetAddress()));
    else if(bName.Index("pfJets_") == 0 && pfJets.find(collectionName) != pfJets.end())
      pfJetPs.push_back(reinterpret_cast<PFJetCollection**>(branch->GetAddress()));
    else if(bName.Index("jptJets_") == 0 && jptJets.find(collectionName) != jptJets.end())
      jptJetPs.push_back(reinterpret_cast<JPTJetCollection**>(branch->GetAddress()));
    else
      continue;
  }
    
  _tree.ResetBranchAddresses();
    
  delete beamSpotP;
  delete verticesP;
  delete tracksP;
  delete superClustersP;
  delete clustersP;
  delete pfParticlesP;
  delete puP;
  delete genParticlesP;
  for(unsigned i(0); i != metPs.size(); ++i) delete metPs[i];
  for(unsigned i(0); i != muonPs.size(); ++i) delete muonPs[i];
  for(unsigned i(0); i != electronPs.size(); ++i) delete electronPs[i];
  for(unsigned i(0); i != photonPs.size(); ++i) delete photonPs[i];
  for(unsigned i(0); i != caloJetPs.size(); ++i) delete caloJetPs[i];
  for(unsigned i(0); i != pfJetPs.size(); ++i) delete pfJetPs[i];
  for(unsigned i(0); i != jptJetPs.size(); ++i) delete jptJetPs[i];

  l1Map.releaseTree(_tree, kTRUE);
  hltMap.releaseTree(_tree, kTRUE);

  if(iTree == -1) inputTree_ = 0;
  else outputTrees_.erase(outputTrees_.begin() + iTree);
}

void
susy::Event::releaseTrees()
{
  if(inputTree_) releaseTree(*inputTree_);
  while(outputTrees_.size() != 0)
    releaseTree(*outputTrees_[0]); // releaseTree removes the tree from the vector
}

void
susy::Event::copyEvent(Event const& _orig)
{
  isRealData             = _orig.isRealData;
  runNumber              = _orig.runNumber;
  eventNumber            = _orig.eventNumber;
  luminosityBlockNumber  = _orig.luminosityBlockNumber;
  bunchCrossing          = _orig.bunchCrossing;
  avgInsRecLumi          = _orig.avgInsRecLumi;
  intgRecLumi            = _orig.intgRecLumi;
  cosmicFlag             = _orig.cosmicFlag;
  rho                    = _orig.rho;
  rhoBarrel              = _orig.rhoBarrel;
  rho25                  = _orig.rho25;
  metFilterBit           = _orig.metFilterBit;

  beamSpot               = _orig.beamSpot;

  l1Map.copy(_orig.l1Map);
  hltMap.copy(_orig.hltMap);

  vertices               = _orig.vertices;
  tracks                 = _orig.tracks;
  superClusters          = _orig.superClusters;
  clusters               = _orig.clusters;
  pfParticles            = _orig.pfParticles;
  pu                     = _orig.pu;
  genParticles           = _orig.genParticles;

  for(std::map<TString, MET>::iterator itr(metMap.begin()); itr != metMap.end(); ++itr){
    std::map<TString, MET>::const_iterator oItr(_orig.metMap.find(itr->first));
    if(oItr != _orig.metMap.end()) itr->second = oItr->second;
    else itr->second.Init();
  }
  for(std::map<TString, MuonCollection>::iterator itr(muons.begin()); itr != muons.end(); ++itr){
    std::map<TString, MuonCollection>::const_iterator oItr(_orig.muons.find(itr->first));
    if(oItr != _orig.muons.end()) itr->second = oItr->second;
    else itr->second.clear();
  }
  for(std::map<TString, ElectronCollection>::iterator itr(electrons.begin()); itr != electrons.end(); ++itr){
    std::map<TString, ElectronCollection>::const_iterator oItr(_orig.electrons.find(itr->first));
    if(oItr != _orig.electrons.end()) itr->second = oItr->second;
    else itr->second.clear();
  }
  for(std::map<TString, PhotonCollection>::iterator itr(photons.begin()); itr != photons.end(); ++itr){
    std::map<TString, PhotonCollection>::const_iterator oItr(_orig.photons.find(itr->first));
    if(oItr != _orig.photons.end()) itr->second = oItr->second;
    else itr->second.clear();
  }
  for(std::map<TString, CaloJetCollection>::iterator itr(caloJets.begin()); itr != caloJets.end(); ++itr){
    std::map<TString, CaloJetCollection>::const_iterator oItr(_orig.caloJets.find(itr->first));
    if(oItr != _orig.caloJets.end()) itr->second = oItr->second;
    else itr->second.clear();
  }
  for(std::map<TString, PFJetCollection>::iterator itr(pfJets.begin()); itr != pfJets.end(); ++itr){
    std::map<TString, PFJetCollection>::const_iterator oItr(_orig.pfJets.find(itr->first));
    if(oItr != _orig.pfJets.end()) itr->second = oItr->second;
    else itr->second.clear();
  }
  for(std::map<TString, JPTJetCollection>::iterator itr(jptJets.begin()); itr != jptJets.end(); ++itr){
    std::map<TString, JPTJetCollection>::const_iterator oItr(_orig.jptJets.find(itr->first));
    if(oItr != _orig.jptJets.end()) itr->second = oItr->second;
    else itr->second.clear();
  }
  for(std::map<TString, Float_t>::iterator itr(gridParams.begin()); itr != gridParams.end(); ++itr){
    std::map<TString, Float_t>::const_iterator oItr(_orig.gridParams.find(itr->first));
    if(oItr != _orig.gridParams.end()) itr->second = oItr->second;
    else itr->second = 0.;
  }
}


