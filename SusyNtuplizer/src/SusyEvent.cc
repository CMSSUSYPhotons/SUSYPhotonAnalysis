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
// $Id: SusyEvent.cc,v 1.30 2013/04/11 16:31:26 yiiyama Exp $
//

#include "SusyEvent.h"

#include "TTree.h"

#include <limits>
#include <stdexcept>
#include <algorithm>

#include <iostream>

// Use for initialization of variables that can have a value 0 and are not always filled
float const BIGVALUE(std::numeric_limits<float>::max());

void susy::PUSummaryInfo::Init() {
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

void susy::Particle::Init() {
  status         = 0;
  charge         = 0;
  motherIndex    = -1;
  pdgId          = 0;
  vertex        *= 0;
  momentum      *= 0;

  mother         = 0;
}

void
susy::Particle::fillRefs(Event const* _evt)
{
  if(motherIndex != -1)
    mother = &_evt->genParticles[motherIndex];
}

void susy::PFParticle::Init() {
  pdgId                       = 0;
  charge                      = 0;
  isPU                        = kFALSE;
  ecalEnergy                  = 0;
  hcalEnergy                  = 0;

  vertex                     *= 0;
  momentum                   *= 0;
}

void susy::MET::Init() {
  sumEt        = 0;
  significance = 0;
  mEt         *= 0;
}

void susy::Vertex::Init() {
  tracksSize = 0;
  sumPt2     = 0;
  chi2       = 0;
  ndof       = 0;
  position  *= 0;
}

void susy::Cluster::Init() {
  nCrystals = 0;
  energy    = 0;
  position *= 0;
}

void susy::SuperCluster::Init() {
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
susy::SuperCluster::fillRefs(Event const* _evt)
{
  if(seedClusterIndex != -1)
    seedCluster = &_evt->clusters[seedClusterIndex];

  basicClusters.assign(basicClusterIndices.size(), 0);
  for(unsigned iC(0); iC != basicClusterIndices.size(); ++iC)
    basicClusters[iC] = &_evt->clusters[basicClusterIndices[iC]];
}

void susy::Track::Init() {
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
susy::Track::fillRefs(Event const* _evt)
{
  if(vertexIndex != -1)
    assignedVertex = &_evt->vertices[vertexIndex];
}

void susy::Photon::Init() {
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
}

void
susy::Photon::fillRefs(Event const* _evt)
{
  if(superClusterIndex != -1)
    superCluster = &_evt->superClusters[superClusterIndex];
}

void susy::Electron::Init() {
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
susy::Electron::fillRefs(Event const* _evt)
{
  if(gsfTrackIndex != -1)
    gsfTrack = &_evt->tracks[gsfTrackIndex];
  if(closestCtfTrackIndex != -1)
    closestCtfTrack = &_evt->tracks[closestCtfTrackIndex];
  if(electronClusterIndex != -1)
    electronCluster = &_evt->clusters[electronClusterIndex];
  if(superClusterIndex != -1)
    superCluster = &_evt->superClusters[superClusterIndex];
}

void susy::Muon::Init() {
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
susy::Muon::fillRefs(Event const* _evt)
{
  if(trackIndex != -1)
    innerTrack = &_evt->tracks[trackIndex];
  if(standAloneTrackIndex != -1)
    outerTrack = &_evt->tracks[standAloneTrackIndex];
  if(combinedTrackIndex != -1)
    globalTrack = &_evt->tracks[combinedTrackIndex];
  if(tpfmsTrackIndex != -1)
    tpfmsTrack = &_evt->tracks[tpfmsTrackIndex];
  if(pickyTrackIndex != -1)
    pickyTrack = &_evt->tracks[pickyTrackIndex];
  if(dytTrackIndex != -1)
    dytTrack = &_evt->tracks[dytTrackIndex];
  if(bestTrackIndex() != -1)
    bestTrack = &_evt->tracks[bestTrackIndex()];
  if(highPtBestTrackIndex() != -1)
    highPtBestTrack = &_evt->tracks[highPtBestTrackIndex()];
}

void susy::CaloJet::Init() {
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
susy::CaloJet::fillRefs(Event const*)
{
}

void susy::PFJet::Init() {
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

  momentum                 *= 0;

  jecScaleFactors.clear();
  jecUncertainty            = 0;

  tracks.clear();
  pfParticles.clear();
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

void susy::JPTJet::Init() {
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
susy::JPTJet::fillRefs(Event const*)
{
}

susy::TriggerMap::const_iterator::const_iterator() :
  maxPos_(0),
  paths_(0),
  prescales_(0),
  decisions_(0),
  pos_(0),
  pair_("", std::pair<UInt_t, Bool_t>(0, kFALSE))
{
}

susy::TriggerMap::const_iterator::const_iterator(TriggerMap const& _map, unsigned _pos) :
  maxPos_(_map.cachedSize_),
  paths_(_map.cachedPaths_),
  prescales_(_map.cachedPrescales_),
  decisions_(_map.cachedDecisions_),
  pos_(_pos)
{
  setPair();
}

susy::TriggerMap::const_iterator&
susy::TriggerMap::const_iterator::operator++()
{
  ++pos_;
  setPair();
  return *this;
}

susy::TriggerMap::const_iterator
susy::TriggerMap::const_iterator::operator++(int)
{
  const_iterator itr(*this);
  pos_++;
  setPair();
  return itr;
}

bool
susy::TriggerMap::const_iterator::operator==(const_iterator const& _rhs) const
{
  return
    paths_ == _rhs.paths_ &&
    prescales_ == _rhs.prescales_ &&
    decisions_ == _rhs.decisions_ &&
    pos_ == _rhs.pos_;
}

void
susy::TriggerMap::const_iterator::setPair()
{
  if(pos_ < maxPos_){
    pair_.first = paths_[pos_];
    pair_.second.first = prescales_[pos_];
    pair_.second.second = decisions_[pos_] == 1;
  }
  else{
    pos_ = maxPos_;
    pair_.first = "";
    pair_.second.first = 0;
    pair_.second.second = kFALSE;
  }
}

susy::TriggerMap::TriggerMap() :
  currentConfig(),
  size_(),
  paths_(),
  prescales_(),
  decisions_(),
  cachedConfigName_(""),
  cachedSize_(0),
  cachedPaths_(0),
  cachedPrescales_(0),
  cachedDecisions_(0)
{
}

void
susy::TriggerMap::Init()
{
  for(std::map<TString, size_t>::iterator sItr(size_.begin()); sItr != size_.end(); ++sItr){
    size_t size(sItr->second);
    UInt_t* prescales(prescales_.find(sItr->first)->second);
    UChar_t* decisions(decisions_.find(sItr->first)->second);

    std::fill_n(prescales, size, 0);
    std::fill_n(decisions, size, 0);
  }
}

std::pair<UInt_t, Bool_t>
susy::TriggerMap::operator[](TString const& _path) const
{
  std::pair<UInt_t, Bool_t> value(0, kFALSE);

  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));
  if(p == cachedPaths_ + cachedSize_) return value;

  TString const& found(*p);

  if((_path.First('*') == _path.Length() - 1 && found.Index(_path.SubString(0, _path.Length() - 1)) == 0) ||
     _path == found){

    unsigned iP(p - cachedPaths_);
    value.first = cachedPrescales_[iP];
    value.second = cachedDecisions_[iP] == 1;
  }
  return value;
}

Bool_t
susy::TriggerMap::pass(TString const& _path) const
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));
  if(p == cachedPaths_ + cachedSize_) return kFALSE;

  TString const& found(*p);

  if((_path.First('*') == _path.Length() - 1 && found.Index(_path.SubString(0, _path.Length() - 1)) == 0) ||
     _path == found)
    return cachedDecisions_[p - cachedPaths_] == 1;
  else
    return kFALSE;
}

UInt_t
susy::TriggerMap::prescale(TString const& _path) const
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));
  if(p == cachedPaths_ + cachedSize_) return 0;

  TString const& found(*p);

  if((_path.First('*') == _path.Length() - 1 && found.Index(_path.SubString(0, _path.Length() - 1)) == 0) ||
     _path == found)
    return cachedPrescales_[p - cachedPaths_];
  else
    return 0;
}

Bool_t
susy::TriggerMap::menuExists(TString const& _menuName) const
{
  return size_.find(_menuName) != size_.end();
}

void
susy::TriggerMap::addMenu(TString const& _menuName, std::vector<TString> const& _sortedFilters)
{
  if(menuExists(_menuName))
    throw std::runtime_error((_menuName + " already exists in TriggerMap").Data());

  size_t size(_sortedFilters.size());
  TString* paths(new TString[size]);
  UInt_t* prescales(new UInt_t[size]);
  UChar_t* decisions(new UChar_t[size]);

  std::copy(_sortedFilters.begin(), _sortedFilters.end(), paths);
  std::fill_n(prescales, size, 0);
  std::fill_n(decisions, size, 0);  

  size_.insert(std::pair<TString, size_t>(_menuName, size));
  paths_.insert(std::pair<TString, TString*>(_menuName, paths));
  prescales_.insert(std::pair<TString, UInt_t*>(_menuName, prescales));
  decisions_.insert(std::pair<TString, UChar_t*>(_menuName, decisions));
}

void
susy::TriggerMap::set(TString const& _path, UInt_t _prescale, Bool_t _decision)
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));
  if(p == cachedPaths_ + cachedSize_ || *p != _path) return;

  UInt_t offset(p - cachedPaths_);

  cachedPrescales_[offset] = _prescale;
  cachedDecisions_[offset] = _decision ? 1 : 0;
}

void
susy::TriggerMap::clear()
{
  currentConfig = "";
  for(std::map<TString, size_t>::iterator itr(size_.begin()); itr != size_.end(); ++itr){
    delete [] paths_[itr->first];
    delete [] prescales_[itr->first];
    delete [] decisions_[itr->first];
  }  
  size_.clear();
  paths_.clear();
  prescales_.clear();
  decisions_.clear();

  cachedConfigName_ = "";
  cachedSize_ = 0;
  cachedPaths_ = 0;
  cachedPrescales_ = 0;
  cachedDecisions_ = 0;
}

size_t
susy::TriggerMap::size() const
{
  setCache_();
  return cachedSize_;
}

UInt_t*
susy::TriggerMap::prescales(TString const& _config) const
{
  std::map<TString, UInt_t*>::const_iterator itr(prescales_.find(_config));
  if(itr != prescales_.end()) return itr->second;
  else return 0;
}

UChar_t*
susy::TriggerMap::decisions(TString const& _config) const
{
  std::map<TString, UChar_t*>::const_iterator itr(decisions_.find(_config));
  if(itr != decisions_.end()) return itr->second;
  else return 0;
}

susy::TriggerMap::const_iterator
susy::TriggerMap::begin() const
{
  setCache_();

  return const_iterator(*this, 0);
}

susy::TriggerMap::iterator
susy::TriggerMap::begin()
{
  setCache_();

  return iterator(*this, 0);
}

susy::TriggerMap::const_iterator
susy::TriggerMap::end() const
{
  setCache_();

  return const_iterator(*this, cachedSize_);
}

susy::TriggerMap::iterator
susy::TriggerMap::end()
{
  setCache_();

  return iterator(*this, cachedSize_);
}

susy::TriggerMap::const_iterator
susy::TriggerMap::find(TString const& _path) const
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));
  if(p == cachedPaths_ + cachedSize_ || *p != _path) return end();

  return const_iterator(*this, p - cachedPaths_);
}

susy::TriggerMap::iterator
susy::TriggerMap::find(TString const& _path)
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));
  if(p == cachedPaths_ + cachedSize_ || *p != _path) return end();

  return iterator(*this, p - cachedPaths_);
}

susy::TriggerMap::const_iterator
susy::TriggerMap::lower_bound(TString const& _path) const
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));

  return const_iterator(*this, p - cachedPaths_);
}

susy::TriggerMap::iterator
susy::TriggerMap::lower_bound(TString const& _path)
{
  setCache_();

  TString const* p(std::lower_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));

  return iterator(*this, p - cachedPaths_);
}

susy::TriggerMap::const_iterator
susy::TriggerMap::upper_bound(TString const& _path) const
{
  setCache_();

  TString const* p(std::upper_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));

  return const_iterator(*this, p - cachedPaths_);
}

susy::TriggerMap::iterator
susy::TriggerMap::upper_bound(TString const& _path)
{
  setCache_();

  TString const* p(std::upper_bound(cachedPaths_, cachedPaths_ + cachedSize_, _path));

  return iterator(*this, p - cachedPaths_);
}

void
susy::TriggerMap::setCache_() const
{
  if(currentConfig == cachedConfigName_) return;
  if(!menuExists(currentConfig))
    throw std::runtime_error((currentConfig + " does not exist in TriggerMap").Data());

  cachedConfigName_ = currentConfig;

  cachedSize_ = size_.find(currentConfig)->second;
  cachedPaths_ = paths_.find(currentConfig)->second;
  cachedPrescales_ = prescales_.find(currentConfig)->second;
  cachedDecisions_ = decisions_.find(currentConfig)->second;
}

susy::Event::Event()
{
  Init();
}

susy::Event::~Event()
{
  Init();

  while(trees_.size() != 0)
    releaseTree(*trees_[0].first);
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
  metFilterBit                = 0;

  beamSpot                   *= 0;

  l1Map.Init();
  hltMap.Init();

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
susy::Event::bindTree(TTree& _tree, Bool_t _readMode/* = kTRUE*/)
{
  if(_readMode){
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

    _tree.SetBranchAddress("l1Config", new TString*(&l1Map.currentConfig));
    _tree.SetBranchAddress("hltConfig", new TString*(&hltMap.currentConfig));
    _tree.SetBranchAddress("beamSpot", new TVector3*(&beamSpot));
    _tree.SetBranchAddress("vertices", new VertexCollection*(&vertices));
    _tree.SetBranchAddress("tracks", new TrackCollection*(&tracks));
    _tree.SetBranchAddress("superClusters", new SuperClusterCollection*(&superClusters));
    _tree.SetBranchAddress("clusters", new ClusterCollection*(&clusters));
    _tree.SetBranchAddress("pfParticles", new PFParticleCollection*(&pfParticles));
    _tree.SetBranchAddress("pu", new PUSummaryInfoCollection*(&pu));
    _tree.SetBranchAddress("genParticles", new ParticleCollection*(&genParticles));

    TObjArray* branches(_tree.GetListOfBranches());

    for(int iB(0); iB != branches->GetEntries(); ++iB){
      TBranch* branch(static_cast<TBranch*>(branches->At(iB)));
      TString bName(branch->GetName());
      TString collectionName(bName(bName.First('_') + 1, bName.Length()));
      void** add(0);

      if(bName.Index("met_") == 0)
        add = reinterpret_cast<void**>(new MET*(&metMap[collectionName]));
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

    for(int iB(0); iB != branches->GetEntries(); ++iB){
      TBranch* branch(static_cast<TBranch*>(branches->At(iB)));
      TString bName(branch->GetName());
      TriggerMap* map(0);

      if(bName.Index("l1Bits_") == 0) map = &l1Map;
      else if(bName.Index("hltBits_") == 0) map = &hltMap;
      else continue;

      TObjArray* leaves(branch->GetListOfLeaves());

      std::vector<TString> paths;
      for(int iL(0); iL != leaves->GetEntries(); ++iL)
        paths.push_back(leaves->At(iL)->GetName());

      TString config(bName(bName.First('_') + 1, bName.Length()));

      map->addMenu(config, paths);
      _tree.SetBranchAddress(bName, map->decisions(config));
      _tree.SetBranchAddress(bName.ReplaceAll("Bits", "Prescales"), map->prescales(config));
    }
  }
  else{
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

    _tree.Branch("l1Config", "TString", new TString*(&l1Map.currentConfig));
    _tree.Branch("hltConfig", "TString", new TString*(&hltMap.currentConfig));
    _tree.Branch("beamSpot", "TVector3", new TVector3*(&beamSpot));
    _tree.Branch("vertices", "std::vector<susy::Vertex>", new VertexCollection*(&vertices));
    _tree.Branch("tracks", "std::vector<susy::Track>", new TrackCollection*(&tracks));
    _tree.Branch("superClusters", "std::vector<susy::SuperCluster>", new SuperClusterCollection*(&superClusters));
    _tree.Branch("clusters", "std::vector<susy::Cluster>", new ClusterCollection*(&clusters));
    _tree.Branch("pfParticles", "std::vector<susy::PFParticle>", new PFParticleCollection*(&pfParticles));
    _tree.Branch("pu", "std::vector<susy::PUSummaryInfo>", new PUSummaryInfoCollection*(&pu));
    _tree.Branch("genParticles", "std::vector<susy::Particle>", new ParticleCollection*(&genParticles));

    for(std::map<TString, MET>::iterator itr(metMap.begin()); itr != metMap.end(); ++itr)
      _tree.Branch("met_" + itr->first, "susy::MET", new MET*(&itr->second));
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
  }

  trees_.push_back(std::pair<TTree*, Bool_t>(&_tree, _readMode));
}

void
susy::Event::releaseTree(TTree& _tree)
{
  for(std::vector<std::pair<TTree*, Bool_t> >::iterator tItr(trees_.begin()); tItr != trees_.end(); ++tItr){
    if(tItr->first == &_tree){
      TString** l1ConfigP(reinterpret_cast<TString**>(_tree.FindBranch("l1Config")->GetAddress()));
      TString** hltConfigP(reinterpret_cast<TString**>(_tree.FindBranch("hltConfig")->GetAddress()));
      TVector3** beamSpotP(reinterpret_cast<TVector3**>(_tree.FindBranch("beamSpot")->GetAddress()));
      VertexCollection** verticesP(reinterpret_cast<VertexCollection**>(_tree.FindBranch("vertices")->GetAddress()));
      TrackCollection** tracksP(reinterpret_cast<TrackCollection**>(_tree.FindBranch("tracks")->GetAddress()));
      SuperClusterCollection** superClustersP(reinterpret_cast<SuperClusterCollection**>(_tree.FindBranch("superClusters")->GetAddress()));
      ClusterCollection** clustersP(reinterpret_cast<ClusterCollection**>(_tree.FindBranch("clusters")->GetAddress()));
      PFParticleCollection** pfParticlesP(reinterpret_cast<PFParticleCollection**>(_tree.FindBranch("pfParticles")->GetAddress()));
      PUSummaryInfoCollection** puP(reinterpret_cast<PUSummaryInfoCollection**>(_tree.FindBranch("pu")->GetAddress()));
      ParticleCollection** genParticlesP(reinterpret_cast<ParticleCollection**>(_tree.FindBranch("genParticles")->GetAddress()));

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

        if(bName.Index("met_") == 0)
          metPs.push_back(reinterpret_cast<MET**>(branch->GetAddress()));
        else if(bName.Index("muons_") == 0)
          muonPs.push_back(reinterpret_cast<MuonCollection**>(branch->GetAddress()));
        else if(bName.Index("electrons_") == 0)
          electronPs.push_back(reinterpret_cast<ElectronCollection**>(branch->GetAddress()));
        else if(bName.Index("photons_") == 0)
          photonPs.push_back(reinterpret_cast<PhotonCollection**>(branch->GetAddress()));
        else if(bName.Index("caloJets_") == 0)
          caloJetPs.push_back(reinterpret_cast<CaloJetCollection**>(branch->GetAddress()));
        else if(bName.Index("pfJets_") == 0)
          pfJetPs.push_back(reinterpret_cast<PFJetCollection**>(branch->GetAddress()));
        else if(bName.Index("jptJets_") == 0)
          jptJetPs.push_back(reinterpret_cast<JPTJetCollection**>(branch->GetAddress()));
        else
          continue;
      }
    
      _tree.ResetBranchAddresses();
    
      delete l1ConfigP;
      delete hltConfigP;
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

      std::vector<std::pair<TTree*, Bool_t> >::iterator cur(tItr);
      --cur;
      trees_.erase(tItr);
      tItr = cur;
    }
  }
}
  
void
susy::Event::setTriggerTable(TString const& _menuName, std::vector<std::string> const& _filters, Bool_t _isHLT/* = kTRUE*/)
{
  UInt_t* prescales(0);
  UChar_t* decisions(0);
  TString trigType(_isHLT ? "hlt" : "l1");

  std::vector<TString> filters(_filters.size());
  std::copy(_filters.begin(), _filters.end(), filters.begin());
  std::sort(filters.begin(), filters.end());

  if(_isHLT){
    hltMap.addMenu(_menuName, filters);
    prescales = hltMap.prescales(_menuName);
    decisions = hltMap.decisions(_menuName);
  }
  else{
    l1Map.addMenu(_menuName, filters);
    prescales = l1Map.prescales(_menuName);
    decisions = l1Map.decisions(_menuName);
  }

  TString psLeaves;
  TString bitLeaves;
  for(unsigned iF(0); iF != filters.size(); ++iF){
    psLeaves += filters[iF] + "/i:";
    bitLeaves += filters[iF] + "/b:";
  }
  psLeaves.Remove(psLeaves.Length() - 1);
  bitLeaves.Remove(bitLeaves.Length() - 1);

  for(unsigned iT(0); iT != trees_.size(); ++iT){
    if(trees_[iT].second) continue; // is a read-only tree
    TTree* tree(trees_[iT].first);

    TBranch* psBranch(tree->Branch(trigType + "Prescales_" + _menuName, prescales, psLeaves));
    TBranch* bitBranch(tree->Branch(trigType + "Bits_" + _menuName, decisions, bitLeaves));

    psBranch->SetFirstEntry(tree->GetEntries());
    bitBranch->SetFirstEntry(tree->GetEntries());
  }
}
