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
// $Id: SusyEvent.cc,v 1.28 2013/03/30 13:24:11 yiiyama Exp $
//

#include "SusyEvent.h"

#include <limits>

// Use for initialization of variables that can have a value 0 and are not always filled
float const BIGVALUE(std::numeric_limits<float>::max());

void susy::PUSummaryInfo::Init() {
  numInteractions = 0;
  zPositions.clear();
  sumPTLowPT.clear();
  sumPTHighPT.clear();
  numTracksLowPT.clear();
  numTracksHighPT.clear();
  instLumi.clear();
  dataMixerRun.clear();
  dataMixerEvt.clear();
  dataMixerLumiSection.clear();
  BX = 0;
  trueNumInteractions = 0;
}

void susy::Particle::Init() {
    status         = 0;
    motherIndex    = -1;
    pdgId          = 0;
    charge         = 0;
    vertex        *= 0;
    momentum      *= 0;
}

void susy::CorrMETData::Init() {
  dmEx          = 0;
  dmEy          = 0;
  dsumEt        = 0;
  dSignificance = 0;
}

void susy::MET::Init() {
  sumEt        = 0;
  significance = 0;
  mEt         *= 0;
  vertex      *= 0;
  mEtCorr.clear();
}

void susy::Vertex::Init() {
  chi2       = 0;
  ndof       = 0;
  tracksSize = 0;
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
  vertex                          *= 0;
  momentum                        *= 0;
  extrapolatedPositions.clear();
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
  vertex                         *= 0;
  momentum                       *= 0;
  MVAregEnergyAndErr.first        = 0;
  MVAregEnergyAndErr.second       = 0;
  MVAcorrMomentum                *= 0;
  idPairs.clear();

}

void susy::Electron::Init() {
  fidBit                           = 0;
  boolPack                         = 0;
  scPixCharge                      = 0;

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

  trackPositions.clear();
  trackMomentums.clear();

  vertex     *= 0;
  momentum   *= 0;
  idPairs.clear();
}

void susy::Muon::Init() {
    type                    = 0;
    bestTrackType           = 0;
    nMatches                = 0;
    nValidHits              = 0;
    nValidTrackerHits       = 0;
    nValidMuonHits          = 0;
    nPixelLayersWithMeasurement = 0;
    nStripLayersWithMeasurement = 0;
    nChambers               = 0;
    nMatchedStations        = 0;
    timeNDof                = 0;
    timeDirection           = 0;
    timeAtIp                = 0;
    timeAtIpError           = 0;
    caloCompatibility       = 0;
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

    idPairs.clear();
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

  bTagDiscriminators.clear();

  vertex                   *= 0;
  momentum                 *= 0;

  tracklist.clear();

  pfParticleList.clear();

  puJetIdDiscriminants.clear();
  puJetIdFlags.clear();

  jecScaleFactors.clear();
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

  vertex                   *= 0;
  momentum                 *= 0;

  jecScaleFactors.clear();
}


void susy::PFParticle::Init() {
  pdgId                       = 0;
  charge                      = 0;
  ecalEnergy                  = 0;
  rawEcalEnergy               = 0;
  hcalEnergy                  = 0;
  rawHcalEnergy               = 0;
  pS1Energy                   = 0;
  pS2Energy                   = 0;

  vertex                     *= 0;
  positionAtECALEntrance     *= 0;
  momentum                   *= 0;
}



void susy::Event::Init() {
    isRealData                  = 0;
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
    metFilterBit_2              = 0;

    beamSpot                   *= 0;

    l1Map.clear();
    hltMap.clear();
    metMap.clear();

    vertices.clear();
    tracks.clear();
    superClusters.clear();
    clusters.clear();
    muons.clear();
    photons.clear();
    electrons.clear();
    caloJets.clear();
    pfJets.clear();
    jptJets.clear();
    pfParticles.clear();

    generalTracks.clear();

    pu.clear();
    simVertices.clear();
    genParticles.clear();
    gridParams.clear();
}

