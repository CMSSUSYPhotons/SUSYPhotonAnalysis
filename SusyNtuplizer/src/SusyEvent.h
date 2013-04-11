// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEvent.h
//
/*

 Description: Objects definitions used for SusyNtuples

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEvent.h,v 1.39 2013/04/03 22:04:53 yiiyama Exp $
//

#ifndef SusyEvent_h
#define SusyEvent_h

#include <vector>
#include <string>
#include <map>
#include <math.h>

#include <TLorentzVector.h>

class TTree;

namespace susy {

  const float etaGapBegin = 1.4442;
  const float etaGapEnd = 1.566;
  const float etaGap = 1.499;
  const float etaMax = 2.5;

  enum MetFilters {
    kCSCBeamHalo,
    kHcalNoise,
    kEcalDeadCellTP,
    kEcalDeadCellBE,
    kHcalLaser,
    kTrackingFailure,
    kEEBadSC,
    kHcalLaser2012,
    kEcalLaserCorr,
    kManyStripClus53X,
    kTooManyStripClus53X,
    kLogErrorTooManyClusters,
    kEERingOfFire,
    kInconsistentMuon,
    kGreedyMuon,
    nMetFilters
  };

  // b-tagging vector positions
  enum BTagDiscriminators {
    kTCHE,   // TrackCountingHighEff
    kTCHP,   // TrackCountingHighPur
    kJP,     // JetProbability
    kJBP,    // JetBProbability
    kSSV,    // SimpleSecondaryVertex
    kCSV,    // CombinedSecondaryVertex
    kCSVMVA, // CombinedSecondaryVertexMVA
    kSE,     // SoftElectron
    kSM,     // SoftMuon
    nBTagDiscriminators
  };

  enum PUJetIdAlgorithms {
    kPUJetIdFull,
    kPUJetIdCutBased,
    kPUJetIdSimple,
    nPUJetIdAlgorithms
  };

  class Event; // forward declaration for fillRefs methods

  class PUSummaryInfo { /*each PUSummaryInfo object holds information for one BX (early, in time,
                          or late)*/

  public:

    PUSummaryInfo()  { Init(); }
    ~PUSummaryInfo() { Init(); }
    void Init();

    //all info below from https://twiki.cern.ch/twiki/bin/view/CMS/PileupInformation
    /*low_cut = 0.1 GeV, high_cut = 0.5 GeV, tracks summed/counted are TrackingParticles from
      simulation truth info*/

    Char_t BX; //to which bunch crossing does this interaction belong?
    UChar_t numInteractions; //the number of pileup interactions that have been added to the event
    Float_t trueNumInteractions;

    std::vector<Float_t> zPositions; /*the true primary vertex position along the z axis for each
                                     added interaction*/
    std::vector<Float_t> sumPTLowPT; /*the sum of the transverse momentum of the tracks originating
                                     from each interaction, where track pT > low_cut*/
    std::vector<Float_t> sumPTHighPT; /*the sum of the transverse momentum of the tracks originating
                                      from each interaction, where track pT > high_cut*/
    std::vector<UShort_t> numTracksLowPT; /*the number of tracks originating from each interaction,
                                       where track pT > low_cut*/
    std::vector<UShort_t> numTracksHighPT; /*the number of tracks originating from each interaction,
                                        where track pT > high_cut*/
    std::vector<Float_t> instLumi; //for PU from DataMixer

    std::vector<UInt_t> dataMixerRun;
    std::vector<UInt_t> dataMixerEvt;
    std::vector<UInt_t> dataMixerLumiSection;

  };


  class Particle {

  public:

    Particle()  { Init(); }
    ~Particle() { Init(); }
    void Init();
    void fillRefs(Event const*);

    UChar_t status;
    Char_t  charge;
    Short_t motherIndex;
    Int_t   pdgId;

    TVector3 vertex;
    TLorentzVector momentum;

    Particle const* mother;  //!

  };

  class PFParticle {

  public:

    PFParticle() { Init(); }
    ~PFParticle() { Init(); }
    void Init();

    Char_t  charge;
    Bool_t  isPU;
    Short_t pdgId;

    Float_t ecalEnergy;
    Float_t hcalEnergy;

    TVector3 vertex;
    TLorentzVector momentum;

  };

  class MET {

  public:

    MET()  { Init(); }
    ~MET() { Init(); }
    void Init();

    Float_t met() const {  return mEt.Mod(); }
    Float_t metX() const { return mEt.X(); }
    Float_t metY() const { return mEt.Y(); }

    Float_t  sumEt;
    Float_t  significance;
    TVector2 mEt;

  };



  class Vertex {

  public:

    Vertex()  { Init(); }
    ~Vertex() { Init(); }
    void Init();
    Bool_t isFake() const { return (chi2 == 0 && ndof == 0 && tracksSize == 0); }

    UShort_t tracksSize;
    Float_t  sumPt2;
    Float_t  chi2;
    Float_t  ndof;
    TVector3 position;

  };


  class Cluster {

  public:

    Cluster()  { Init(); }
    ~Cluster() { Init(); }
    void Init();

    UChar_t  nCrystals;
    Float_t  energy;
    TVector3 position;

  };


  class SuperCluster {

  public:

    SuperCluster()  { Init(); }
    ~SuperCluster() { Init(); }
    void Init();
    void fillRefs(Event const*);

    Short_t  seedClusterIndex;                       // index in vector<Cluster> below
    Float_t  energy;
    Float_t  preshowerEnergy;
    Float_t  phiWidth;
    Float_t  etaWidth;
    TVector3 position;
    std::vector<UShort_t> basicClusterIndices;

    const Cluster* seedCluster;                 //!
    std::vector<const Cluster*> basicClusters;  //!
  };


  class Track {

  public:

    Track()  { Init(); }
    ~Track() { Init(); }
    void Init();
    void fillRefs(Event const*);

    // derived quantities
    Float_t normChi2() const { return (ndof != 0) ? chi2/ndof : chi2*1e6; }
    Float_t qoverp() const { return charge/momentum.P(); }
    Float_t lambda() const { return M_PI/2 - momentum.Theta(); }
    Float_t dsz() const { return vertex.Z()*momentum.Pt()/momentum.P() - (vertex.X()*momentum.Px()+vertex.Y()*momentum.Py())/momentum.Pt() * momentum.Pz()/momentum.P(); }
    Float_t dz() const { return vertex.Z() - (vertex.X()*momentum.Px()+vertex.Y()*momentum.Py())/momentum.Pt() * (momentum.Pz()/momentum.Pt()); }
    Float_t dxy() const { return (-vertex.X()*momentum.Py() + vertex.Y()*momentum.Px())/momentum.Pt(); }
    Float_t d0() const { return -dxy(); }
    Float_t phi() const { return momentum.Phi(); }
    Bool_t loose() const {         return (quality & ( 0x1 << 0)); }
    Bool_t tight() const {         return (quality & ( 0x1 << 1)); }
    Bool_t highPurity() const {    return (quality & ( 0x1 << 2)); }
    Bool_t confirmed() const {     return (quality & ( 0x1 << 3)); }
    Bool_t goodIterative() const { return (confirmed() || highPurity()); }

    UChar_t        algorithm;
    UChar_t        quality;
    UChar_t        numberOfValidHits;
    UChar_t        numberOfValidTrackerHits;
    UChar_t        numberOfValidMuonHits;
    UChar_t        numberOfValidPixelHits;
    UChar_t        numberOfValidStripHits;

    Short_t        vertexIndex;

    Float_t        chi2;
    Float_t        ndof;
    Float_t        charge;
    Float_t        error[5]; // (qoverp, lambda, phi, dxy, dsz) for general tracks, (qoverp, lambda, phi) for GSF tracks.
    Float_t        ptError;
    TVector3       vertex;   // position of the point of reference for momentum calculation (not the assigned PV)
    TLorentzVector momentum;

    const Vertex*  assignedVertex; //!
  };


  class Photon {

  public:

    Photon()  { Init(); }
    ~Photon() { Init(); }
    void Init();
    void fillRefs(Event const*);

    // fiducial bits
    Bool_t isEB() const        { return (fidBit & (0x1 << 0)); }
    Bool_t isEE() const        { return (fidBit & (0x1 << 1)); }
    Bool_t isEBEtaGap() const  { return (fidBit & (0x1 << 2)); }
    Bool_t isEBPhiGap() const  { return (fidBit & (0x1 << 3)); }
    Bool_t isEERingGap() const { return (fidBit & (0x1 << 4)); }
    Bool_t isEEDeeGap() const  { return (fidBit & (0x1 << 5)); }
    Bool_t isEBEEGap() const   { return (fidBit & (0x1 << 6)); }
    Bool_t isPF() const        { return (fidBit & (0x1 << 7)); }

    Float_t hcalTowerSumEtConeDR04() const { return (hcalDepth1TowerSumEtConeDR04+hcalDepth2TowerSumEtConeDR04); }
    Float_t hcalTowerSumEtConeDR03() const { return (hcalDepth1TowerSumEtConeDR03+hcalDepth2TowerSumEtConeDR03); }
    Float_t r1x5() const { return ((e5x5 > 0) ? e1x5/e5x5 : 0); }
    Float_t r2x5() const { return ((e5x5 > 0) ? e2x5/e5x5 : 0); }

    Int_t          fidBit;
    Int_t          nPixelSeeds;
    Bool_t         passelectronveto;
    Float_t        hadronicOverEm;
    Float_t        hadTowOverEm; //2012 hOverE
    Float_t        hadronicDepth1OverEm;
    Float_t        hadronicDepth2OverEm;
    Float_t        e1x2;
    Float_t        e1x5;
    Float_t        e2x5;
    Float_t        e3x3;
    Float_t        e5x5;
    Float_t        maxEnergyXtal;
    Float_t        sigmaEtaEta;
    Float_t        sigmaIetaIeta;
    Float_t        sigmaIphiIphi;
    Float_t        r9;

    Float_t        ecalRecHitSumEtConeDR04;
    Float_t        hcalDepth1TowerSumEtConeDR04;
    Float_t        hcalDepth2TowerSumEtConeDR04;
    Float_t        hcalIsoConeDR04_2012;
    Float_t        trkSumPtSolidConeDR04;
    Float_t        trkSumPtHollowConeDR04;
    UChar_t        nTrkSolidConeDR04;
    UChar_t        nTrkHollowConeDR04;

    Float_t        ecalRecHitSumEtConeDR03;
    Float_t        hcalDepth1TowerSumEtConeDR03;
    Float_t        hcalDepth2TowerSumEtConeDR03;
    Float_t        hcalIsoConeDR03_2012;
    Float_t        trkSumPtSolidConeDR03;
    Float_t        trkSumPtHollowConeDR03;
    UChar_t        nTrkSolidConeDR03;
    UChar_t        nTrkHollowConeDR03;

    // calculated from alternative code
    Float_t        chargedHadronIso;
    Float_t        neutralHadronIso;
    Float_t        photonIso;

    // read from IsoDeposit
    Float_t        chargedHadronIsoDeposit;
    Float_t        neutralHadronIsoDeposit;
    Float_t        photonIsoDeposit;

    Float_t        seedTime; // seed timing

    // MIP Variables

    Float_t        mipChi2;
    Float_t        mipTotEnergy;
    Float_t        mipSlope;
    Float_t        mipIntercept;
    Int_t          mipNhitCone;
    Bool_t         mipIsHalo;

    // Conversion info
    Bool_t         convInfo;
    Float_t        convDist;
    Float_t        convDcot;
    Float_t        convVtxChi2;
    Float_t        convVtxNdof;
    TVector3       convVertex;
    Float_t        convDxy;
    Float_t        convDz;
    Float_t        convLxy;
    Float_t        convLz;
    Float_t        convZofPVFromTracks;
    Int_t          convTrackChargeProd;
    Int_t          convTrack1nHit;
    Int_t          convTrack2nHit;
    Float_t        convTrack1chi2;
    Float_t        convTrack2chi2;
    Float_t        convTrack1pT;
    Float_t        convTrack2pT;
    Float_t        convTrack1InnerZ;
    Float_t        convTrack2InnerZ;
    Float_t        convTrack1InnerX;
    Float_t        convTrack2InnerX;
    Float_t        convTrack1InnerY;
    Float_t        convTrack2InnerY;
    Float_t        convTrack1Signedd0;
    Float_t        convTrack2Signedd0;

    Short_t        superClusterIndex;
    Float_t        superClusterPreshowerEnergy;
    Float_t        superClusterPhiWidth;
    Float_t        superClusterEtaWidth;
    TVector3       caloPosition;

    std::pair<Double_t, Double_t> MVAregEnergyAndErr;
    TLorentzVector MVAcorrMomentum;

    TLorentzVector momentum;

    const SuperCluster* superCluster; //!

  };


  class Electron {

  public:

    Electron()  { Init(); }
    ~Electron() { Init(); }
    void Init();
    void fillRefs(Event const*);

    // fiducial bits
    Bool_t isEB() const        { return (fidBit & (0x1 << 0)); }
    Bool_t isEE() const        { return (fidBit & (0x1 << 1)); }
    Bool_t isEBEEGap() const   { return (fidBit & (0x1 << 2)); }
    Bool_t isEBEtaGap() const  { return (fidBit & (0x1 << 3)); }
    Bool_t isEBPhiGap() const  { return (fidBit & (0x1 << 4)); }
    Bool_t isEEDeeGap() const  { return (fidBit & (0x1 << 5)); }
    Bool_t isEERingGap() const { return (fidBit & (0x1 << 6)); }
    Bool_t isEBGap() const     { return (isEBEtaGap() || isEBPhiGap()); }
    Bool_t isEEGap() const     { return (isEEDeeGap() || isEERingGap()); }
    Bool_t isGap() const       { return (isEBGap() || isEEGap() || isEBEEGap()); }

    // boolean variables packed in boolPack
    Bool_t isGsfCtfScPixChargeConsistent() const { return (boolPack & (0x1 << 0)); }
    Bool_t isGsfScPixChargeConsistent() const    { return (boolPack & (0x1 << 1)); }
    Bool_t isGsfCtfChargeConsistent() const      { return (boolPack & (0x1 << 2)); }
    Bool_t ecalDrivenSeed() const                { return (boolPack & (0x1 << 3)); }
    Bool_t trackerDrivenSeed() const             { return (boolPack & (0x1 << 4)); }
    Bool_t passingCutBasedPreselection() const   { return (boolPack & (0x1 << 5)); }
    Bool_t passingMvaPreselection() const        { return (boolPack & (0x1 << 6)); }
    Bool_t ambiguous() const                     { return (boolPack & (0x1 << 7)); }
    Bool_t isEcalEnergyCorrected() const         { return (boolPack & (0x1 << 8)); }
    Bool_t isEnergyScaleCorrected() const        { return (boolPack & (0x1 << 9)); }
    Bool_t convFlags() const                     { return (boolPack & (0x1 << 10)); }
    Bool_t isPF() const                          { return (boolPack & (0x1 << 11)); }
    Bool_t ecalDriven() const                    { return (ecalDrivenSeed() && passingCutBasedPreselection()); }

    Float_t hcalOverEcal() const { return (hcalDepth1OverEcal + hcalDepth2OverEcal); }
    Float_t dr03HcalTowerSumEt() const { return (dr03HcalDepth1TowerSumEt + dr03HcalDepth2TowerSumEt); }
    Float_t dr04HcalTowerSumEt() const { return (dr04HcalDepth1TowerSumEt + dr04HcalDepth2TowerSumEt); }

    UChar_t        fidBit;
    UShort_t       boolPack;
    Char_t         scPixCharge;

    Float_t        eSuperClusterOverP;
    Float_t        eSeedClusterOverP;
    Float_t        eSeedClusterOverPout;
    Float_t        eEleClusterOverPout;
    Float_t        deltaEtaSuperClusterTrackAtVtx;
    Float_t        deltaEtaSeedClusterTrackAtCalo;
    Float_t        deltaEtaEleClusterTrackAtCalo;
    Float_t        deltaPhiSuperClusterTrackAtVtx;
    Float_t        deltaPhiSeedClusterTrackAtCalo;
    Float_t        deltaPhiEleClusterTrackAtCalo;

    Float_t        shFracInnerHits;

    Float_t        sigmaEtaEta;
    Float_t        sigmaIetaIeta;
    Float_t        sigmaIphiIphi;
    Float_t        e1x5;
    Float_t        e2x5Max;
    Float_t        e5x5;
    Float_t        r9;
    Float_t        hcalDepth1OverEcal;          // hadronic energy on depth1 / em enrgy
    Float_t        hcalDepth2OverEcal;          // hadronic energy on depth2 / em enrgy
    Float_t        hcalOverEcalBc;   //2012 hOverE

    Float_t        dr03TkSumPt;
    Float_t        dr03EcalRecHitSumEt;
    Float_t        dr03HcalDepth1TowerSumEt;
    Float_t        dr03HcalDepth2TowerSumEt;
    Float_t        dr03HcalDepth1TowerSumEtBc;  //2012 iso with new hOverE
    Float_t        dr03HcalDepth2TowerSumEtBc;  //2012 iso with new hOverE

    Float_t        dr04TkSumPt;
    Float_t        dr04EcalRecHitSumEt;
    Float_t        dr04HcalDepth1TowerSumEt;
    Float_t        dr04HcalDepth2TowerSumEt;
    Float_t        dr04HcalDepth1TowerSumEtBc;  //2012 iso with new hOverE
    Float_t        dr04HcalDepth2TowerSumEtBc;  //2012 iso with new hOverE

    // Conversion info
    Float_t        convDist;
    Float_t        convDcot;
    Float_t        convRadius;

    Float_t        chargedHadronIso;
    Float_t        neutralHadronIso;
    Float_t        photonIso;

    // MVA ID calculated within PF sequence at reconstruction
    Int_t          mvaStatus;
    Float_t        mva;

    // MVA ID calculated using EGammaAnalysisTools
    Float_t        mvaTrig;
    Float_t        mvaNonTrig;

    Char_t         bremClass;
    Float_t        fbrem;

    Float_t        ecalEnergy;                  // corrected
    Float_t        ecalEnergyError;             // correction error
    Float_t        trackMomentumError;

    Short_t        gsfTrackIndex;
    Short_t        closestCtfTrackIndex;
    Short_t        electronClusterIndex;
    Short_t        superClusterIndex;

    Int_t          nMissingHits;
    Bool_t         passConversionVeto;


    TVector3       trackPositionAtVtx;
    TVector3       trackPositionAtCalo;
    TLorentzVector trackMomentumAtVtx;
    TLorentzVector trackMomentumAtCalo;
    TLorentzVector trackMomentumOut;
    TLorentzVector trackMomentumAtEleClus;
    TLorentzVector trackMomentumAtVtxWithConstraint;

    TVector3       vertex;
    TLorentzVector momentum;

    const Track*        gsfTrack;        //!
    const Track*        closestCtfTrack; //!
    const Cluster*      electronCluster; //!
    const SuperCluster* superCluster;    //!
  };



  class Muon {

  public:

    Muon()  { Init(); }
    ~Muon() { Init(); }
    void Init();
    void fillRefs(Event const*);

    // muon type
    Bool_t isGlobalMuon() const     { return (type & (0x1 << 1)); }
    Bool_t isTrackerMuon() const    { return (type & (0x1 << 2)); }
    Bool_t isStandAloneMuon() const { return (type & (0x1 << 3)); }
    Bool_t isCaloMuon() const       { return (type & (0x1 << 4)); }
    Bool_t isPFMuon() const         { return (type & (0x1 << 5)); }
    Bool_t tmLastStationLoose() const      { return (qualityFlags & (0x1 << 0)); }
    Bool_t tmLastStationTight() const      { return (qualityFlags & (0x1 << 1)); }
    Bool_t tmOneStationLoose()  const      { return (qualityFlags & (0x1 << 2)); }
    Bool_t tmOneStationTight()  const      { return (qualityFlags & (0x1 << 3)); }
    Bool_t tmLastStationLowPtLoose() const { return (qualityFlags & (0x1 << 4)); }
    Bool_t tmLastStationLowPtTight() const { return (qualityFlags & (0x1 << 5)); }
    Short_t bestTrackIndex() const {
      switch(bestTrackType){
      case 1: return trackIndex; case 2: return standAloneTrackIndex;
      case 3: return combinedTrackIndex; case 4: return tpfmsTrackIndex;
      case 5: return pickyTrackIndex; case 6: return dytTrackIndex;
      default: return -1;
      }
    }
    Short_t highPtBestTrackIndex() const {
      switch(highPtBestTrackType){
      case 1: return trackIndex; case 2: return standAloneTrackIndex;
      case 3: return combinedTrackIndex; case 4: return tpfmsTrackIndex;
      case 5: return pickyTrackIndex; case 6: return dytTrackIndex;
      default: return -1;
      }
    }
    UInt_t nTrackerLayersWithMeasurement() const {
      return nPixelLayersWithMeasurement + nStripLayersWithMeasurement;
    }

    UChar_t        type;
    UChar_t        bestTrackType;
    UChar_t        highPtBestTrackType;    // best high-Pt track type from muon::tevOptimized

    UChar_t        qualityFlags;           // results of various muon::isGoodMuon calls

    // Using SegmentAndTrackArbitration for all
    UChar_t        nChambers;              // number of muon chambers the track traversed through (regardless of segment existence)
    UChar_t        nMatches;               // number of muon chambers with matched segments (<= nChambers)
    UChar_t        stationMask;            // bits 0-3 -> DT stations 1-4, bits 4-7 -> CSC stations 1-4. Number of non-zero bits = nMatchedStations
    UChar_t        nMatchedStations;       // number of muon stations with matched segments (<= nMatches)
    UChar_t        nValidHits;             // values from combinedMuon (= globalTrack)
    UChar_t        nValidTrackerHits;      // *
    UChar_t        nValidMuonHits;         // *
    UChar_t        nPixelLayersWithMeasurement; // *
    UChar_t        nStripLayersWithMeasurement; // *

    UChar_t        timeNDof;         // null value implies timing measurement is invalid
    Char_t         timeDirection;
    Float_t        timeAtIp;
    Float_t        timeAtIpError;
    Float_t        caloCompatibility;
    Float_t        segmentCompatibility;
    Float_t        emEnergy;
    Float_t        hadEnergy;
    Float_t        trackIsoR03;
    Float_t        ecalIsoR03;
    Float_t        hcalIsoR03;
    Float_t        trackIsoR05;
    Float_t        ecalIsoR05;
    Float_t        hcalIsoR05;

    Float_t        sumChargedHadronPt03;
    Float_t        sumChargedParticlePt03;
    Float_t        sumNeutralHadronEt03;
    Float_t        sumPhotonEt03;
    Float_t        sumNeutralHadronEtHighThreshold03;
    Float_t        sumPhotonEtHighThreshold03;
    Float_t        sumPUPt03;

    Float_t        sumChargedHadronPt04;
    Float_t        sumChargedParticlePt04;
    Float_t        sumNeutralHadronEt04;
    Float_t        sumPhotonEt04;
    Float_t        sumNeutralHadronEtHighThreshold04;
    Float_t        sumPhotonEtHighThreshold04;
    Float_t        sumPUPt04;

    Short_t        trackIndex;             // tracker only
    Short_t        standAloneTrackIndex;   // muon detector only
    Short_t        combinedTrackIndex;     // combined
    Short_t        tpfmsTrackIndex;
    Short_t        pickyTrackIndex;
    Short_t        dytTrackIndex;
    TLorentzVector momentum;

    const Track* innerTrack;       //!
    const Track* outerTrack;       //!
    const Track* globalTrack;      //!
    const Track* tpfmsTrack;       //!
    const Track* pickyTrack;       //!
    const Track* dytTrack;         //!
    const Track* bestTrack;        //!
    const Track* highPtBestTrack;  //!
  };


  class CaloJet {

  public:

    CaloJet()  { Init(); }
    ~CaloJet() { Init(); }
    void Init();
    void fillRefs(Event const*);

    // Basic Jet Info
    Float_t        partonFlavour;
    Float_t        jetCharge;
    Float_t        etaMean;
    Float_t        phiMean;
    Float_t        etaEtaMoment;
    Float_t        etaPhiMoment;
    Float_t        phiPhiMoment;
    Float_t        maxDistance;
    Float_t        jetArea;
    Float_t        pileup;
    UChar_t        nPasses;
    UChar_t        nConstituents;

    // CaloJet info
    Float_t        maxEInEmTowers;
    Float_t        maxEInHadTowers;
    Float_t        energyFractionHadronic;
    Float_t        emEnergyFraction;
    Float_t        hadEnergyInHB;
    Float_t        hadEnergyInHO;
    Float_t        hadEnergyInHE;
    Float_t        hadEnergyInHF;
    Float_t        emEnergyInEB;
    Float_t        emEnergyInEE;
    Float_t        emEnergyInHF;
    Float_t        towersArea;
    UChar_t        n90;
    UChar_t        n60;

    // Jet ID info
    Float_t        fHPD;
    Float_t        fRBX;
    Float_t        n90Hits;
    Float_t        fSubDetector1;
    Float_t        fSubDetector2;
    Float_t        fSubDetector3;
    Float_t        fSubDetector4;
    Float_t        restrictedEMF;
    UChar_t        nHCALTowers;
    UChar_t        nECALTowers;
    Float_t        approximatefHPD;
    Float_t        approximatefRBX;
    UChar_t        hitsInN90;
    UChar_t        numberOfHits2RPC;
    UChar_t        numberOfHits3RPC;
    UChar_t        numberOfHitsRPC;

    TVector3       vertex;
    TLorentzVector momentum; // uncorrected momentum
    TLorentzVector detectorP4;

    std::map<TString, Float_t> jecScaleFactors;
    Float_t        jecUncertainty;
  };


  class PFJet {

  public:

    PFJet()  { Init(); }
    ~PFJet() { Init(); }
    void Init();
    void fillRefs(Event const*);

    Bool_t passPuJetIdLoose(unsigned algo)  const { return algo < nPUJetIdAlgorithms ? ( puJetIdFlags[algo] & (1 << 2) ) != 0 : kFALSE; }
    Bool_t passPuJetIdMedium(unsigned algo) const { return algo < nPUJetIdAlgorithms ? ( puJetIdFlags[algo] & (1 << 1) ) != 0 : kFALSE; }
    Bool_t passPuJetIdTight(unsigned algo)  const { return algo < nPUJetIdAlgorithms ? ( puJetIdFlags[algo] & (1 << 0) ) != 0 : kFALSE; }

    // Basic Jet Info
    Int_t          phyDefFlavour;
    Int_t          algDefFlavour;
    Float_t        jetCharge;
    Float_t        etaMean;
    Float_t        phiMean;
    Float_t        etaEtaMoment;
    Float_t        etaPhiMoment;
    Float_t        phiPhiMoment;
    Float_t        maxDistance;
    Float_t        jetArea;
    Float_t        pileup;
    UChar_t        nPasses;
    UChar_t        nConstituents;

    // PF Jet Info
    Float_t        chargedHadronEnergy;
    Float_t        neutralHadronEnergy;
    Float_t        photonEnergy;
    Float_t        electronEnergy;
    Float_t        muonEnergy;
    Float_t        HFHadronEnergy;
    Float_t        HFEMEnergy;
    Float_t        chargedEmEnergy;
    Float_t        chargedMuEnergy;
    Float_t        neutralEmEnergy;
    UChar_t        chargedHadronMultiplicity;
    UChar_t        neutralHadronMultiplicity;
    UChar_t        photonMultiplicity;
    UChar_t        electronMultiplicity;
    UChar_t        muonMultiplicity;
    UChar_t        HFHadronMultiplicity;
    UChar_t        HFEMMultiplicity;
    UChar_t        chargedMultiplicity;
    UChar_t        neutralMultiplicity;

    // Should contain ntuple indices of tracks associated with this jet
    std::vector<UShort_t> tracklist;

    // List of constituent PFParticles
    std::vector<UShort_t> pfParticleList;

    // Pileup Jet Id info
    Float_t        puJetIdDiscriminants[nPUJetIdAlgorithms];
    Int_t          puJetIdFlags[nPUJetIdAlgorithms];

    Float_t        bTagDiscriminators[nBTagDiscriminators];

    TLorentzVector momentum; // uncorrected momentum

    std::map<TString, Float_t> jecScaleFactors;
    Float_t        jecUncertainty;

    std::vector<const Track*> tracks;           //!
    std::vector<const PFParticle*> pfParticles; //!
  };


  class JPTJet {

  public:

    JPTJet()  { Init(); }
    ~JPTJet() { Init(); }
    void Init();
    void fillRefs(Event const*);

    // Basic Jet Info
    Float_t        partonFlavour;
    Float_t        jetCharge;
    Float_t        etaMean;
    Float_t        phiMean;
    Float_t        etaEtaMoment;
    Float_t        etaPhiMoment;
    Float_t        phiPhiMoment;
    Float_t        maxDistance;
    Float_t        jetArea;
    Float_t        pileup;
    UChar_t        nPasses;
    UChar_t        nConstituents;

    Float_t        chargedHadronEnergy;
    Float_t        neutralHadronEnergy;
    Float_t        chargedEmEnergy;
    Float_t        neutralEmEnergy;
    UChar_t        chargedMultiplicity;
    UChar_t        muonMultiplicity;
    UChar_t        elecMultiplicity;
    Float_t        getZSPCor;

    TLorentzVector momentum; // uncorrected momentum

    std::map<TString, Float_t> jecScaleFactors;
    Float_t        jecUncertainty;
  };


  // Class TriggerMap used to be just a typedef of std::map<TString, std::pair<Int_t, UChar_t> >.
  // Interface of the class is thus defined to be the least disruptive to the existing analysis code,
  // i.e. to resemble that of the std::map.

  class TriggerMap {
  public:

    class const_iterator {
      friend class TriggerMap;
    public:
      typedef std::pair<TString, std::pair<UInt_t, Bool_t> > value_type;

      const_iterator();
      const_iterator& operator++();
      const_iterator operator++(int);
      value_type const* operator->() const { return &pair_; }
      value_type const& operator*() const { return pair_; }
      bool operator==(const_iterator const&) const;
      bool operator!=(const_iterator const& _rhs) const { return !operator==(_rhs); }
    protected:
      size_t maxPos_;
      TString const* paths_;
      UInt_t const* prescales_;
      UChar_t const* decisions_;
      UInt_t pos_;
      value_type pair_;
      const_iterator(TriggerMap const&, unsigned);
      void setPair();
    };

    // No non-const iterator functionalities are needed for our purpose
    typedef const_iterator iterator;

    TriggerMap();
    ~TriggerMap() { clear(); }

    void Init();

    std::pair<UInt_t, Bool_t> operator[](TString const&) const;
    Bool_t pass(TString const&) const;
    UInt_t prescale(TString const&) const;

    Bool_t menuExists(TString const&) const;
    void addMenu(TString const&, std::vector<TString> const&);
    void set(TString const&, UInt_t, Bool_t);
    void clear();
    size_t size() const;

    UInt_t* prescales(TString const& _config) const;
    UChar_t* decisions(TString const& _config) const;

    TString currentConfig;

    const_iterator begin() const;
    const_iterator end() const;
    iterator begin();
    iterator end();
    const_iterator find(TString const&) const;
    iterator find(TString const&);
    const_iterator lower_bound(TString const&) const;
    iterator lower_bound(TString const&);
    const_iterator upper_bound(TString const&) const;
    iterator upper_bound(TString const&);

  private:
    void setCache_() const;
    std::map<TString, size_t> size_;
    std::map<TString, TString*> paths_;
    std::map<TString, UInt_t*> prescales_;
    std::map<TString, UChar_t*> decisions_;
    mutable TString cachedConfigName_;
    mutable size_t cachedSize_;
    mutable TString const* cachedPaths_;
    mutable UInt_t* cachedPrescales_;
    mutable UChar_t* cachedDecisions_;
  };

  typedef std::vector<PUSummaryInfo> PUSummaryInfoCollection;
  typedef std::vector<Vertex> VertexCollection;
  typedef std::vector<Track> TrackCollection;
  typedef std::vector<SuperCluster> SuperClusterCollection;
  typedef std::vector<Cluster> ClusterCollection;
  typedef std::vector<Particle> ParticleCollection;
  typedef std::vector<Muon> MuonCollection;
  typedef std::vector<Electron> ElectronCollection;
  typedef std::vector<Photon> PhotonCollection;
  typedef std::vector<CaloJet> CaloJetCollection;
  typedef std::vector<PFJet> PFJetCollection;
  typedef std::vector<JPTJet> JPTJetCollection;
  typedef std::vector<PFParticle> PFParticleCollection;

  class Event {

  public:

    Event();
    ~Event();

    // Initialize members
    void Init();
    void fillRefs();

    // Set up the data members according to the Tree content
    // In case of TChain, files need to be added first
    // Set the second argument to kTRUE (default) for read access
    void bindTree(TTree&, Bool_t = kTRUE);
    void releaseTree(TTree&);
    void setTriggerTable(TString const&, std::vector<std::string> const&, Bool_t = kTRUE);

    Bool_t passMetFilter(UInt_t filterIndex) const { return metFilterBit & (0x1 << filterIndex); }
    // JetMET recommended met filters
    Bool_t passMetFilters() const {
      return passMetFilter(kCSCBeamHalo) && passMetFilter(kHcalNoise) && passMetFilter(kEcalDeadCellTP) &&
        passMetFilter(kHcalLaser) && passMetFilter(kTrackingFailure) && passMetFilter(kEEBadSC);
    }

    UChar_t                                        isRealData;
    UChar_t                                        cosmicFlag;             // empty for now
    UInt_t                                         runNumber;
    UInt_t                                         eventNumber;
    UInt_t                                         luminosityBlockNumber;
    UShort_t                                       bunchCrossing;
    Int_t                                          metFilterBit;

    Float_t                                        avgInsRecLumi;
    Float_t                                        intgRecLumi;
    Float_t                                        rho;                    // from kt6PFJets
    Float_t                                        rhoBarrel;              // from kt6PFJetsRhoBarrelOnly
    Float_t                                        rho25;                  // from kt6PFJetsRho25
                                                   
    TVector3                                       beamSpot;
                                                   
    TriggerMap                                     l1Map;
    TriggerMap                                     hltMap;
                                                   
    VertexCollection                               vertices;
    TrackCollection                                tracks;
    SuperClusterCollection                         superClusters;          // only selected super clusters associated with photons and electrons
    ClusterCollection                              clusters;               // only selected basic clusters associated with super clusters
    PFParticleCollection                           pfParticles;            //

    std::map<TString, susy::MET>                   metMap;                 //
    std::map<TString, susy::MuonCollection>        muons;                  //
    std::map<TString, susy::ElectronCollection>    electrons;              //
    std::map<TString, susy::PhotonCollection>      photons;                //
    std::map<TString, susy::CaloJetCollection>     caloJets;               //
    std::map<TString, susy::PFJetCollection>       pfJets;                 //
    std::map<TString, susy::JPTJetCollection>      jptJets;                // not filled

    // generated information. Valid only for isRealData == 0, i.e. MC
    PUSummaryInfoCollection                        pu;                     // PU summary info
    ParticleCollection                             genParticles;
    std::map<TString, Float_t>                     gridParams;             // pairs of parameter name and value

  private:
    // keep a pointer to each tree that is bound to this object. Bool is true for read-only trees
    std::vector<std::pair<TTree*, Bool_t> >        trees_;                 //
  };

} // namespace susy

#endif
