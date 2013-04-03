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
// $Id: SusyEvent.h,v 1.38 2013/03/30 13:24:12 yiiyama Exp $
//

#ifndef SusyEvent_h
#define SusyEvent_h

#include <vector>
#include <map>
#include <math.h>

#include <TLorentzVector.h>


namespace susy {

  const float etaGapBegin = 1.4442;
  const float etaGapEnd = 1.566;
  const float etaGap = 1.499;
  const float etaMax = 2.5;

  // b-tagging vector positions -- same as susyNtuplizer_cfi.py order
  const unsigned int kTCHE = 0;
  const unsigned int kTCHP = 1;
  const unsigned int kJP = 2;
  const unsigned int kJBP = 3;
  const unsigned int kSSV = 4;
  const unsigned int kCSV = 5;
  const unsigned int kCSVMVA = 6;
  const unsigned int kSE = 7;
  const unsigned int kSM = 8;

  // pileup jet id vector positions -- same as susyNtuplizer_cfi.py order
  const unsigned int kFull = 0;
  const unsigned int kCutBased = 1;
  const unsigned int kSimple = 2;

  enum PFIsoTypes {
    kChargedHadron,
    kNeutralHadron,
    kPhoton,
    nPFIsoTypes
  };

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

  class PUSummaryInfo { /*each PUSummaryInfo object holds information for one BX (early, in time,
                          or late)*/

  public:

    PUSummaryInfo()  { Init(); }
    ~PUSummaryInfo() { Init(); }
    void Init();

    //all info below from https://twiki.cern.ch/twiki/bin/view/CMS/PileupInformation
    /*low_cut = 0.1 GeV, high_cut = 0.5 GeV, tracks summed/counted are TrackingParticles from
      simulation truth info*/
    int numInteractions; //the number of pileup interactions that have been added to the event
    std::vector<float> zPositions; /*the true primary vertex position along the z axis for each
                                     added interaction*/
    std::vector<float> sumPTLowPT; /*the sum of the transverse momentum of the tracks originating
                                     from each interaction, where track pT > low_cut*/
    std::vector<float> sumPTHighPT; /*the sum of the transverse momentum of the tracks originating
                                      from each interaction, where track pT > high_cut*/
    std::vector<int> numTracksLowPT; /*the number of tracks originating from each interaction,
                                       where track pT > low_cut*/
    std::vector<int> numTracksHighPT; /*the number of tracks originating from each interaction,
                                        where track pT > high_cut*/
    std::vector<float> instLumi; //for PU from DataMixer
    std::vector<unsigned int> dataMixerRun;
    std::vector<unsigned int> dataMixerEvt;
    std::vector<unsigned int> dataMixerLumiSection;
    int BX; /*to which bunch crossing does this interaction belong?  New in 3_11_3 and 4_1_3 or
              later*/
    float trueNumInteractions;

  };


  class Particle {

  public:

    Particle()  { Init(); }
    ~Particle() { Init(); }
    void Init();

    UChar_t status;
    Char_t  charge;
    Short_t   motherIndex;
    Int_t   pdgId;
    TVector3 vertex;
    TLorentzVector momentum;

  };


  class CorrMETData {

  public:
    CorrMETData() { Init(); }
    ~CorrMETData() { Init(); }

    void Init();

    Float_t  dmEx;             // for uncorrection, correctedEx - dmEx
    Float_t  dmEy;             // for uncorrection, correctedEy - dmEy
    Float_t  dsumEt;           // for uncorrection, correctedSumEt - dsumEt
    Float_t  dSignificance;    // for uncorrection, correctedSig - dSignificance

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
    TVector3 vertex;
    std::vector<susy::CorrMETData>  mEtCorr;

  };



  class Vertex {

  public:

    Vertex()  { Init(); }
    ~Vertex() { Init(); }
    void Init();
    bool isFake() { return (chi2 == 0 && ndof == 0 && tracksSize == 0); }

    Float_t  chi2;
    Float_t  ndof;
    UShort_t  tracksSize;
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

    Short_t  seedClusterIndex; // index in vector<Cluster> below
    Float_t  energy;
    Float_t  preshowerEnergy;
    Float_t  phiWidth;
    Float_t  etaWidth;
    TVector3 position;
    std::vector<UShort_t> basicClusterIndices;

  };


  class Track {

  public:

    Track()  { Init(); }
    ~Track() { Init(); }
    void Init();

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
    Float_t        error[5]; // qoverp, lambda, phi, dxy, dsz
    TVector3       vertex;
    TLorentzVector momentum;
    std::map<TString,TVector3>  extrapolatedPositions;

  };


  class Photon {

  public:

    Photon()  { Init(); }
    ~Photon() { Init(); }
    void Init();

    // fiducial bits
    bool isEB()           { return (fidBit & (0x1 << 0)); }
    bool isEE()           { return (fidBit & (0x1 << 1)); }
    bool isEBEtaGap()     { return (fidBit & (0x1 << 2)); }
    bool isEBPhiGap()     { return (fidBit & (0x1 << 3)); }
    bool isEERingGap()    { return (fidBit & (0x1 << 4)); }
    bool isEEDeeGap()     { return (fidBit & (0x1 << 5)); }
    bool isEBEEGap()      { return (fidBit & (0x1 << 6)); }
    bool isPF()           { return (fidBit & (0x1 << 7)); }

    Float_t hcalTowerSumEtConeDR04() { return (hcalDepth1TowerSumEtConeDR04+hcalDepth2TowerSumEtConeDR04); }
    Float_t hcalTowerSumEtConeDR03() { return (hcalDepth1TowerSumEtConeDR03+hcalDepth2TowerSumEtConeDR03); }
    Float_t r1x5() { return ((e5x5 > 0) ? e1x5/e5x5 : 0); }
    Float_t r2x5() { return ((e5x5 > 0) ? e2x5/e5x5 : 0); }

    Int_t          fidBit;
    Int_t          nPixelSeeds;
    bool           passelectronveto;
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
    bool           mipIsHalo;

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

    std::pair<double,double> MVAregEnergyAndErr;
    TLorentzVector MVAcorrMomentum;

    TVector3 vertex; // photon vertex when reconstructed.
    TLorentzVector momentum;
    std::map<TString,UChar_t> idPairs;

  };


  class Electron {

  public:

    Electron()  { Init(); }
    ~Electron() { Init(); }
    void Init();

    // fiducial bits
    Bool_t isEB() {        return (fidBit & (0x1 << 0)); }
    Bool_t isEE() {        return (fidBit & (0x1 << 1)); }
    Bool_t isEBEEGap() {   return (fidBit & (0x1 << 2)); }
    Bool_t isEBEtaGap() {  return (fidBit & (0x1 << 3)); }
    Bool_t isEBPhiGap() {  return (fidBit & (0x1 << 4)); }
    Bool_t isEEDeeGap() {  return (fidBit & (0x1 << 5)); }
    Bool_t isEERingGap() { return (fidBit & (0x1 << 6)); }
    Bool_t isEBGap() {     return (isEBEtaGap() || isEBPhiGap()); }
    Bool_t isEEGap() {     return (isEEDeeGap() || isEERingGap()); }
    Bool_t isGap() {       return (isEBGap() || isEEGap() || isEBEEGap()); }

    // boolean variables packed in boolPack
    Bool_t isGsfCtfScPixChargeConsistent() { return (boolPack & (0x1 << 0)); }
    Bool_t isGsfScPixChargeConsistent() {    return (boolPack & (0x1 << 1)); }
    Bool_t isGsfCtfChargeConsistent() {      return (boolPack & (0x1 << 2)); }
    Bool_t ecalDrivenSeed() {                return (boolPack & (0x1 << 3)); }
    Bool_t trackerDrivenSeed() {             return (boolPack & (0x1 << 4)); }
    Bool_t passingCutBasedPreselection() {   return (boolPack & (0x1 << 5)); }
    Bool_t passingMvaPreselection() {        return (boolPack & (0x1 << 6)); }
    Bool_t ambiguous() {                     return (boolPack & (0x1 << 7)); }
    Bool_t isEcalEnergyCorrected() {         return (boolPack & (0x1 << 8)); }
    Bool_t isEnergyScaleCorrected() {        return (boolPack & (0x1 << 9)); }
    Bool_t convFlags() {                     return (boolPack & (0x1 << 10)); }
    Bool_t isPF() {                          return (boolPack & (0x1 << 11)); }
    Bool_t ecalDriven() {                    return (ecalDrivenSeed() && passingCutBasedPreselection()); }

    Float_t hcalOverEcal() { return (hcalDepth1OverEcal + hcalDepth2OverEcal); }
    Float_t dr03HcalTowerSumEt() { return (dr03HcalDepth1TowerSumEt + dr03HcalDepth2TowerSumEt); }
    Float_t dr04HcalTowerSumEt() { return (dr04HcalDepth1TowerSumEt + dr04HcalDepth2TowerSumEt); }

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
    Int_t          mvaStatus;
    Float_t        mva;

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

    // AtVtx, AtCalo
    std::map<TString,TVector3> trackPositions;
    // AtVtx, AtCalo, Out, AtEleClus, AtVtxWithConstraint
    std::map<TString,TLorentzVector> trackMomentums;

    TVector3       vertex;
    TLorentzVector momentum;
    std::map<TString,Float_t> idPairs;

  };



  class Muon {

  public:

    Muon()  { Init(); }
    ~Muon() { Init(); }
    void Init();

    // muon type
    Bool_t isGlobalMuon() {     return (type & (0x1 << 1)); }
    Bool_t isTrackerMuon() {    return (type & (0x1 << 2)); }
    Bool_t isStandAloneMuon() { return (type & (0x1 << 3)); }
    Bool_t isCaloMuon() {       return (type & (0x1 << 4)); }
    Bool_t isPFMuon() {         return (type & (0x1 << 5)); }
    Short_t bestTrackIndex() { switch(bestTrackType){
      case 1: return trackIndex; case 2: return standAloneTrackIndex;
      case 3: return combinedTrackIndex; case 4: return tpfmsTrackIndex;
      case 5: return pickyTrackIndex; case 6: return dytTrackIndex;
      default: return -1;
      } }

    // arbitration type is default
    UChar_t        type;
    UChar_t        bestTrackType;
    UChar_t        nMatches;               // number of muon chambers with matched segments (<= nChambers)
    UChar_t        nValidHits;             // *
    UChar_t        nValidTrackerHits;      // *
    UChar_t        nValidMuonHits;         // *
    UChar_t        nPixelLayersWithMeasurement; // *
    UChar_t        nStripLayersWithMeasurement; // * values from combinedMuon (= globalTrack)
    UChar_t        nChambers;              // number of muon chambers the track traversed through (regardless of segment existence)
    UChar_t        nMatchedStations;       // number of muon stations with matched segments (<= nMatches)
    UChar_t        timeNDof;         // null value implies timing measurement is invalid
    Char_t         timeDirection;
    Float_t        timeAtIp;
    Float_t        timeAtIpError;
    Float_t        caloCompatibility;
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

    std::map<TString, UChar_t> idPairs;

  };


  class CaloJet {

  public:

    CaloJet()  { Init(); }
    ~CaloJet() { Init(); }
    void Init();

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
  };


  class PFJet {

  public:

    PFJet()  { Init(); }
    ~PFJet() { Init(); }
    void Init();

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
    std::vector<Float_t> puJetIdDiscriminants;
    std::vector<Int_t> puJetIdFlags;

    Bool_t passPuJetIdLoose(unsigned int type) const { return type < puJetIdFlags.size() ? ( puJetIdFlags[type] & (1 << 2) ) != 0 : false ; }
    Bool_t passPuJetIdMedium(unsigned int type) const { return type < puJetIdFlags.size() ? ( puJetIdFlags[type] & (1 << 1) ) != 0 : false ; }
    Bool_t passPuJetIdTight(unsigned int type) const { return type < puJetIdFlags.size() ? ( puJetIdFlags[type] & (1 << 0) ) != 0 : false ; }

    // IMPORTANT: This vector of float stores btag-discriminator variables from various collections
    // which defined in susyNtuplizer_cfi.py file. The order of variables are strongly
    // dependent on the order of collections in python file.
    std::vector<Float_t> bTagDiscriminators;

    TVector3       vertex;
    TLorentzVector momentum; // uncorrected momentum

    std::map<TString, Float_t> jecScaleFactors;
  };


  class JPTJet {

  public:

    JPTJet()  { Init(); }
    ~JPTJet() { Init(); }
    void Init();

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

    TVector3       vertex;
    TLorentzVector momentum; // uncorrected momentum

    std::map<TString, Float_t> jecScaleFactors;
  };


  class PFParticle {

  public:
    PFParticle() { Init(); }
    ~PFParticle() { Init(); }
    void Init();

    Int_t pdgId;
    Char_t  charge;
    Float_t ecalEnergy;
    Float_t rawEcalEnergy;
    Float_t hcalEnergy;
    Float_t rawHcalEnergy;
    Float_t pS1Energy;
    Float_t pS2Energy;

    TVector3 vertex;
    TVector3 positionAtECALEntrance;
    TLorentzVector momentum;

  };



  typedef std::vector<susy::Electron> ElectronCollection;
  typedef std::vector<susy::Photon> PhotonCollection;
  typedef std::vector<susy::CaloJet> CaloJetCollection;
  typedef std::vector<susy::PFJet> PFJetCollection;
  typedef std::vector<susy::JPTJet> JPTJetCollection;
  typedef std::vector<susy::PFParticle> PFParticleCollection;
  typedef std::vector<susy::PUSummaryInfo> PUSummaryInfoCollection;
  typedef std::vector<susy::Muon> MuonCollection;
  typedef std::map<TString, std::pair<Int_t, UChar_t> > TriggerMap;



  class Event {

  public:

    Event()  { Init(); }
    ~Event() { Init(); }

    // Initialize members
    void Init();

    bool passCSCBeamHalo()     const { return metFilterBit & (0x1 << 0); }
    bool passHcalNoise()       const { return metFilterBit & (0x1 << 1); }
    bool passEcalDeadCellTP()  const { return metFilterBit & (0x1 << 2); }
    bool passEcalDeadCellBE()  const { return metFilterBit & (0x1 << 3); }
    bool passHcalLaser()       const { return metFilterBit & (0x1 << 4); }
    bool passTrackingFailure() const { return metFilterBit & (0x1 << 5); }
    bool passEEBadSC()         const { return metFilterBit & (0x1 << 6); }

    bool passEERingOfFire()     const { return metFilterBit_2 & (0x1 << 0); }
    bool passInconsistentMuon() const { return metFilterBit_2 & (0x1 << 1); }
    bool passGreedyMuon()       const { return metFilterBit_2 & (0x1 << 2); }
    bool passHcalLaser2012()    const { return metFilterBit_2 & (0x1 << 3); }
    bool passEcalLaserCorr()    const { return metFilterBit_2 & (0x1 << 4); }
    bool passManyStripClus()    const { return metFilterBit_2 & (0x1 << 5); }
    bool passTooManyStripClus() const { return metFilterBit_2 & (0x1 << 6); }
    bool passLogErrorTooManyClusters() const { return metFilterBit_2 & (0x1 << 7); }

    bool passTrkPOGFilters() const { return passManyStripClus() && passTooManyStripClus() && passLogErrorTooManyClusters(); }

    // JetMET recommended met filters
    bool passMetFilters() const { return passCSCBeamHalo() && passHcalNoise() &&
    passEcalDeadCellTP() && passHcalLaser() && passTrackingFailure() && passEEBadSC(); }

    // Members are made as public intentionally for easy access

    UChar_t                                     isRealData;
    Int_t                                       runNumber;
    ULong_t                                     eventNumber;
    Int_t                                       luminosityBlockNumber;
    Int_t                                       bunchCrossing;
    Float_t                                     avgInsRecLumi;
    Float_t                                     intgRecLumi;
    UChar_t                                     cosmicFlag; // empty for now
    Float_t                                     rho; // from kt6PFJets
    Float_t                                     rhoBarrel; // from kt6PFJetsRhoBarrelOnly
    Float_t                                     rho25; // from kt6PFJetsRho25
    Int_t                                       metFilterBit;
    Int_t                                       metFilterBit_2;

    TVector3                                    beamSpot;

    susy::TriggerMap                            l1Map;  // <name, <prescale, bit> >
    susy::TriggerMap                            hltMap; // <name, <prescale, bit> >
    std::map<TString,susy::MET>                 metMap;

    std::vector<susy::Vertex>                   vertices;
    std::vector<susy::Track>                    tracks;          // only selected tracks associated with objects directly and photons with dR<0.4
    std::vector<susy::SuperCluster>             superClusters;   // only selected super clusters associated with objects
    std::vector<susy::Cluster>                  clusters;        // only selected basic clusters associated with super clusters
    std::map<TString,susy::MuonCollection>      muons;
    std::map<TString,susy::ElectronCollection>  electrons;
    std::map<TString,susy::PhotonCollection>    photons;
    std::map<TString,susy::CaloJetCollection>   caloJets;
    std::map<TString,susy::PFJetCollection>     pfJets;
    std::map<TString,susy::JPTJetCollection>    jptJets;         // dropped for 2011B analysis
    std::map<TString,susy::PFParticleCollection> pfParticles;

    // optional collections
    std::vector<susy::Track>                    generalTracks;   // not stored by default

    // generated information. Valid only for isRealData == 0, i.e. MC
    susy::PUSummaryInfoCollection               pu; //PU summary info
    std::vector<TVector3>                       simVertices; // Geant vertex, primary only, dropped for 2011B analysis
    std::vector<susy::Particle>                 genParticles;
    std::map<TString, Float_t>                  gridParams; // pairs of parameter name and value

  };




} // namespace susy

#endif
