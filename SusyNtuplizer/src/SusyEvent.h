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
// $Id: SusyEvent.h,v 1.13 2011/06/03 16:58:47 dwjang Exp $
//

#ifndef SusyEvent_h
#define SusyEvent_h

#include <vector>
#include <map>
#include <math.h>

#include <TLorentzVector.h>


namespace susy {

  const float etaGapBegin = 1.442;
  const float etaGapEnd = 1.556;
  const float etaGap = 1.499;
  const float etaMax = 2.5;


  class Particle {

  public:

    Particle()  { Init(); }
    ~Particle() { Init(); }
    void Init();

    UChar_t status;
    Int_t   motherId;
    Int_t   pdgId;
    Char_t  charge;
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
    UChar_t  tracksSize;
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
    std::vector<Int_t> basicClusterIndices;
    
  };


  class Track {

  public:

    Track()  { Init(); }
    ~Track() { Init(); }
    void Init();

    // derived quantities
    float normChi2() const { return (ndof != 0) ? chi2/ndof : chi2*1e6; }
    float qoverp() const { return charge/momentum.P(); }
    float lambda() const { return M_PI/2 - momentum.Theta(); }
    float dsz() const { return vertex.Z()*momentum.Pt()/momentum.P() - (vertex.X()*momentum.Px()+vertex.Y()*momentum.Py())/momentum.Pt() * momentum.Pz()/momentum.P(); }
    float dz() const { return vertex.Z() - (vertex.X()*momentum.Px()+vertex.Y()*momentum.Py())/momentum.Pt() * (momentum.Pz()/momentum.Pt()); }
    float dxy() const { return (-vertex.X()*momentum.Py() + vertex.Y()*momentum.Px())/momentum.Pt(); }
    float d0() const { return -dxy(); }
    float phi() const { return momentum.Phi(); }
    bool loose() const {         return (quality & ( 0x1 << 0)); }
    bool tight() const {         return (quality & ( 0x1 << 1)); }
    bool highPurity() const {    return (quality & ( 0x1 << 2)); }
    bool confirmed() const {     return (quality & ( 0x1 << 3)); }
    bool goodIterative() const { return (confirmed() || highPurity()); }

    Int_t          algorithm;
    Int_t          quality;
    UChar_t        numberOfValidHits;
    UChar_t        numberOfValidTrackerHits;
    UChar_t        numberOfValidMuonHits;
    UChar_t        numberOfValidPixelHits;
    UChar_t        numberOfValidStripHits;
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
    Float_t        hadronicOverEm;
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
    Float_t        r9;

    Float_t        ecalRecHitSumEtConeDR04;
    Float_t        hcalDepth1TowerSumEtConeDR04;
    Float_t        hcalDepth2TowerSumEtConeDR04;
    Float_t        trkSumPtSolidConeDR04;
    Float_t        trkSumPtHollowConeDR04;
    UChar_t        nTrkSolidConeDR04;
    UChar_t        nTrkHollowConeDR04;

    Float_t        ecalRecHitSumEtConeDR03;
    Float_t        hcalDepth1TowerSumEtConeDR03;
    Float_t        hcalDepth2TowerSumEtConeDR03;
    Float_t        trkSumPtSolidConeDR03;
    Float_t        trkSumPtHollowConeDR03;
    UChar_t        nTrkSolidConeDR03;
    UChar_t        nTrkHollowConeDR03;

    Float_t        chargedHadronIso;
    Float_t        neutralHadronIso;
    Float_t        photonIso;

    Float_t        seedTime; // seed timing

    // cluster shape variables (barrel only)
    Float_t        sMaj;
    Float_t        sMin;
    Float_t        alpha;
    Float_t        roundness;
    Float_t        angle;

    // Conversion info
    Float_t        convDist;
    Float_t        convDcot;
    Float_t        convVtxChi2;
    UChar_t        convVtxNdof;
    TVector3       convVertex;

    Short_t        superClusterIndex;
    Float_t        superClusterPreshowerEnergy;
    Float_t        superClusterPhiWidth;
    Float_t        superClusterEtaWidth;
    TVector3       caloPosition;

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
    bool isEB() {        return (fidBit & (0x1 << 0)); }
    bool isEE() {        return (fidBit & (0x1 << 1)); }
    bool isEBEEGap() {   return (fidBit & (0x1 << 2)); }
    bool isEBEtaGap() {  return (fidBit & (0x1 << 3)); }
    bool isEBPhiGap() {  return (fidBit & (0x1 << 4)); }
    bool isEEDeeGap() {  return (fidBit & (0x1 << 5)); }
    bool isEERingGap() { return (fidBit & (0x1 << 6)); }
    bool isEBGap() {     return (isEBEtaGap() || isEBPhiGap()); }
    bool isEEGap() {     return (isEEDeeGap() || isEERingGap()); }
    bool isGap() {       return (isEBGap() || isEEGap() || isEBEEGap()); }

    // boolean variables packed in boolPack
    bool isGsfCtfScPixChargeConsistent() { return (boolPack & (0x1 << 0)); }
    bool isGsfScPixChargeConsistent() {    return (boolPack & (0x1 << 1)); }
    bool isGsfCtfChargeConsistent() {      return (boolPack & (0x1 << 2)); }
    bool ecalDrivenSeed() {                return (boolPack & (0x1 << 3)); }
    bool trackerDrivenSeed() {             return (boolPack & (0x1 << 4)); }
    bool passingCutBasedPreselection() {   return (boolPack & (0x1 << 5)); }
    bool passingMvaPreselection() {        return (boolPack & (0x1 << 6)); }
    bool ambiguous() {                     return (boolPack & (0x1 << 7)); }
    bool isEcalEnergyCorrected() {         return (boolPack & (0x1 << 8)); }
    bool isEnergyScaleCorrected() {        return (boolPack & (0x1 << 9)); }
    bool convFlags() {                     return (boolPack & (0x1 << 10)); }
    bool isPF() {                          return (boolPack & (0x1 << 11)); }
    bool ecalDriven() {                    return (ecalDrivenSeed() && passingCutBasedPreselection()); }

    Float_t hcalOverEcal() { return (hcalDepth1OverEcal + hcalDepth2OverEcal); }
    Float_t dr03HcalTowerSumEt() { return (dr03HcalDepth1TowerSumEt + dr03HcalDepth2TowerSumEt); }
    Float_t dr04HcalTowerSumEt() { return (dr04HcalDepth1TowerSumEt + dr04HcalDepth2TowerSumEt); }

    Int_t          fidBit;
    Int_t          boolPack;
    Int_t          scPixCharge;

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
    Float_t        e1x5;
    Float_t        e2x5Max;
    Float_t        e5x5;
    Float_t        hcalDepth1OverEcal;          // hadronic energy on depth1 / em enrgy
    Float_t        hcalDepth2OverEcal;          // hadronic energy on depth2 / em enrgy

    Float_t        dr03TkSumPt;
    Float_t        dr03EcalRecHitSumEt;
    Float_t        dr03HcalDepth1TowerSumEt;
    Float_t        dr03HcalDepth2TowerSumEt;

    Float_t        dr04TkSumPt;
    Float_t        dr04EcalRecHitSumEt;
    Float_t        dr04HcalDepth1TowerSumEt;
    Float_t        dr04HcalDepth2TowerSumEt;

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
    bool isGlobalMuon() {     return (type & (0x1 << 1)); }
    bool isTrackerMuon() {    return (type & (0x1 << 2)); }
    bool isStandAloneMuon() { return (type & (0x1 << 3)); }
    bool isCaloMuon() {       return (type & (0x1 << 4)); }

    UChar_t        type;
    UChar_t        nMatches;
    UChar_t        nValidHits;
    UChar_t        nValidTrackerHits;
    UChar_t        nValidMuonHits;
    UChar_t        nChambers;
    UChar_t        timeNDof;
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

    Short_t        trackIndex;             // tracker only
    Short_t        standAloneTrackIndex;   // muon detector only
    Short_t        combinedTrackIndex;     // combined
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
  typedef std::map<TString, std::pair<Int_t, UChar_t> > TriggerMap;



  class Event {

  public:

    Event()  { Init(); }
    ~Event() { Init(); }

    // Initialize members
    void Init();

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

    TVector3                                    beamSpot;

    susy::TriggerMap                            l1Map;  // <name, <prescale, bit> >
    susy::TriggerMap                            hltMap; // <name, <prescale, bit> >
    std::map<TString,susy::MET>                 metMap;

    std::vector<susy::Vertex>                   vertices;
    std::vector<susy::Track>                    tracks;          // only selected tracks associated with objects directly and photons with dR<0.4
    std::vector<susy::SuperCluster>             superClusters;   // only selected super clusters associated with objects
    std::vector<susy::Cluster>                  clusters;        // only selected basic clusters associated with super clusters
    std::vector<susy::Muon>                     muons;
    std::map<TString,susy::ElectronCollection>  electrons;
    std::map<TString,susy::PhotonCollection>    photons;
    std::map<TString,susy::CaloJetCollection>   caloJets;
    std::map<TString,susy::PFJetCollection>     pfJets;
    std::map<TString,susy::JPTJetCollection>    jptJets;
    std::map<TString,susy::PFParticleCollection> pfParticles;

    // optional collections
    std::vector<susy::Track>                    generalTracks;   // not stored by default

    // generated information. Valid only for isRealData == 0, i.e. MC
    std::vector<TVector3>                       simVertices; // Geant vertex, primary only
    std::vector<susy::Particle>                 genParticles;
    std::map<TString, Float_t>                  gridParams; // pairs of parameter name and value

  };




} // namespace susy

#endif
