// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyNtuplizer
// 
/**\class SusyNtuplizer SusyNtuplizer.cc SusyAnalysis/SusyNtuplizer/src/SusyNtuplizer.cc

 Description: Ntuple maker for SUSY analysis

 Implementation:
     Putting all header files in src area violates CMS coding rule,
     but it is convenient for later use in standalone mode.
*/
//
// Original Author:  Dongwook Jang
// $Id: SusyNtuplizer.cc,v 1.1 2011/03/24 23:46:27 dwjang Exp $
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"


// simple geometry
#include "RecoParticleFlow/PFProducer/interface/PFGeometry.h"

// for track extrapolation
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

// for conversion finder related
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"


// for ecal rechit related
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Calibration/IsolatedParticles/interface/eECALMatrix.h"


// system include files
#include <memory>
#include <string>
#include <vector>
#include <algorithm>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include "SusyEvent.h"

const int nJECFlavours = 5;
const std::string jecFlavours[nJECFlavours] = {"gluon", "uds", "charm", "bottom", "none"};


// Class definition
class SusyNtuplizer : public edm::EDAnalyzer {

public:
  explicit SusyNtuplizer(const edm::ParameterSet&);
  ~SusyNtuplizer();

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  std::string name() { return "SusyNtuplizer"; }
  void fillGenInfos(const edm::Event&, const edm::EventSetup&);
  void fillTrack(const reco::TrackRef& in, susy::Track& out);
  void fillTrack(const reco::GsfTrackRef& in, susy::Track& out);
  void fillCluster(const reco::CaloClusterPtr& in, susy::Cluster& out);
  void fillCluster(const reco::SuperClusterRef& in, susy::SuperCluster& out, int& basicClusterIndex); // basicClusterIndex will be incremented here
  void fillParticle(const reco::GenParticle* in, susy::Particle& out, int igen);
  void fillParticle(const reco::PFCandidateRef& in, susy::Particle& out);
  void fillExtrapolations(const reco::TrackRef& ttk, std::map<TString,TVector3>& positions);
  bool sameGenParticles(const reco::GenParticle* gp, const reco::GenParticle* gp2);

  // ----------member data ---------------------------

  // InputTags

  edm::InputTag lumiSummaryTag_;
  edm::InputTag l1GTReadoutTag_;
  edm::InputTag l1GTObjectMapTag_;
  edm::InputTag hltCollectionTag_;
  edm::InputTag vtxCollectionTag_;
  edm::InputTag trackCollectionTag_;
  edm::InputTag muonCollectionTag_;
  std::vector<std::string> muonIDCollectionTags_;
  edm::InputTag electronCollectionTag_;
  std::vector<std::string> electronIDCollectionTags_;
  edm::InputTag photonCollectionTag_;
  std::vector<std::string> photonIDCollectionTags_;
  edm::InputTag pfElectronCollectionTag_;
  edm::InputTag pfParticleCollectionTag_;
  edm::InputTag tauCollectionTag_;
  edm::InputTag genCollectionTag_;
  std::vector<std::string> caloJetCollectionTags_;
  std::vector<std::string> pfJetCollectionTags_;
  std::vector<std::string> jptJetCollectionTags_;
  std::vector<std::string> metCollectionTags_;

  edm::ESHandle<MagneticField> magneticField_;
  PropagatorWithMaterial* propagator_;
  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder_;

  PFGeometry pfGeom_;

  // debugLevel
  // 0 : default (no printout from this module)
  // 1 : minimal (function level printing)
  // 2 : print the size of collections
  // 3 : print the values of objects in the collection
  int debugLevel_;

  // flag for storing generated informations in ntuples
  // default : false
  bool storeGenInfos_;

  // flag for storing tau collection in the ntuple
  // default : false
  bool storeTauColl_;

  // flag for storing PF electron collection in the ntuple
  bool storePFElectronColl_;
  // default : true

  // flag for storing PF collection in the ntuple
  bool storePFParticleColl_;
  // default : true

  // input RECO mode
  // false : default - reading from AOD
  //         no extrapolation info will be saved in ntuples
  // true : reading from RECO
  bool recoMode_;

  std::string outputFileName_;

  susy::Event* susyEvent_;
  TTree*       susyTree_;

};


// Constructor - passing parameters, memory allocation to ntuple variables
SusyNtuplizer::SusyNtuplizer(const edm::ParameterSet& iConfig) {

  lumiSummaryTag_            = iConfig.getParameter<edm::InputTag>("lumiSummaryTag");
  l1GTReadoutTag_            = iConfig.getParameter<edm::InputTag>("l1GTReadoutTag");
  l1GTObjectMapTag_          = iConfig.getParameter<edm::InputTag>("l1GTObjectMapTag");
  hltCollectionTag_          = iConfig.getParameter<edm::InputTag>("hltCollectionTag");
  vtxCollectionTag_          = iConfig.getParameter<edm::InputTag>("vtxCollectionTag");
  trackCollectionTag_        = iConfig.getParameter<edm::InputTag>("trackCollectionTag");
  muonCollectionTag_         = iConfig.getParameter<edm::InputTag>("muonCollectionTag");
  muonIDCollectionTags_      = iConfig.getParameter<std::vector<std::string> >("muonIDCollectionTags");
  electronCollectionTag_     = iConfig.getParameter<edm::InputTag>("electronCollectionTag");
  electronIDCollectionTags_  = iConfig.getParameter<std::vector<std::string> >("electronIDCollectionTags");
  photonCollectionTag_       = iConfig.getParameter<edm::InputTag>("photonCollectionTag");
  photonIDCollectionTags_    = iConfig.getParameter<std::vector<std::string> >("photonIDCollectionTags");
  pfElectronCollectionTag_   = iConfig.getParameter<edm::InputTag>("pfElectronCollectionTag");
  pfParticleCollectionTag_   = iConfig.getParameter<edm::InputTag>("pfParticleCollectionTag");
  tauCollectionTag_          = iConfig.getParameter<edm::InputTag>("tauCollectionTag");
  genCollectionTag_          = iConfig.getParameter<edm::InputTag>("genCollectionTag");
  caloJetCollectionTags_     = iConfig.getParameter<std::vector<std::string> >("caloJetCollectionTags");
  pfJetCollectionTags_       = iConfig.getParameter<std::vector<std::string> >("pfJetCollectionTags");
  jptJetCollectionTags_      = iConfig.getParameter<std::vector<std::string> >("jptJetCollectionTags");
  metCollectionTags_         = iConfig.getParameter<std::vector<std::string> >("metCollectionTags");

  debugLevel_ = iConfig.getParameter<int>("debugLevel");
  storeGenInfos_ = iConfig.getParameter<bool>("storeGenInfos");
  storeTauColl_ = iConfig.getParameter<bool>("storeTauColl");
  storePFElectronColl_ = iConfig.getParameter<bool>("storePFElectronColl");
  storePFParticleColl_ = iConfig.getParameter<bool>("storePFParticleColl");
  recoMode_ = iConfig.getParameter<bool>("recoMode");
  outputFileName_ = iConfig.getParameter<std::string>("outputFileName");

  susyEvent_ = new susy::Event;

  TFile* outF = new TFile(outputFileName_.c_str(),"RECREATE");
  //  TH1::AddDirectory(kTRUE);
  if(outF) {
    susyTree_ = new TTree("susyTree","SUSY Event");
    susyTree_->Branch("susyEvent","susy::Event",&susyEvent_);
    susyTree_->SetAutoSave(10000000);     // 10 MB
    //    susyTree_->SetMaxTreeSize(500000000); // 500 MB
  }

  if(debugLevel_ > 0) std::cout << name() << " : ctor" << std::endl;

}


// Destructor - free up memory used
SusyNtuplizer::~SusyNtuplizer() {

  if(debugLevel_ > 0) std::cout << name() << " : dtor" << std::endl;

  if(susyEvent_) delete susyEvent_;

}

// ------------ method called once each job just before starting event loop  ------------
// Open output file and setup tree structures
void SusyNtuplizer::beginJob() {

  if(debugLevel_ > 0) std::cout << name() << " : beginJob" << std::endl;

}

// ------------ method called once each job just after ending the event loop  ------------
// Write the tree to the output
void SusyNtuplizer::endJob() {

  if(debugLevel_ > 0) std::cout << name() << " : endJob" << std::endl;

  susyTree_->GetCurrentFile()->cd();
  //  susyTree_->GetCurrentFile()->Write();
  susyTree_->Write();
  susyTree_->GetCurrentFile()->Close();

}


// ------------ method called to for each event  ------------
// fill the tree variables
void SusyNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if(debugLevel_ > 0) std::cout << name() << " : initial event setup" << std::endl;

  int trackIndex = 0;
  int clusterIndex = 0;
  int superClusterIndex = 0;

  if(recoMode_){
    try {
      iSetup.get<IdealMagneticFieldRecord>().get(magneticField_);
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "IdealMagneticFieldRecord is not available!!! " << e.what();
    }

    try {
      iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",transientTrackBuilder_);
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "TransientTrackRecord is not available!!! " << e.what();
    }

    propagator_ = new PropagatorWithMaterial(alongMomentum,0.000511,&*magneticField_);
  }

  susyEvent_->Init();

  if(debugLevel_ > 0) std::cout << name() << ", fill event info" << std::endl;

  susyEvent_->isRealData = iEvent.isRealData() ? 1 : 0;
  susyEvent_->runNumber = iEvent.id().run();
  susyEvent_->eventNumber = iEvent.id().event();
  susyEvent_->luminosityBlockNumber = iEvent.getLuminosityBlock().luminosityBlock();

//   edm::Handle<LumiSummary> lsH;

//   try {
//     iEvent.getByLabel(lumiSummaryTag_, lsH);
//     susyEvent_->avgInsRecLumi = lsH->avgInsRecLumi();
//     susyEvent_->intgRecLumi = lsH->intgRecLumi();
//   }
//   catch(cms::Exception& e) {
//     edm::LogError(name()) << "LumiSummary is not available!!! " << e.what();
//   }

  if(debugLevel_ > 1) std::cout << name() << ", run " << iEvent.id().run()
				<< ", event " << iEvent.id().event()
				<< ", isRealData " << iEvent.isRealData()
				<< ", lumiBlock " << iEvent.getLuminosityBlock().luminosityBlock() << std::endl;

  if(debugLevel_ > 0) std::cout << name() << ", fill L1 map" << std::endl;

  edm::Handle<L1GlobalTriggerReadoutRecord>   gtRecord;
  edm::Handle<L1GlobalTriggerObjectMapRecord> gtOMRec;
  
  try {
    iEvent.getByLabel(l1GTReadoutTag_, gtRecord);
    iEvent.getByLabel(l1GTObjectMapTag_, gtOMRec);

    const DecisionWord dWord = gtRecord->decisionWord();
 
    std::vector<L1GlobalTriggerObjectMap>::const_iterator iter_begin = gtOMRec->gtObjectMap().begin();
    std::vector<L1GlobalTriggerObjectMap>::const_iterator iter_end = gtOMRec->gtObjectMap().end();
    std::vector<L1GlobalTriggerObjectMap>::const_iterator iter = iter_begin;

    for( ; iter != iter_end; iter++) {
      susyEvent_->l1Map[TString(iter->algoName().c_str())] = UChar_t(dWord[iter->algoBitNumber()]);
      if(debugLevel_ > 1) std::cout << iter->algoName() << " : " << dWord[iter->algoBitNumber()] << std::endl;
    }
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "L1TriggerCollection is not available!!! " << e.what();
  }


  if(debugLevel_ > 0) std::cout << name() << ", fill HLT map" << std::endl;

  edm::Handle<edm::TriggerResults> hltH;
  try {
    iEvent.getByLabel(hltCollectionTag_,hltH);
    int nHlt = hltH->size();
    const edm::TriggerNames& hltTriggerNames = iEvent.triggerNames(*hltH);
    if(nHlt != int(hltTriggerNames.size())) edm::LogError(name()) << "TriggerPathName size mismatches !!! ";
    for(int i=0; i<nHlt; i++) {
      susyEvent_->hltMap[TString(hltTriggerNames.triggerName(i).c_str())] = UChar_t(hltH->accept(i));
      if(debugLevel_ > 1) std::cout << hltTriggerNames.triggerName(i) << " : " << hltH->accept(i) << std::endl;
    }
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "TriggerResults is not available!!! " << e.what();
  }


  if(debugLevel_ > 0) std::cout << name() << ", fill beam spot" << std::endl;

  edm::Handle<reco::BeamSpot> bsh;
  try {
    iEvent.getByType(bsh);
    susyEvent_->beamSpot.SetXYZ(bsh->position().x(),bsh->position().y(),bsh->position().z());
    if(debugLevel_ > 1) std::cout << "beamSpot : " << bsh->position().x() << ", "
				  << bsh->position().y() << ", " << bsh->position().z() << std::endl;
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "BeamSpot is not available!!! " << e.what();
  }


  if(debugLevel_ > 0) std::cout << name() << ", fill vertex - the first entry is the primary vertex" << std::endl;

  edm::Handle<reco::VertexCollection> vtxH;
  try {
    iEvent.getByLabel(vtxCollectionTag_,vtxH);
    for(reco::VertexCollection::const_iterator it = vtxH->begin();
	it != vtxH->end(); it++){
      TVector3 vtx(it->x(),it->y(),it->z());
      susyEvent_->vertices.push_back(vtx);
      if(debugLevel_ > 1) std::cout << "vertex : " << it->x() << ", " << it->y() << ", " << it->z() << std::endl;
    }
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "VertexCollection is not available!!! " << e.what();
  }



  if(debugLevel_ > 0) std::cout << name() << ", get general track collection" << std::endl;

  edm::Handle<reco::TrackCollection> trackH;
  try {
    iEvent.getByLabel(trackCollectionTag_,trackH);
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "TrackCollection is not available!!! " << e.what();
  }


  // This is issued because ConversionFinder requires it.
  edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
  edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
  edm::ESHandle<CaloGeometry> cgH;
  const CaloGeometry* geo = 0;
  edm::ESHandle<CaloTopology> ctH;
  const CaloTopology* caloTopology = 0;

  if(recoMode_) {
    
    if(debugLevel_ > 0) std::cout << name() << ", get ecal rechits" << std::endl;
    
    
    try {
      iEvent.getByLabel("ecalRecHit","EcalRecHitsEB",barrelRecHitsHandle);
      iEvent.getByLabel("ecalRecHit","EcalRecHitsEE",endcapRecHitsHandle);
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "EcalRecHitCollection is not available!!! " << e.what();
    }
    
    
    if(debugLevel_ > 0) std::cout << name() << ", get calo geometry record." << std::endl;
    
    try{
      iSetup.get<CaloGeometryRecord>().get(cgH);
      geo = cgH.product();
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "CaloGeometryRecord is not available!!! " << e.what();
    }

    
    if(debugLevel_ > 0) std::cout << name() << ", get calo topology record." << std::endl;
    
    try {
      iSetup.get<CaloTopologyRecord>().get(ctH); 
      caloTopology = ctH.product();
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "CaloTopologyRecord is not available!!! " << e.what();
    }

  }// if(recoMode_


  if(debugLevel_ > 0) std::cout << name() << ", fill all kinds of met collections" << std::endl;


  // Met from several met collections
  int nMetColl = metCollectionTags_.size();

  for(int iMet=0; iMet<nMetColl; iMet++) {

    //    edm::Handle<edm::View<pat::MET> > metH;
    edm::Handle<edm::View<reco::MET> > metH;
    try {
      iEvent.getByLabel(edm::InputTag(metCollectionTags_[iMet]),metH);
      if(debugLevel_ > 1){
	std::cout << "size of MET coll : " << metH->size() << std::endl;
      }
      susy::MET susyMet;
      susyMet.mEt.Set(metH->begin()->p4().px(),metH->begin()->p4().py());
      susyMet.vertex.SetXYZ(metH->begin()->vx(),metH->begin()->vy(),metH->begin()->vz());
      susyMet.sumEt = metH->begin()->sumEt();
      susyMet.significance = metH->begin()->significance();
      if(debugLevel_ > 2) {
	std::cout << "met, metX, metY, sumEt, significance : " << susyMet.mEt.Mod() << ", " << susyMet.mEt.X() << ", "
		  << susyMet.mEt.Y() << ", " << susyMet.sumEt << ", " << susyMet.significance << std::endl;
      }
      int nCorr = metH->begin()->mEtCorr().size();
      for(int i=0; i<nCorr; i++){
	susy::CorrMETData corr;
	corr.dmEx          = metH->begin()->dmEx()[i];
	corr.dmEy          = metH->begin()->dmEy()[i];
	corr.dsumEt        = metH->begin()->dsumEt()[i];
	corr.dSignificance = metH->begin()->dSignificance()[i];
	susyMet.mEtCorr.push_back(corr);
	if(debugLevel_ > 2) {
	  std::cout << "dmEx, dmEy, dsumEt, dSignificance : "
		    << corr.dmEx << ", " << corr.dmEy << ", "
		    << corr.dsumEt << ", " << corr.dSignificance << std::endl;
	}
      }//for
      susyEvent_->metMap[TString(metCollectionTags_[iMet].c_str())] = susyMet;
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "pat::MET " << metCollectionTags_[iMet] << " is not available!!! " << e.what();
    }

  }// for iMet


  if(debugLevel_ > 0) std::cout << name() << ", fill photon" << std::endl;
  
  const int nPhoIDC = photonIDCollectionTags_.size();
  std::vector< const edm::ValueMap<Bool_t>* > phoIds;
  
  for(int i=0; i<nPhoIDC; i++) {
    edm::Handle<edm::ValueMap<Bool_t> > phoIDCH;
    try {
      iEvent.getByLabel("PhotonIDProd",photonIDCollectionTags_[i], phoIDCH);
      phoIds.push_back( phoIDCH.product() );
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << photonIDCollectionTags_[i] << " is not available!!! " << e.what();
    }
  }

  edm::Handle<reco::PhotonCollection> photonH;
  try {
    iEvent.getByLabel(photonCollectionTag_,photonH);
    if(debugLevel_ > 1) std::cout << "size of PhotonCollection : " << photonH->size() << std::endl;
    int ipho = 0;
    for(reco::PhotonCollection::const_iterator it = photonH->begin();
	it != photonH->end(); it++){

      reco::PhotonRef phoRef(photonH,ipho++);

      susy::Photon pho;

      // pack fiducial bits
      pho.fidBit |= (it->isEB()        << 0);
      pho.fidBit |= (it->isEE()        << 1);
      pho.fidBit |= (it->isEBEtaGap()  << 2);
      pho.fidBit |= (it->isEBPhiGap()  << 3);
      pho.fidBit |= (it->isEERingGap() << 4);
      pho.fidBit |= (it->isEEDeeGap()  << 5);
      pho.fidBit |= (it->isEBEEGap()   << 6);

      pho.nPixelSeeds                       = it->electronPixelSeeds().size();
      pho.hadronicOverEm                    = it->hadronicOverEm();
      pho.hadronicDepth1OverEm              = it->hadronicDepth1OverEm();
      pho.hadronicDepth2OverEm              = it->hadronicDepth2OverEm();

      pho.e1x5                              = it->e1x5();
      pho.e2x5                              = it->e2x5();
      pho.e3x3                              = it->e3x3();
      pho.e5x5                              = it->e5x5();
      pho.maxEnergyXtal                     = it->maxEnergyXtal();
      pho.sigmaEtaEta                       = it->sigmaEtaEta();
      pho.sigmaIetaIeta                     = it->sigmaIetaIeta();
      pho.r1x5                              = it->r1x5();
      pho.r2x5                              = it->r2x5();
      pho.r9                                = it->r9();

      pho.ecalRecHitSumEtConeDR04           = it->ecalRecHitSumEtConeDR04();
      pho.hcalTowerSumEtConeDR04            = it->hcalTowerSumEtConeDR04();
      pho.hcalDepth1TowerSumEtConeDR04      = it->hcalDepth1TowerSumEtConeDR04();
      pho.hcalDepth2TowerSumEtConeDR04      = it->hcalDepth2TowerSumEtConeDR04();
      pho.trkSumPtSolidConeDR04             = it->trkSumPtSolidConeDR04();
      pho.trkSumPtHollowConeDR04            = it->trkSumPtHollowConeDR04();
      pho.nTrkSolidConeDR04                 = it->nTrkSolidConeDR04();
      pho.nTrkHollowConeDR04                = it->nTrkHollowConeDR04();

      pho.ecalRecHitSumEtConeDR03           = it->ecalRecHitSumEtConeDR03();
      pho.hcalTowerSumEtConeDR03            = it->hcalTowerSumEtConeDR03();
      pho.hcalDepth1TowerSumEtConeDR03      = it->hcalDepth1TowerSumEtConeDR03();
      pho.hcalDepth2TowerSumEtConeDR03      = it->hcalDepth2TowerSumEtConeDR03();
      pho.trkSumPtSolidConeDR03             = it->trkSumPtSolidConeDR03();
      pho.trkSumPtHollowConeDR03            = it->trkSumPtHollowConeDR03();
      pho.nTrkSolidConeDR03                 = it->nTrkSolidConeDR03();
      pho.nTrkHollowConeDR03                = it->nTrkHollowConeDR03();

      pho.chargedHadronIso                  = it->chargedHadronIso();
      pho.neutralHadronIso                  = it->neutralHadronIso();
      pho.photonIso                         = it->photonIso();

      if(it->superCluster().isNonnull()) {
	if(recoMode_) {
	  double e1000 = spr::eECALmatrix(it->superCluster()->seed()->seed(),barrelRecHitsHandle,endcapRecHitsHandle,geo,caloTopology,1,0,0,0);
	  double e0100 = spr::eECALmatrix(it->superCluster()->seed()->seed(),barrelRecHitsHandle,endcapRecHitsHandle,geo,caloTopology,0,1,0,0);
	  double e0010 = spr::eECALmatrix(it->superCluster()->seed()->seed(),barrelRecHitsHandle,endcapRecHitsHandle,geo,caloTopology,0,0,1,0);
	  double e0001 = spr::eECALmatrix(it->superCluster()->seed()->seed(),barrelRecHitsHandle,endcapRecHitsHandle,geo,caloTopology,0,0,0,1);
	  pho.e1x2 = std::max( std::max(e1000,e0100), std::max(e0010, e0001) );
	}

	pho.superClusterPreshowerEnergy = it->superCluster()->preshowerEnergy();
	pho.superClusterPhiWidth = it->superCluster()->phiWidth();
	pho.superClusterEtaWidth = it->superCluster()->etaWidth();
	susy::SuperCluster superCluster;
	fillCluster(it->superCluster(), superCluster, clusterIndex);
	superCluster.index = superClusterIndex;
	pho.superClusterIndex = superClusterIndex;
	susyEvent_->superClusters.push_back(superCluster);
	superClusterIndex++;
      }

      pho.caloPosition.SetXYZ(it->caloPosition().x(),it->caloPosition().y(),it->caloPosition().z());
      pho.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());

      for(int k=0; k<nPhoIDC; k++){
	pho.idPairs[ TString(photonIDCollectionTags_[k].c_str()) ] = (*phoIds[k])[phoRef];
      }// for id

      // conversion ID
      pho.dist   = 0;
      pho.dcot   = 0;
      pho.radius = 0;

      susyEvent_->photons.push_back(pho);

      if(debugLevel_ > 2) std::cout << "pt, e, hadEm : " << it->pt()
				     << ", " << it->energy()
				     << ", " << it->hadronicOverEm() << std::endl;
    }// for it
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "reco::Photon is not available!!! " << e.what();
  }



  if(debugLevel_ > 0) std::cout << name() << ", fill electron" << std::endl;

  const int nEleIDC = electronIDCollectionTags_.size();
  std::vector< const edm::ValueMap<Float_t>* > eleIds;
  
  for(int i=0; i<nEleIDC; i++) {
    edm::Handle<edm::ValueMap<Float_t> > eleIDCH;
    try {
      iEvent.getByLabel(electronIDCollectionTags_[i], eleIDCH);
      eleIds.push_back( eleIDCH.product() );
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << electronIDCollectionTags_[i] << " is not available!!! " << e.what();
    }
  }

  edm::Handle<reco::GsfElectronCollection > electronH;
  try {
    iEvent.getByLabel(electronCollectionTag_,electronH);
    if(debugLevel_ > 1) std::cout << "size of ElectronCollection : " << electronH->size() << std::endl;
    int iele = 0;
    for(reco::GsfElectronCollection::const_iterator it = electronH->begin();
	it != electronH->end(); it++){

      reco::GsfElectronRef eleRef(electronH,iele++);

      susy::Electron ele;
      // fiducial bits
      ele.fidBit  = 0;
      ele.fidBit |= it->isEB()        << 0;
      ele.fidBit |= it->isEE()        << 1;
      ele.fidBit |= it->isEBEEGap()   << 2;
      ele.fidBit |= it->isEBEtaGap()  << 3;
      ele.fidBit |= it->isEBPhiGap()  << 4;
      ele.fidBit |= it->isEEDeeGap()  << 5;
      ele.fidBit |= it->isEERingGap() << 6;

      ele.scPixCharge = it->scPixCharge();
      ele.boolPack = 0;
      ele.boolPack |= it->isGsfCtfScPixChargeConsistent() << 0;
      ele.boolPack |= it->isGsfScPixChargeConsistent()    << 1;
      ele.boolPack |= it->isGsfCtfChargeConsistent()      << 2;
      ele.boolPack |= it->ecalDrivenSeed()                << 3;
      ele.boolPack |= it->trackerDrivenSeed()             << 4;
      ele.boolPack |= it->ecalDriven()                    << 5;
      ele.boolPack |= it->passingCutBasedPreselection()   << 6;
      ele.boolPack |= it->passingMvaPreselection()        << 7;
      ele.boolPack |= it->ambiguous()                     << 8;
      ele.boolPack |= it->isEcalEnergyCorrected()         << 9;
      ele.boolPack |= it->isMomentumCorrected()           << 10;
      ele.boolPack |= it->convFlags()                     << 11;

      ele.eSuperClusterOverP             = it->eSuperClusterOverP();
      ele.eSeedClusterOverP              = it->eSeedClusterOverP();
      ele.eSeedClusterOverPout           = it->eSeedClusterOverPout();
      ele.eEleClusterOverPout            = it->eEleClusterOverPout();
      ele.deltaEtaSuperClusterTrackAtVtx = it->deltaEtaSuperClusterTrackAtVtx();
      ele.deltaEtaSeedClusterTrackAtCalo = it->deltaEtaSeedClusterTrackAtCalo();
      ele.deltaEtaEleClusterTrackAtCalo  = it->deltaEtaEleClusterTrackAtCalo();
      ele.deltaPhiSuperClusterTrackAtVtx = it->deltaPhiSuperClusterTrackAtVtx();
      ele.deltaPhiSeedClusterTrackAtCalo = it->deltaPhiSeedClusterTrackAtCalo();
      ele.deltaPhiEleClusterTrackAtCalo  = it->deltaPhiEleClusterTrackAtCalo();

      ele.trackPositions["AtVtx"] = TVector3(it->trackPositionAtVtx().X(),it->trackPositionAtVtx().Y(),it->trackPositionAtVtx().Z());
      ele.trackPositions["AtCalo"] = TVector3(it->trackPositionAtCalo().X(),it->trackPositionAtCalo().Y(),it->trackPositionAtCalo().Z());
      ele.trackMomentums["AtVtx"] = TLorentzVector(it->trackMomentumAtVtx().X(),it->trackMomentumAtVtx().Y(),
						   it->trackMomentumAtVtx().Z(),it->trackMomentumAtVtx().R());
      ele.trackMomentums["AtCalo"] = TLorentzVector(it->trackMomentumAtCalo().X(),it->trackMomentumAtCalo().Y(),
						    it->trackMomentumAtCalo().Z(),it->trackMomentumAtCalo().R());
      ele.trackMomentums["Out"] = TLorentzVector(it->trackMomentumOut().X(),it->trackMomentumOut().Y(),
						 it->trackMomentumOut().Z(),it->trackMomentumOut().R());
      ele.trackMomentums["AtEleClus"] = TLorentzVector(it->trackMomentumAtEleClus().X(),it->trackMomentumAtEleClus().Y(),
						       it->trackMomentumAtEleClus().Z(),it->trackMomentumAtEleClus().R());
      ele.trackMomentums["AtVtxWithConstraint"] = TLorentzVector(it->trackMomentumAtVtxWithConstraint().X(),it->trackMomentumAtVtxWithConstraint().Y(),
								 it->trackMomentumAtVtxWithConstraint().Z(),it->trackMomentumAtVtxWithConstraint().R());

      ele.vertex.SetXYZ(it->vx(),it->vy(),it->vz());
      ele.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());

      ele.shFracInnerHits = it->shFracInnerHits();

      ele.sigmaEtaEta                       = it->sigmaEtaEta();
      ele.sigmaIetaIeta                     = it->sigmaIetaIeta();
      ele.e1x5                              = it->e1x5();
      ele.e2x5Max                           = it->e2x5Max();
      ele.e5x5                              = it->e5x5();
      ele.hcalDepth1OverEcal                = it->hcalDepth1OverEcal();
      ele.hcalDepth2OverEcal                = it->hcalDepth2OverEcal();
      ele.hcalOverEcal                      = it->hcalOverEcal();

      ele.dr03TkSumPt              = it->dr03TkSumPt();
      ele.dr03EcalRecHitSumEt      = it->dr03EcalRecHitSumEt();
      ele.dr03HcalDepth1TowerSumEt = it->dr03HcalDepth1TowerSumEt();
      ele.dr03HcalDepth2TowerSumEt = it->dr03HcalDepth2TowerSumEt();
      ele.dr03HcalTowerSumEt       = it->dr03HcalTowerSumEt();
 
      ele.dr04TkSumPt              = it->dr04TkSumPt();
      ele.dr04EcalRecHitSumEt      = it->dr04EcalRecHitSumEt();
      ele.dr04HcalDepth1TowerSumEt = it->dr04HcalDepth1TowerSumEt();
      ele.dr04HcalDepth2TowerSumEt = it->dr04HcalDepth2TowerSumEt();
      ele.dr04HcalTowerSumEt       = it->dr04HcalTowerSumEt();
 
      ele.convDist   = it->convDist();
      ele.convDcot   = it->convDcot();
      ele.convRadius = it->convRadius();

      ele.mva        = it->mva();

      //enum Classification { UNKNOWN=-1, GOLDEN=0, BIGBREM=1, OLDNARROW=2, SHOWERING=3, GAP=4 } ;
      ele.bremClass                         = int(it->classification());
      ele.fbrem                             = it->fbrem();

      ele.ecalEnergy                        = it->ecalEnergy();
      ele.ecalEnergyError                   = it->ecalEnergyError();
      ele.trackMomentumError                = it->trackMomentumError();
      ele.electronMomentumError             = it->electronMomentumError();

      susy::Cluster cluster;
      fillCluster(it->electronCluster(),cluster);
      cluster.index = clusterIndex;
      ele.electronClusterIndex = clusterIndex;
      susyEvent_->clusters.push_back(cluster);
      clusterIndex++;

      susy::SuperCluster superCluster;
      fillCluster(it->superCluster(), superCluster, clusterIndex);
      superCluster.index = superClusterIndex;
      ele.superClusterIndex = superClusterIndex;
      susyEvent_->superClusters.push_back(superCluster);
      superClusterIndex++;

      if(it->closestCtfTrackRef().isNonnull()){
	susy::Track track;
	fillTrack(it->closestCtfTrackRef(), track);
	track.index = trackIndex;
	ele.closestCtfTrackIndex = trackIndex;
	susyEvent_->tracks.push_back(track);
	trackIndex++;
      }
      if(it->gsfTrack().isNonnull()) {
	susy::Track track;
	fillTrack(it->gsfTrack(), track);
	track.index = trackIndex;
	track.extrapolatedPositions["ECALInnerWall"] = TVector3(it->trackPositionAtCalo().X(),it->trackPositionAtCalo().Y(),it->trackPositionAtCalo().Z());
	ele.gsfTrackIndex = trackIndex;
	susyEvent_->tracks.push_back(track);
	trackIndex++;
      }// if(it->gsfTrack().isNonnull())


      for(int k=0; k<nEleIDC; k++){
	ele.idPairs[ TString(electronIDCollectionTags_[k].c_str()) ] = (*eleIds[k])[eleRef];
      }// for id

      susyEvent_->electrons.push_back(ele);

      if(debugLevel_ > 2) std::cout << "pt, e, hadEm : " << it->pt()
				     << ", " << it->energy()
				     << ", " << it->hadronicOverEm() << std::endl;
    }// for it
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "reco::Electron is not available!!! " << e.what();
  }




  if(debugLevel_ > 0) std::cout << name() << ", fill muon" << std::endl;

  const int nMuIDC = muonIDCollectionTags_.size();
  std::vector< const edm::ValueMap<Bool_t>* > muIds;
  
  for(int i=0; i<nMuIDC; i++) {
    edm::Handle<edm::ValueMap<Bool_t> > muIDCH;
    try {
      iEvent.getByLabel(muonIDCollectionTags_[i], muIDCH);
      muIds.push_back( muIDCH.product() );
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << muonIDCollectionTags_[i] << " is not available!!! " << e.what();
    }
  }

  edm::Handle<reco::MuonCollection > muonH;
  try {
    iEvent.getByLabel(muonCollectionTag_,muonH);
    if(debugLevel_ > 1) std::cout << "size of MuonCollection : " << muonH->size() << std::endl;
    int imu = 0;
    for(reco::MuonCollection::const_iterator it = muonH->begin();
	it != muonH->end(); it++){

      reco::MuonRef muRef(muonH,imu++);

      susy::Muon mu;
      mu.type                               = it->type();

      mu.caloCompatibility                  = it->caloCompatibility();
      mu.nMatches                           = it->numberOfMatches();
      mu.nChambers                          = it->numberOfChambers();
      if(it->combinedMuon().isNonnull()){
	mu.nValidHits                         = it->combinedMuon()->hitPattern().numberOfValidHits();
	mu.nValidTrackerHits                  = it->combinedMuon()->hitPattern().numberOfValidTrackerHits();
	mu.nValidMuonHits                     = it->combinedMuon()->hitPattern().numberOfValidMuonHits();
      }
      mu.emEnergy                           = it->calEnergy().em;
      mu.hadEnergy                          = it->calEnergy().had;
      mu.trackIsoR03                        = it->isolationR03().sumPt;
      mu.ecalIsoR03                         = it->isolationR03().emEt;
      mu.hcalIsoR03                         = it->isolationR03().hadEt;
      mu.trackIsoR05                        = it->isolationR05().sumPt;
      mu.ecalIsoR05                         = it->isolationR05().emEt;
      mu.hcalIsoR05                         = it->isolationR05().hadEt;

      if(it->isTimeValid()){
	mu.timeNDof                         = it->time().nDof;
	mu.timeDirection                    = Int_t(it->time().direction());
	if(it->time().direction() == reco::MuonTime::OutsideIn){
	  mu.timeAtIp                       = it->time().timeAtIpOutIn;
	  mu.timeAtIpError                  = it->time().timeAtIpOutInErr;
	}
	else if(it->time().direction() == reco::MuonTime::InsideOut){
	  mu.timeAtIp                       = it->time().timeAtIpInOut;
	  mu.timeAtIpError                  = it->time().timeAtIpInOutErr;
	}
      }// if(it->isTimeValid())

      if(it->track().isNonnull()){
	susy::Track track;
	fillTrack(it->track(), track);
	track.index = trackIndex;
	mu.trackIndex = trackIndex;
	susyEvent_->tracks.push_back(track);
	trackIndex++;
      }
      if(it->standAloneMuon().isNonnull()){
	susy::Track track;
	fillTrack(it->standAloneMuon(), track);
	track.index = trackIndex;
	mu.standAloneTrackIndex = trackIndex;
	susyEvent_->tracks.push_back(track);
	trackIndex++;
      }
      if(it->combinedMuon().isNonnull()){
	susy::Track track;
	fillTrack(it->combinedMuon(), track);
	track.index = trackIndex;
	mu.combinedTrackIndex = trackIndex;
	susyEvent_->tracks.push_back(track);
	trackIndex++;
      }

      mu.momentum.SetXYZT(it->p4().px(),it->p4().py(),it->p4().pz(),it->p4().e());

      for(int k=0; k<nMuIDC; k++){
	mu.idPairs[ TString(muonIDCollectionTags_[k].c_str()) ] = (*muIds[k])[muRef];
      }// for id

      susyEvent_->muons.push_back(mu);
      if(debugLevel_ > 2) std::cout << "type, emE, hadE, pt : " << it->type()
				    << ", " << it->calEnergy().em
				    << ", " << it->calEnergy().had
				    << ", " << it->pt() << std::endl;
    }
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "reco::Muon is not available!!! " << e.what();
  }



  if(debugLevel_ > 0) std::cout << name() << ", fill calojet collections" << std::endl;

  int nJetColl = caloJetCollectionTags_.size();

  for(int iJetC=0; iJetC < nJetColl; iJetC++) {
    susy::CaloJetCollection jetCollection;
    TString key = TString(caloJetCollectionTags_[iJetC].c_str()).ReplaceAll("CaloJets","");

    edm::Handle<edm::ValueMap<reco::JetID> > jetIDH;
    try {
      iEvent.getByLabel((key+"JetID").Data(), jetIDH);
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << caloJetCollectionTags_[iJetC] << " of JetID is not available!!! " << e.what();
    }

    edm::Handle<reco::CaloJetCollection> jetH;
    try {
      iEvent.getByLabel(edm::InputTag(caloJetCollectionTags_[iJetC]),jetH);
      if(debugLevel_ > 1) std::cout << "size of " << key << " JetCollection : " << jetH->size() << std::endl;
      int ijet = 0;
      for(reco::CaloJetCollection::const_iterator it = jetH->begin();
	  it != jetH->end(); it++){

	reco::CaloJetRef jetRef(jetH,ijet++);

	if(it->pt() < 10.0) continue;

	susy::CaloJet jet;

	// Basic Jet
	//	jet.partonFlavour = it->partonFlavour();
	//	jet.charge = it->jetCharge();
	jet.etaMean         = it->etaPhiStatistics().etaMean;
	jet.phiMean         = it->etaPhiStatistics().phiMean;
	jet.etaEtaMoment    = it->etaPhiStatistics().etaEtaMoment;
	jet.etaPhiMoment    = it->etaPhiStatistics().etaPhiMoment;
	jet.phiPhiMoment    = it->etaPhiStatistics().phiPhiMoment;
	jet.maxDistance     = it->maxDistance();
	jet.jetArea         = it->jetArea();
	jet.pileup          = it->pileup();
	jet.nPasses         = it->nPasses();
	jet.nConstituents   = it->nConstituents();

	// CaloJet specific
	jet.maxEInEmTowers                = it->maxEInEmTowers();
	jet.maxEInHadTowers               = it->maxEInHadTowers();
	jet.energyFractionHadronic        = it->energyFractionHadronic();
	jet.emEnergyFraction              = it->emEnergyFraction();
	jet.hadEnergyInHB                 = it->hadEnergyInHB();
	jet.hadEnergyInHO                 = it->hadEnergyInHO();
	jet.hadEnergyInHE                 = it->hadEnergyInHE();
	jet.emEnergyInHF                  = it->emEnergyInHF();
	jet.emEnergyInEB                  = it->emEnergyInEB();
	jet.emEnergyInEE                  = it->emEnergyInEE();
	jet.hadEnergyInHF                 = it->hadEnergyInHF();
	jet.towersArea                    = it->towersArea();
	jet.n90                           = it->n90();
	jet.n60                           = it->n60();
	
	jet.vertex.SetXYZ(it->vx(),it->vy(),it->vz());
	jet.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());
	jet.detectorP4.SetXYZT(it->detectorP4().px(),it->detectorP4().py(),
			       it->detectorP4().pz(),it->detectorP4().energy());

	// accessing Jet ID information
	const reco::JetID& jetID = (*jetIDH)[jetRef];
	jet.fHPD                          = jetID.fHPD;
	jet.fRBX                          = jetID.fRBX;
	jet.n90Hits                       = jetID.n90Hits;
	jet.fSubDetector1                 = jetID.fSubDetector1;
	jet.fSubDetector2                 = jetID.fSubDetector2;
	jet.fSubDetector3                 = jetID.fSubDetector3;
	jet.fSubDetector4                 = jetID.fSubDetector4;
	jet.restrictedEMF                 = jetID.restrictedEMF;
	jet.nHCALTowers                   = jetID.nHCALTowers;
	jet.nECALTowers                   = jetID.nECALTowers;
	jet.approximatefHPD               = jetID.approximatefHPD;
	jet.approximatefRBX               = jetID.approximatefRBX;
	jet.hitsInN90                     = jetID.hitsInN90;
	jet.numberOfHits2RPC              = jetID.numberOfHits2RPC;
	jet.numberOfHits3RPC              = jetID.numberOfHits3RPC;
	jet.numberOfHitsRPC               = jetID.numberOfHitsRPC;
	
	jetCollection.push_back(jet);
	if(debugLevel_ > 2) std::cout << "pt, e : " << it->pt() << ", " << it->energy() << std::endl;

      }// for it
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "reco::CaloJet is not available!!! " << e.what();
    }

    susyEvent_->caloJets[key] = jetCollection;

  }// for iJetC



  if(debugLevel_ > 0) std::cout << name() << ", fill pfjet collections" << std::endl;

  nJetColl = pfJetCollectionTags_.size();
  for(int iJetC=0; iJetC < nJetColl; iJetC++) {
    susy::PFJetCollection jetCollection;
    TString key = TString(pfJetCollectionTags_[iJetC].c_str()).ReplaceAll("PFJets","");

    edm::Handle<edm::ValueMap<reco::JetID> > jetIDH;
    try {
      iEvent.getByLabel((key+"JetID").Data(), jetIDH);
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << pfJetCollectionTags_[iJetC] << " of JetID is not available!!! " << e.what();
    }
    edm::Handle<reco::PFJetCollection> jetH;
    try {
      iEvent.getByLabel(edm::InputTag(pfJetCollectionTags_[iJetC]),jetH);
      if(debugLevel_ > 1) std::cout << "size of " << key << " JetCollection : " << jetH->size() << std::endl;
      int ijet = 0;
      for(reco::PFJetCollection::const_iterator it = jetH->begin();
	  it != jetH->end(); it++){

	reco::PFJetRef jetRef(jetH,ijet++);

	if(it->pt() < 10.0) continue;

	susy::PFJet jet;

	// Basic Jet
	//	jet.partonFlavour = it->partonFlavour();
	//	jet.charge = it->jetCharge();
	jet.etaMean         = it->etaPhiStatistics().etaMean;
	jet.phiMean         = it->etaPhiStatistics().phiMean;
	jet.etaEtaMoment    = it->etaPhiStatistics().etaEtaMoment;
	jet.etaPhiMoment    = it->etaPhiStatistics().etaPhiMoment;
	jet.phiPhiMoment    = it->etaPhiStatistics().phiPhiMoment;
	jet.maxDistance     = it->maxDistance();
	jet.jetArea         = it->jetArea();
	jet.pileup          = it->pileup();
	jet.nPasses         = it->nPasses();
	jet.nConstituents   = it->nConstituents();

	
	jet.vertex.SetXYZ(it->vx(),it->vy(),it->vz());
	jet.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());

	jet.chargedHadronEnergy = it->chargedHadronEnergy();
	jet.neutralHadronEnergy = it->neutralHadronEnergy();
	jet.photonEnergy = it->photonEnergy();
	jet.electronEnergy = it->electronEnergy();
	jet.muonEnergy = it->muonEnergy();
	jet.HFHadronEnergy = it->HFHadronEnergy();
	jet.HFEMEnergy = it->HFEMEnergy();
	jet.chargedHadronMultiplicity = it->chargedHadronMultiplicity();
	jet.neutralHadronMultiplicity = it->neutralHadronMultiplicity();
	jet.photonMultiplicity = it->photonMultiplicity();
	jet.electronMultiplicity = it->electronMultiplicity();
	jet.muonMultiplicity = it->muonMultiplicity();
	jet.HFHadronMultiplicity = it->HFHadronMultiplicity();
	jet.HFEMMultiplicity = it->HFEMMultiplicity();
	jet.chargedEmEnergy = it->chargedEmEnergy();
	jet.chargedMuEnergy = it->chargedMuEnergy();
	jet.neutralEmEnergy = it->neutralEmEnergy();
	jet.chargedMultiplicity = it->chargedMultiplicity();
	jet.neutralMultiplicity = it->neutralMultiplicity();

	// accessing Jet ID information
	const reco::JetID& jetID = (*jetIDH)[jetRef];
	jet.fHPD                          = jetID.fHPD;
	jet.fRBX                          = jetID.fRBX;
	jet.n90Hits                       = jetID.n90Hits;
	jet.fSubDetector1                 = jetID.fSubDetector1;
	jet.fSubDetector2                 = jetID.fSubDetector2;
	jet.fSubDetector3                 = jetID.fSubDetector3;
	jet.fSubDetector4                 = jetID.fSubDetector4;
	jet.restrictedEMF                 = jetID.restrictedEMF;
	jet.nHCALTowers                   = jetID.nHCALTowers;
	jet.nECALTowers                   = jetID.nECALTowers;
	jet.approximatefHPD               = jetID.approximatefHPD;
	jet.approximatefRBX               = jetID.approximatefRBX;
	jet.hitsInN90                     = jetID.hitsInN90;
	jet.numberOfHits2RPC              = jetID.numberOfHits2RPC;
	jet.numberOfHits3RPC              = jetID.numberOfHits3RPC;
	jet.numberOfHitsRPC               = jetID.numberOfHitsRPC;

	jetCollection.push_back(jet);
	if(debugLevel_ > 2) std::cout << "pt, e : " << it->pt() << ", " << it->energy() << std::endl;

      }// for it
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "reco::Jet is not available!!! " << e.what();
    }

    susyEvent_->pfJets[key] = jetCollection;

  }// for PFJet



  if(debugLevel_ > 0) std::cout << name() << ", fill jptjet collections" << std::endl;

  nJetColl = jptJetCollectionTags_.size();
  for(int iJetC=0; iJetC < nJetColl; iJetC++) {
    susy::JPTJetCollection jetCollection;
    TString key(jptJetCollectionTags_[iJetC].c_str());

    edm::Handle<reco::JPTJetCollection> jetH;
    try {
      iEvent.getByLabel(edm::InputTag(jptJetCollectionTags_[iJetC]),jetH);
      if(debugLevel_ > 1) std::cout << "size of " << key << " JetCollection : " << jetH->size() << std::endl;
      int ijet = 0;
      for(reco::JPTJetCollection::const_iterator it = jetH->begin();
	  it != jetH->end(); it++){

	reco::JPTJetRef jetRef(jetH,ijet++);

	if(it->pt() < 10.0) continue;

	susy::JPTJet jet;

	// Basic Jet
	//	jet.partonFlavour = it->partonFlavour();
	//	jet.charge = it->jetCharge();
	jet.etaMean         = it->etaPhiStatistics().etaMean;
	jet.phiMean         = it->etaPhiStatistics().phiMean;
	jet.etaEtaMoment    = it->etaPhiStatistics().etaEtaMoment;
	jet.etaPhiMoment    = it->etaPhiStatistics().etaPhiMoment;
	jet.phiPhiMoment    = it->etaPhiStatistics().phiPhiMoment;
	jet.maxDistance     = it->maxDistance();
	jet.jetArea         = it->jetArea();
	jet.pileup          = it->pileup();
	jet.nPasses         = it->nPasses();
	jet.nConstituents   = it->nConstituents();

	jet.vertex.SetXYZ(it->vx(),it->vy(),it->vz());
	jet.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());

	jet.chargedHadronEnergy = it->chargedHadronEnergy();
	jet.neutralHadronEnergy = it->neutralHadronEnergy();
	jet.chargedEmEnergy     = it->chargedEmEnergy();
	jet.neutralEmEnergy     = it->neutralEmEnergy();
	jet.chargedMultiplicity = it->chargedMultiplicity();
	jet.muonMultiplicity    = it->muonMultiplicity();
	jet.elecMultiplicity    = it->elecMultiplicity();
	jet.getZSPCor           = it->getZSPCor();

	jetCollection.push_back(jet);
	if(debugLevel_ > 2) std::cout << "pt, e : " << it->pt() << ", " << it->energy() << std::endl;
      }// for it
    }
    catch(cms::Exception& e) {
      edm::LogError(name()) << "reco::Jet is not available!!! " << e.what();
    }

    susyEvent_->jptJets[key] = jetCollection;

  }// for JPTJet




  if( storeGenInfos_ && ! iEvent.isRealData() ) {
    if(debugLevel_ > 0) std::cout << name() << ", fill generated informations" << std::endl;
    fillGenInfos(iEvent, iSetup);
  }



  if(debugLevel_ > 0) std::cout << name() << ", fill tree at last" << std::endl;

  susyTree_->Fill();


  // end of event cleaning procedure - delete any variables created by "new" operator in this procedure
  //  if(recoMode_ && propagator_) delete propagator_;

}


void SusyNtuplizer::fillTrack(const reco::TrackRef& in, susy::Track& out) {

  if(in.isNull()) return;

  try {
    out.algorithm = Int_t(in->algo());
    out.quality = in->qualityMask();

    out.chi2 = in->chi2();
    out.ndof = in->ndof();
    out.charge = in->charge();
    for(int i=0; i<reco::Track::dimension; i++) out.error[i] = in->error(i);
    out.vertex.SetXYZ(in->vx(),in->vy(),in->vz());
    out.momentum.SetXYZT(in->px(),in->py(),in->pz(),in->p());

    out.statusCode = 0x0; // normal track flag

    if(recoMode_){
      out.statusCode |= ((in->innerOk() ? 0x1 : 0) << 1);
      out.statusCode |= ((in->outerOk() ? 0x1 : 0) << 2);

      out.nHits = in->numberOfValidHits();
      out.hitPositions["inner"] = TVector3(in->innerPosition().x(),in->innerPosition().y(),in->innerPosition().z());
      out.hitPositions["outer"] = TVector3(in->outerPosition().x(),in->outerPosition().y(),in->outerPosition().z());
      out.extrapolatedPositions.clear();
      fillExtrapolations(in,out.extrapolatedPositions);
    }// if(recoMode_)

  }// try
  catch(cms::Exception& e) {
    edm::LogError(name()) << " TrackRef is not available!!! " << e.what();
  }
 
}


void SusyNtuplizer::fillTrack(const reco::GsfTrackRef& in, susy::Track& out) {

  if(in.isNull()) return;

  try {
    out.algorithm = Int_t(in->algo());
    out.quality = in->qualityMask();

    out.chi2 = in->chi2();
    out.ndof = in->ndof();
    out.charge = in->chargeMode();
    for(int i=0; i<reco::GsfTrack::dimensionMode; i++) out.error[i] = in->errorMode(i);
    out.vertex.SetXYZ(in->vx(),in->vy(),in->vz());
    out.momentum.SetXYZT(in->pxMode(),in->pyMode(),in->pzMode(),in->pMode());

    out.statusCode = 0x1; // gsfTrack flag

    if(recoMode_){
      out.statusCode |= ((in->innerOk() ? 0x1 : 0) << 1);
      out.statusCode |= ((in->outerOk() ? 0x1 : 0) << 2);

      out.nHits = in->numberOfValidHits();
      out.hitPositions["inner"] = TVector3(in->innerPosition().x(),in->innerPosition().y(),in->innerPosition().z());
      out.hitPositions["outer"] = TVector3(in->outerPosition().x(),in->outerPosition().y(),in->outerPosition().z());
      out.extrapolatedPositions.clear();
    }// if(recoMode_)

  }// try
  catch(cms::Exception& e) {
    edm::LogError(name()) << " GsfTrackRef is not available!!! " << e.what();
  }

}


void SusyNtuplizer::fillCluster(const reco::CaloClusterPtr& in, susy::Cluster& out) {

  if(in.isNull()) return;

  try {
    out.energy = in->energy();
    out.position.SetXYZ(in->x(),in->y(),in->z());
    out.nCrystals = in->size();
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << " CaloClusterPtr is not available!!! " << e.what();
  }

}


void SusyNtuplizer::fillCluster(const reco::SuperClusterRef& in, susy::SuperCluster& out, int& basicClusterIndex) {

  if(in.isNull()) return;

  try {
    out.energy = in->energy();
    out.preshowerEnergy = in->preshowerEnergy();
    out.phiWidth = in->phiWidth();
    out.etaWidth = in->etaWidth();
    out.position.SetXYZ(in->x(),in->y(),in->z());

    for(reco::CaloClusterPtrVector::const_iterator it = in->clustersBegin();
	it != in->clustersEnd(); it++){
      susy::Cluster cluster;
      fillCluster(*it,cluster);
      if(in->seed() == *it) out.seedClusterIndex = basicClusterIndex;
      cluster.index = basicClusterIndex;
      out.basicClusterIndices.push_back(basicClusterIndex);
      susyEvent_->clusters.push_back(cluster);
      basicClusterIndex++;
    }
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << " SuperClusterRef is not available!!! " << e.what();
  }
}


void SusyNtuplizer::fillParticle(const reco::PFCandidateRef& in, susy::Particle& out) {

  if(in.isNull()) return;

  try {
    out.status = in->status();
    out.pdgId = in->pdgId();
    out.charge = in->charge();
    out.vertex.SetXYZ(in->vx(),in->vy(),in->vz());
    out.momentum.SetXYZT(in->px(),in->py(),in->pz(),in->energy());
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << " PFCandidateRef is not available!!! " << e.what();
  }

}


void SusyNtuplizer::fillParticle(const reco::GenParticle* in, susy::Particle& out, int igen) {

  if(in == 0) return;

  try {
    out.index = igen;
    out.status = in->status();
    out.pdgId = in->pdgId();
    out.charge = in->charge();
    out.vertex.SetXYZ(in->vx(),in->vy(),in->vz());
    out.momentum.SetXYZT(in->px(),in->py(),in->pz(),in->energy());
  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << " PFCandidateRef is not available!!! " << e.what();
  }

}


void SusyNtuplizer::fillExtrapolations(const reco::TrackRef& rtrk, std::map<TString,TVector3>& positions) {

  if(rtrk.isNull()) return;

  try {
    reco::TransientTrack  ttk = transientTrackBuilder_->build(*rtrk);
    const TrajectoryStateOnSurface tsos = ttk.innermostMeasurementState();
    if(!tsos.isValid()) return;
    float eta = tsos.globalPosition().eta();
    TrajectoryStateOnSurface state;
    TVector3 position(0.0,0.0,0.0);

    if(std::abs(eta) < susy::etaGap){
      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.barrelBound(PFGeometry::ECALInnerWall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["ECALInnerWall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.barrelBound(PFGeometry::HCALInnerWall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["HCALInnerWall"] = position;
    }
    else if(eta > susy::etaGap && eta < susy::etaMax) {
      state = propagator_->propagate(tsos, pfGeom_.positiveEndcapDisk(PFGeometry::ECALInnerWall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["ECALInnerWall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.positiveEndcapDisk(PFGeometry::HCALInnerWall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["HCALInnerWall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.positiveEndcapDisk(PFGeometry::PS1Wall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["PS1Wall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.positiveEndcapDisk(PFGeometry::PS2Wall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["PS2Wall"] = position;
    }
    else if(eta < -susy::etaGap && eta > -susy::etaMax) {
      state = propagator_->propagate(tsos, pfGeom_.negativeEndcapDisk(PFGeometry::ECALInnerWall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["ECALInnerWall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.negativeEndcapDisk(PFGeometry::HCALInnerWall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["HCALInnerWall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.negativeEndcapDisk(PFGeometry::PS1Wall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["PS1Wall"] = position;

      position *= 0;
      state = propagator_->propagate(tsos, pfGeom_.negativeEndcapDisk(PFGeometry::PS2Wall));
      if(state.isValid()) position.SetXYZ(state.globalPosition().x(),state.globalPosition().y(),state.globalPosition().z());
      positions["PS2Wall"] = position;
    }

  }
  catch(cms::Exception& e) {
    edm::LogError(name()) << "Something wrong in fillExtrapolations!!! " << e.what();
  }

}



void SusyNtuplizer::fillGenInfos(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // GenParticle : "genParticles"
  // GenParticles will be stored with status != 3 && pt > 5 GeV only

  edm::Handle<reco::GenParticleCollection> gpH;
  try {
    iEvent.getByLabel(genCollectionTag_,gpH);
    reco::GenParticleCollection::const_iterator it_begin = gpH->begin();
    reco::GenParticleCollection::const_iterator it_end   = gpH->end();

    int igen = 0;
    for(reco::GenParticleCollection::const_iterator it = it_begin; it != it_end; it++, igen++){
      //      if(it->status() == 3) continue;
      //      if(it->pt() < 5.0) continue;
      const reco::GenParticle* gp = &*it;
      susy::Particle part;
      fillParticle(gp,part,igen);
      susyEvent_->genParticles.push_back(part);
    }// for

    // fill motherIndex in genParticles
    igen = 0;
    for(reco::GenParticleCollection::const_iterator it = it_begin; it != it_end; it++, igen++){
      const reco::GenParticle* gp = &*it;
      susy::Particle* part = &(susyEvent_->genParticles[igen]);
      //      const Candidate * mom = gp->mother();
      const reco::GenParticle* mom = (const reco::GenParticle*) gp->mother();
      if(mom == 0) {
	part->motherIndex = -1;
	continue;
      }
      int jgen = 0;
      for(reco::GenParticleCollection::const_iterator it2 = it_begin; it2 != it_end; it2++, jgen++){
	const reco::GenParticle* gp2 = &*it2;
	if(sameGenParticles(mom,gp2)) {
	  part->motherIndex = jgen;
	  break;
	}
      }
    }// for


  }// try
  catch (cms::Exception& e) {
    edm::LogError(name()) << "reco::GenParticle is not available!!! " << e.what();

  }

  // SimVertex info is not available in SUSYPAT sample.

}


bool SusyNtuplizer::sameGenParticles(const reco::GenParticle* gp, const reco::GenParticle* gp2) {

  if( gp->pdgId() != gp2->pdgId() ) return false;
  if( gp->status() != gp2->status() ) return false;
  if( gp->charge() != gp2->charge() ) return false;
  TVector3 vtx(gp->vx() - gp2->vx(), gp->vy() - gp2->vy(), gp->vz() - gp2->vz());
  if(vtx.Mag() > 1.0e-6) return false;
  TVector3 p(gp->px() - gp2->px(),gp->py() - gp2->py(),gp->pz() - gp2->pz());
  if(p.Mag() > 1.0e-6) return false;
  if(gp->energy() - gp2->energy() > 1.0e-6) return false;

  return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(SusyNtuplizer);
