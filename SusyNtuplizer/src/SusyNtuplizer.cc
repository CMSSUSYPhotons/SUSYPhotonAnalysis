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
// $Id: SusyNtuplizer.cc,v 1.49 2013/04/01 09:53:25 yiiyama Exp $
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
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

// simple geometry
#include "RecoParticleFlow/PFTracking/interface/PFGeometry.h"

// for track extrapolation
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

//for conversion safe electron veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

// for ecal rechit related
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Calibration/IsolatedParticles/interface/eECALMatrix.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

// Jet Energy Correction
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

// pileup summary info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// b-tagging info
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

// Pileup Jet Id
#include "CMGTools/External/interface/PileupJetIdentifier.h"

// PFIsolation
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"

//Photon SC energy MVA regression
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"

//HCAL laser events 2012 filter
#include "PhysicsTools/Utilities/interface/EventFilterFromListStandAlone.h"

// system include files
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <utility>
#include <set>
#include <list>
#include <map>
#include <cstdlib>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

// Class definition
class SusyNtuplizer : public edm::EDAnalyzer {

public:
  explicit SusyNtuplizer(edm::ParameterSet const&);
  ~SusyNtuplizer();

private:
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void beginJob();
  void analyze(edm::Event const&, edm::EventSetup const&);
  void endJob();

  std::string name() const { return "SusyNtuplizer"; }

  void fillLumiSummary(edm::Event const&, edm::EventSetup const&);
  void fillBeamSpot(edm::Event const&, edm::EventSetup const&);
  void fillRhos(edm::Event const&, edm::EventSetup const&);
  void fillTriggerMaps(edm::Event const&, edm::EventSetup const&);
  void fillVertices(edm::Event const&, edm::EventSetup const&);
  void fillGeneralTracks(edm::Event const&, edm::EventSetup const&);
  void fillMetFilters(edm::Event const&, edm::EventSetup const&);
  void fillPFParticles(edm::Event const&, edm::EventSetup const&);
  void fillGenInfo(edm::Event const&, edm::EventSetup const&);
  void fillGenParticles(edm::Event const&, edm::EventSetup const&);
  void fillTriggerEvent(edm::Event const&, edm::EventSetup const&);

  void fillMet(edm::Event const&, edm::EventSetup const&);
  void fillMuons(edm::Event const&, edm::EventSetup const&);
  void fillElectrons(edm::Event const&, edm::EventSetup const&);
  void fillPhotons(edm::Event const&, edm::EventSetup const&);
  void fillCaloJets(edm::Event const&, edm::EventSetup const&);
  void fillPFJets(edm::Event const&, edm::EventSetup const&);
  void fillJPTJets(edm::Event const&, edm::EventSetup const&);

  unsigned fillTrack(reco::TrackRef const&);
  unsigned fillGsfTrack(reco::GsfTrackRef const&);
  unsigned fillSuperCluster(reco::SuperClusterRef const&);
  unsigned fillCluster(reco::CaloClusterPtr const&);
  unsigned fillPFParticle(reco::PFCandidatePtr const&);

  unsigned fillTrackCommon(edm::Ptr<reco::Track> const&, bool&);

  void fillTracksAround(reco::Candidate const&, double, edm::Handle<reco::TrackCollection> const&);

  void finalize();

  // ----------member data ---------------------------

  // InputTags

  std::string lumiSummaryTag_;
  std::string vtxCollectionTag_;
  std::string trackCollectionTag_;
  std::string pfCandidateCollectionTag_;
  std::string genCollectionTag_;
  std::string puSummaryInfoTag_;
  std::string triggerEventTag_;
  std::vector<std::string> muonCollectionTags_;
  std::vector<std::string> electronCollectionTags_;
  std::vector<std::string> photonCollectionTags_;
  std::vector<std::string> caloJetCollectionTags_;
  std::vector<std::string> pfJetCollectionTags_;
  std::vector<std::string> jptJetCollectionTags_;
  std::vector<std::string> metCollectionTags_;
  std::vector<std::string> bTagCollectionTags_;
  std::map<std::string, std::vector<std::string> > muonIdCollectionTags_;
  std::map<std::string, std::vector<std::string> > electronIdCollectionTags_;
  std::map<std::string, std::vector<std::string> > photonIdCollectionTags_;
  std::map<std::string, std::vector<std::string> > photonIsoDepTags_;
  std::map<std::string, std::vector<std::string> > electronIsoDepTags_;
  std::map<std::string, std::vector<std::string> > puJetIdCollectionTags_;
  std::string photonSCRegressionWeights_;

  PropagatorWithMaterial* propagator_;
  TransientTrackBuilder const* transientTrackBuilder_;

  // for HLT prescales
  HLTConfigProvider* hltConfig_;

  // for L1 menu
  L1GtUtils* l1GtUtils_;

  // PFIsolator
  PFIsolationEstimator* isolator03_;

  // text-based event veto for HcalLaser2012 MET filter
  EventFilterFromListStandAlone* hcalLaser2012Filter_;

  // debugLevel
  // 0 : default (no printout from this module)
  // 1 : minimal (function level printing)
  // 2 : print the size of collections
  // 3 : print the values of objects in the collection
  int debugLevel_;

  // flag for storing L1 information
  // default : true
  bool storeL1Info_;

  // flag for storing HLT information
  // default : true
  bool storeHLTInfo_;

  // flag for storing generated informations in ntuples
  // default : false
  bool storeGenInfo_;

  // flag for storing generalTracks collection
  // default : false
  bool storeGeneralTracks_;

  // flag for storing parton flavour matches for ak5pf jets
  // default : true (turned off in runOverAOD.py)
  bool storePFJetPartonMatches_;

  // input RECO mode
  // false : default - reading from AOD
  //         no extrapolation info will be saved in ntuples
  // true : reading from RECO
  bool recoMode_;

  // flag for whether or not MC is fastsim
  // certain collections are not produced by FamosSequence so we skip them
  // default: false (turned on in runOverAOD.py)
  bool isFastSim_;

  // flag for recording TriggerEvent into a separate tree (susyTriggers)
  bool storeTriggerEvents_;

  // electronThreshold
  double electronThreshold_;

  // muonThreshold
  double muonThreshold_;

  // photonThreshold
  double photonThreshold_;

  // jetThreshold will be applied on corrected jets' pt
  double jetThreshold_;

  // PFParticleThreshold
  double pfParticleThreshold_;

  // vector of ak5pf and pdgId matches
  // JetPartonMatcher matches, and we store them in here to match to specific susy::PFJet later
  std::vector< std::pair<reco::Jet, int> > physicsDefinitionMatches;
  std::vector< std::pair<reco::Jet, int> > algorithmicDefinitionMatches;

  //Photon SC energy MVA regression
  EGEnergyCorrector scEnergyCorrector_;

  // since we deal with both reco::Track and reco::GsfTrack, the track reference needs to be polymorphic
  typedef std::map<edm::Ptr<reco::Track>, unsigned> TrackStore;
  typedef std::map<reco::SuperClusterRef, unsigned> SuperClusterStore;
  typedef std::map<reco::CaloClusterPtr, unsigned> CaloClusterStore;
  typedef std::map<reco::PFCandidatePtr, unsigned> PFCandidateStore;
  struct ProductStore {
    void clear() { tracks.clear(); superClusters.clear(); basicClusters.clear(); pfCandidates.clear(); }
    TrackStore tracks;
    SuperClusterStore superClusters;
    CaloClusterStore basicClusters;
    PFCandidateStore pfCandidates;
  } productStore_;

  susy::Event* susyEvent_;
  TTree*       susyTree_;

  susy::TriggerEvent* triggerEvent_;
};


// Constructor - passing parameters, memory allocation to ntuple variables
SusyNtuplizer::SusyNtuplizer(const edm::ParameterSet& iConfig) :
  lumiSummaryTag_(iConfig.getParameter<std::string>("lumiSummaryTag")),
  vtxCollectionTag_(iConfig.getParameter<std::string>("vtxCollectionTag")),
  trackCollectionTag_(iConfig.getParameter<std::string>("trackCollectionTag")),
  pfCandidateCollectionTag_(iConfig.getParameter<std::string>("pfCandidateCollectionTag")),
  genCollectionTag_(iConfig.getParameter<std::string>("genCollectionTag")),
  puSummaryInfoTag_(iConfig.getParameter<std::string>("puSummaryInfoTag")),
  triggerEventTag_(iConfig.getParameter<std::string>("triggerEventTag")),
  muonCollectionTags_(iConfig.getParameter<std::vector<std::string> >("muonCollectionTags")),
  electronCollectionTags_(iConfig.getParameter<std::vector<std::string> >("electronCollectionTags")),
  photonCollectionTags_(iConfig.getParameter<std::vector<std::string> >("photonCollectionTags")),
  caloJetCollectionTags_(iConfig.getParameter<std::vector<std::string> >("caloJetCollectionTags")),
  pfJetCollectionTags_(iConfig.getParameter<std::vector<std::string> >("pfJetCollectionTags")),
  jptJetCollectionTags_(iConfig.getParameter<std::vector<std::string> >("jptJetCollectionTags")),
  metCollectionTags_(iConfig.getParameter<std::vector<std::string> >("metCollectionTags")),
  bTagCollectionTags_(iConfig.getParameter<std::vector<std::string> >("bTagCollectionTags")),
  muonIdCollectionTags_(),
  electronIdCollectionTags_(),
  photonIdCollectionTags_(),
  photonIsoDepTags_(),
  electronIsoDepTags_(),
  puJetIdCollectionTags_(),
  photonSCRegressionWeights_(""),
  propagator_(0),
  transientTrackBuilder_(0),
  hltConfig_(0),
  l1GtUtils_(0),
  isolator03_(0),
  hcalLaser2012Filter_(0),
  debugLevel_(iConfig.getParameter<int>("debugLevel")),
  storeL1Info_(iConfig.getParameter<bool>("storeL1Info")),
  storeHLTInfo_(iConfig.getParameter<bool>("storeHLTInfo")),
  storeGenInfo_(iConfig.getParameter<bool>("storeGenInfo")),
  storeGeneralTracks_(iConfig.getParameter<bool>("storeGeneralTracks")),
  storePFJetPartonMatches_(iConfig.getParameter<bool>("storePFJetPartonMatches")),
  recoMode_(iConfig.getParameter<bool>("recoMode")),
  isFastSim_(iConfig.getParameter<bool>("isFastSim")),
  storeTriggerEvents_(iConfig.getParameter<bool>("storeTriggerEvents")),
  electronThreshold_(iConfig.getParameter<double>("electronThreshold")),
  muonThreshold_(iConfig.getParameter<double>("muonThreshold")),
  photonThreshold_(iConfig.getParameter<double>("photonThreshold")),
  jetThreshold_(iConfig.getParameter<double>("jetThreshold")),
  pfParticleThreshold_(iConfig.getParameter<double>("pfParticleThreshold")),
  physicsDefinitionMatches(),
  algorithmicDefinitionMatches(),
  scEnergyCorrector_(),
  productStore_(),
  susyEvent_(0),
  susyTree_(0),
  triggerEvent_(0)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "ctor";

  if(iConfig.existsAs<std::string>("photonSCRegressionWeights"))
    photonSCRegressionWeights_ = iConfig.getParameter<std::string>("photonSCRegressionWeights");
  else
    photonSCRegressionWeights_ = iConfig.getParameter<edm::FileInPath>("photonSCRegressionWeights").fullPath();

  // check if the file exists
  TFile* dummyFile(TFile::Open(photonSCRegressionWeights_.c_str()));
  if(!dummyFile || dummyFile->IsZombie())
    throw cms::Exception("IOError") << "Photon SC MVA regression weight file " << photonSCRegressionWeights_ << " cannot be opened";
  delete dummyFile;

  edm::ParameterSet const& muonIdTags(iConfig.getParameterSet("muonIdTags"));
  for(std::vector<std::string>::iterator tItr(photonCollectionTags_.begin()); tItr != photonCollectionTags_.end(); ++tItr){
    if(muonIdTags.existsAs<std::vector<std::string> >(*tItr))
      muonIdCollectionTags_[*tItr] = muonIdTags.getParameter<std::vector<std::string> >(*tItr);
  }

  edm::ParameterSet const& electronIdTags(iConfig.getParameterSet("electronIdTags"));
  for(std::vector<std::string>::iterator tItr(photonCollectionTags_.begin()); tItr != photonCollectionTags_.end(); ++tItr){
    if(electronIdTags.existsAs<std::vector<std::string> >(*tItr))
      electronIdCollectionTags_[*tItr] = electronIdTags.getParameter<std::vector<std::string> >(*tItr);
  }

  edm::ParameterSet const& photonIdTags(iConfig.getParameterSet("photonIdTags"));
  for(std::vector<std::string>::iterator tItr(photonCollectionTags_.begin()); tItr != photonCollectionTags_.end(); ++tItr){
    if(photonIdTags.existsAs<std::vector<std::string> >(*tItr))
      photonIdCollectionTags_[*tItr] = photonIdTags.getParameter<std::vector<std::string> >(*tItr);
  }

  edm::ParameterSet const& photonIsoDepTags(iConfig.getParameterSet("photonIsoDepTags"));
  for(std::vector<std::string>::iterator tItr(photonCollectionTags_.begin()); tItr != photonCollectionTags_.end(); ++tItr){
    if(!photonIsoDepTags.existsAs<edm::ParameterSet>(*tItr)) continue;
    edm::ParameterSet const& tagsPSet(photonIsoDepTags.getParameterSet(*tItr));

    std::vector<std::string>& tags(photonIsoDepTags_[*tItr]);
    tags.resize(susy::nPFIsoTypes);

    tags[susy::kChargedHadron] = tagsPSet.getParameter<std::string>("chargedHadron");
    tags[susy::kNeutralHadron] = tagsPSet.getParameter<std::string>("neutralHadron");
    tags[susy::kPhoton] = tagsPSet.getParameter<std::string>("photon");
  }

  edm::ParameterSet const& electronIsoDepTags(iConfig.getParameterSet("electronIsoDepTags"));
  for(std::vector<std::string>::iterator tItr(electronCollectionTags_.begin()); tItr != electronCollectionTags_.end(); ++tItr){
    if(!electronIsoDepTags.existsAs<edm::ParameterSet>(*tItr)) continue;
    edm::ParameterSet const& tagsPSet(electronIsoDepTags.getParameterSet(*tItr));

    std::vector<std::string>& tags(electronIsoDepTags_[*tItr]);
    tags.resize(susy::nPFIsoTypes);

    tags[susy::kChargedHadron] = tagsPSet.getParameter<std::string>("chargedHadron");
    tags[susy::kNeutralHadron] = tagsPSet.getParameter<std::string>("neutralHadron");
    tags[susy::kPhoton] = tagsPSet.getParameter<std::string>("photon");
  }

  edm::ParameterSet const& puJetIdTags(iConfig.getParameterSet("puJetIdTags"));
  for(std::vector<std::string>::iterator tItr(pfJetCollectionTags_.begin()); tItr != pfJetCollectionTags_.end(); ++tItr){
    if(puJetIdTags.existsAs<std::vector<std::string> >(*tItr))
      puJetIdCollectionTags_[*tItr] = puJetIdTags.getParameter<std::vector<std::string> >(*tItr);
  }

  TString outputFileName(iConfig.getParameter<std::string>("outputFileName"));

  TFile* outF(TFile::Open(outputFileName, "RECREATE"));
  if(!outF || outF->IsZombie())
    throw cms::Exception("IOError") << "Cannot create file " << outputFileName;

  outF->cd();    
  susyEvent_ = new susy::Event;
  susyTree_ = new TTree("susyTree", "SUSY Event");
  susyTree_->Branch("susyEvent", "susy::Event", &susyEvent_);
  susyTree_->SetAutoSave(10000000); // 10M events

  if(storeTriggerEvents_){
    TString triggerFileName(iConfig.getParameter<std::string>("triggerFileName"));

    triggerEvent_ = new susy::TriggerEvent;
    triggerEvent_->bookTrees(triggerFileName);
  }
}

SusyNtuplizer::~SusyNtuplizer()
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "dtor";
}

// ------------ method called once each job just before starting event loop  ------------
void
SusyNtuplizer::beginJob()
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "beginJob";

  if(photonCollectionTags_.size() != 0){
    // PFIsolator init
    isolator03_ = new PFIsolationEstimator;
    isolator03_->initializePhotonIsolation(true); // bool applyVeto
    isolator03_->setConeSize(0.3);
  }

  l1GtUtils_ = new L1GtUtils;

  hltConfig_ = new HLTConfigProvider;

  std::string hcalLaserFilterFile(std::getenv("CMSSW_BASE"));
  hcalLaserFilterFile += "/src/EventFilter/HcalRawToDigi/data/HCALLaser2012AllDatasets.txt.gz";
  gzFile dummyFP(gzopen(hcalLaserFilterFile.c_str(), "r"));
  if(dummyFP != 0){
    gzclose(dummyFP);
    hcalLaser2012Filter_ = new EventFilterFromListStandAlone(hcalLaserFilterFile);
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void
SusyNtuplizer::endJob()
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "endJob";

  finalize();
}

// ---- method called once each job just before starting event loop  ---
void
SusyNtuplizer::beginRun(edm::Run const& iRun, edm::EventSetup const& _eventSetup)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "beginRun";

  try{
    if(storeL1Info_)
      l1GtUtils_->getL1GtRunCache(iRun, _eventSetup, true, false); // use event setup, do not use L1TriggerMenuLite

    if(storeHLTInfo_){
      //intialize HLTConfigProvider
      bool menuChanged;
      if(!hltConfig_->init(iRun, _eventSetup, "HLT", menuChanged))
        throw cms::Exception("RuntimeError") << "HLTConfigProvider::init() returned non 0";

      if(menuChanged) edm::LogInfo(name()) << "beginRun: HLT configuration changed to " << hltConfig_->tableName();
    }

    if(photonCollectionTags_.size() != 0){
      // initialize Photon SC energy MVA regression
      // EventSetup is not used if the third argument is false
      scEnergyCorrector_.Initialize(_eventSetup, photonSCRegressionWeights_, false);
    }
  }
  catch(...){
    finalize();
    throw;
  }
}

// ------------ method called to for each event  ------------
// fill the tree variables
void
SusyNtuplizer::analyze(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "analyze";

  // Event processing is in one big try block to exit gracefully in case of exception
  try{

    if(recoMode_) {
      if(debugLevel_ > 0) edm::LogInfo(name()) << "analyze: setup track extrapolation tools";

      edm::ESHandle<MagneticField> fieldHndl;
      edm::ESHandle<TransientTrackBuilder> builderHndl;
      _eventSetup.get<IdealMagneticFieldRecord>().get(fieldHndl);
      _eventSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builderHndl);

      propagator_ = new PropagatorWithMaterial(alongMomentum, 0.000511, fieldHndl.product());
      transientTrackBuilder_ = builderHndl.product();
    }

    if(!_event.isRealData() && storePFJetPartonMatches_){
      if(debugLevel_ > 0) edm::LogInfo(name()) << "analyze: fill ak5pf jet parton matches";

      physicsDefinitionMatches.clear();
      algorithmicDefinitionMatches.clear();

      edm::Handle<reco::JetMatchedPartonsCollection> matchCollH;
      _event.getByLabel("flavourByRef", matchCollH);

      for(reco::JetMatchedPartonsCollection::const_iterator p = matchCollH->begin(); p != matchCollH->end(); p++) {
        const reco::Jet *aJet = (*p).first.get();
        const reco::MatchedPartons aMatch = (*p).second;

        const reco::GenParticleRef thePhyDef = aMatch.physicsDefinitionParton();
        if(thePhyDef.isNonnull()) {
          physicsDefinitionMatches.push_back( make_pair(*aJet, (int)(thePhyDef.get()->pdgId())) );
        }

        const reco::GenParticleRef theAlgDef = aMatch.algoDefinitionParton();
        if(theAlgDef.isNonnull()) {
          algorithmicDefinitionMatches.push_back( make_pair(*aJet, (int)(theAlgDef.get()->pdgId())) );
        }
      }
    } // if !isRealData

    // reset map from RECO objects to susy objects
    productStore_.clear();

    // initialize susyEvent object
    susyEvent_->Init();

    if(debugLevel_ > 0) edm::LogInfo(name()) << "analyze: fill event info" << std::endl;

    susyEvent_->isRealData = _event.isRealData() ? 1 : 0;
    susyEvent_->runNumber = _event.id().run();
    susyEvent_->eventNumber = _event.id().event();
    susyEvent_->luminosityBlockNumber = _event.luminosityBlock();
    susyEvent_->bunchCrossing = _event.bunchCrossing();

    if(debugLevel_ > 1) edm::LogInfo(name()) << "analyze: run " << _event.id().run()
                                             << ", event " << _event.id().event()
                                             << ", isRealData " << _event.isRealData()
                                             << ", lumiBlock " << _event.luminosityBlock();

    fillLumiSummary(_event, _eventSetup);

    fillBeamSpot(_event, _eventSetup);

    fillRhos(_event, _eventSetup);

    fillTriggerMaps(_event, _eventSetup);

    fillGeneralTracks(_event, _eventSetup);

    fillMetFilters(_event, _eventSetup);

    fillPFParticles(_event, _eventSetup);

    fillMet(_event, _eventSetup);

    fillPhotons(_event, _eventSetup);

    fillElectrons(_event, _eventSetup);

    fillMuons(_event, _eventSetup);

    fillCaloJets(_event, _eventSetup);

    fillPFJets(_event, _eventSetup);

    //  fillJPTJets(_event, _eventSetup);

    fillVertices(_event, _eventSetup);

    fillGenInfo(_event, _eventSetup);

    fillGenParticles(_event, _eventSetup);

    fillTriggerEvent(_event, _eventSetup);

  }
  catch(...){
    delete propagator_;
    finalize();
    throw;
  }

  delete propagator_;
  propagator_ = 0;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "analyze: fill tree at last";

  susyTree_->Fill();
}

void
SusyNtuplizer::fillLumiSummary(edm::Event const& _event, edm::EventSetup const&)
{
  // lumiSummary only available in data
  if(susyEvent_->isRealData == 0) return;

  if(lumiSummaryTag_ == "") return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillLumiSummary";

  edm::LuminosityBlock const& lumiBlock(_event.getLuminosityBlock());
  edm::Handle<LumiSummary> lsH;

  lumiBlock.getByLabel(edm::InputTag(lumiSummaryTag_), lsH);
  susyEvent_->avgInsRecLumi = lsH->avgInsRecLumi();
  susyEvent_->intgRecLumi = lsH->intgRecLumi();
}

void
SusyNtuplizer::fillBeamSpot(edm::Event const& _event, edm::EventSetup const&)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillBeamSpot";

  edm::Handle<reco::BeamSpot> bsh;
  _event.getByLabel("offlineBeamSpot", bsh);
  susyEvent_->beamSpot.SetXYZ(bsh->position().x(), bsh->position().y(), bsh->position().z());

  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillBeamSpot: beamSpot : " << bsh->position().x() << ", "
                                           << bsh->position().y() << ", " << bsh->position().z();
}

void
SusyNtuplizer::fillRhos(edm::Event const& _event, edm::EventSetup const&)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillRhos";

  edm::Handle<double> rhoH;

  _event.getByLabel("kt6PFJets", "rho", rhoH);
  susyEvent_->rho = *rhoH;

  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillRhos: rho calculated from kt6PFJets = " << susyEvent_->rho;

  _event.getByLabel("kt6PFJetsRhoBarrelOnly", "rho", rhoH);
  susyEvent_->rhoBarrel = *rhoH;

  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillRhos: rho calculated from kt6PFJetsRhoBarrelOnly = " << susyEvent_->rhoBarrel;

  _event.getByLabel("kt6PFJetsRho25", "rho", rhoH);
  susyEvent_->rho25 = *rhoH;

  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillRhos: rho calculated from kt6PFJetsRho25 = " << susyEvent_->rho25;
}

void
SusyNtuplizer::fillTriggerMaps(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillTriggerMaps";

  // L1 Info only available in data and FullSim
  if(storeL1Info_ && (susyEvent_->isRealData || !isFastSim_)){
    if(debugLevel_ > 0) edm::LogInfo(name()) << "fillTriggerMaps: L1 map";

    // Get and cache L1 menu
    l1GtUtils_->getL1GtRunCache(_event, _eventSetup, true, false); // use event setup, do not use L1TriggerMenuLite

    int err;
    L1GtTriggerMenu const* menu(l1GtUtils_->ptrL1TriggerMenuEventSetup(err));
    if(err != 0)
      throw cms::Exception("RuntimeError") << "L1GtUtils failed to return the trigger menu";

    AlgorithmMap const& l1GtAlgorithms(menu->gtAlgorithmMap());

    for(CItAlgo iAlgo(l1GtAlgorithms.begin()); iAlgo != l1GtAlgorithms.end(); ++iAlgo){
      std::string const& algoName(iAlgo->second.algoName());

      if(iAlgo->second.algoBitNumber() >= int(L1GlobalTriggerReadoutSetup::NumberPhysTriggers)){
        edm::LogError("DataIntegrity") << "L1 physics algorithm '" << algoName << "' has bit number " << iAlgo->second.algoBitNumber() << " >= " << L1GlobalTriggerReadoutSetup::NumberPhysTriggers << "\n"
                                       << "Skipping";
        continue;
      }

      L1GtUtils::TriggerCategory category;
      int bit;
      if(!l1GtUtils_->l1AlgoTechTrigBitNumber(algoName, category, bit)){
        edm::LogError("DataIntegrity") << "L1 physics algorithm '" << algoName << "' not found in the L1 menu\n"
                                       << "Skipping";
        continue;
      }
      if(category != L1GtUtils::AlgorithmTrigger){
        edm::LogError("DataIntegrity") << "L1 physics algorithm '" << algoName << "' does not have category 'AlgorithmTrigger' from 'L1GtUtils'\n"
                                       << "Skipping";
        continue;
      }

      bool decisionBeforeMask;
      bool decisionAfterMask;
      int prescale;
      int mask;
      int error(l1GtUtils_->l1Results(_event, algoName, decisionBeforeMask, decisionAfterMask, prescale, mask));
      if(error != 0){
        edm::LogError("DataIntegrity") << "L1 physics algorithm '" << algoName << "' decision has error code " << error << " from 'L1GtUtils'\n"
                                       << "Skipping";
        continue;
      }

      susyEvent_->l1Map[algoName] = std::pair<Int_t, UChar_t>(prescale, UChar_t(decisionBeforeMask));
    }

  } // L1 only if data or fullsim

  if(!storeHLTInfo_) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillTriggerMaps: HLT map";

  edm::Handle<edm::TriggerResults> hltH;
  _event.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), hltH);

  unsigned nHlt(hltH->size());

  edm::TriggerNames const& hltTriggerNames(_event.triggerNames(*hltH));
  if(nHlt != hltTriggerNames.size())
    throw cms::Exception("RuntimeError") << "TriggerPathName size mismatches !!!";

  // loop over hlt paths
  for(unsigned i(0); i < nHlt; ++i){
    int prescale(hltConfig_->prescaleValue(_event, _eventSetup, hltTriggerNames.triggerName(i)));

    susyEvent_->hltMap[hltTriggerNames.triggerName(i)] = std::pair<Int_t, UChar_t>(prescale, UChar_t(hltH->accept(i)));

    if(debugLevel_ > 2) edm::LogInfo(name()) << hltTriggerNames.triggerName(i) << " : " << hltH->accept(i);
  }
}

void
SusyNtuplizer::fillVertices(edm::Event const& _event, edm::EventSetup const&)
{
  if(vtxCollectionTag_ == "") return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillVertices";

  edm::Handle<reco::VertexCollection> vtxH;
  _event.getByLabel(edm::InputTag(vtxCollectionTag_), vtxH);
  reco::VertexCollection const& vertices(*vtxH);

  susyEvent_->vertices.resize(vertices.size());

  for(unsigned iV(0); iV != vertices.size(); ++iV){
    reco::Vertex const& recoVtx(vertices[iV]);
    susy::Vertex& susyVtx(susyEvent_->vertices[iV]);

    susyVtx.chi2       = recoVtx.chi2();
    susyVtx.ndof       = recoVtx.ndof();
    susyVtx.tracksSize = recoVtx.tracksSize();
    susyVtx.position.SetXYZ(recoVtx.x(), recoVtx.y(), recoVtx.z());

    if(debugLevel_ > 2) edm::LogInfo(name()) << "vtx" << iV << " : " << recoVtx.x() << ", " << recoVtx.y() << ", " << recoVtx.z();

    for(reco::Vertex::trackRef_iterator trkItr(recoVtx.tracks_begin()); trkItr != recoVtx.tracks_end(); ++trkItr){
      TrackStore::const_iterator itr(productStore_.tracks.find(edm::refToPtr(trkItr->castTo<edm::Ref<reco::TrackCollection> >())));
      if(itr != productStore_.tracks.end())
        susyEvent_->tracks[itr->second].vertexIndex = iV;
    }
  }
}

void
SusyNtuplizer::fillGeneralTracks(edm::Event const& _event, edm::EventSetup const&)
{
  if(trackCollectionTag_ == "") return;

  if(!storeGeneralTracks_) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillGeneralTracks";

  edm::Handle<reco::TrackCollection> trackH;
  _event.getByLabel(edm::InputTag(trackCollectionTag_), trackH);

  unsigned iTrk(0);
  for(reco::TrackCollection::const_iterator tItr(trackH->begin()); tItr != trackH->end(); ++tItr, ++iTrk){
    if(tItr->pt() < 1.0) continue;
    fillTrack(reco::TrackRef(trackH, iTrk));
  }
}

void
SusyNtuplizer::fillMetFilters(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillMetFilters";

  // names must be in the exact same order as enum susy::MetFilters
  std::string filterNames[] = {
    "BeamHaloSummary",
    "HBHENoiseFilterResultProducer:HBHENoiseFilterResult",
    "EcalDeadCellTriggerPrimitiveFilter",
    "EcalDeadCellBoundaryEnergyFilter",
    "hcalLaserEventFilter",
    "trackingFailureFilter",
    "eeBadScFilter",
    "", // HcalLaser1012 is taken from an external list
    "ecalLaserCorrFilter",
    "manystripclus53X",
    "toomanystripclus53X",
    "logErrorTooManyClusters",
    "eeNoiseFilter",
    "inconsistentMuonPFCandidateFilter",
    "greedyMuonPFCandidateFilter"
  };

  std::set<unsigned> notForFastSim;
  notForFastSim.insert(susy::kHcalNoise);

  std::set<unsigned> irregular;
  irregular.insert(susy::kCSCBeamHalo);
  irregular.insert(susy::kHcalLaser2012);

  std::vector<bool> pass(susy::nMetFilters, false);

  edm::Handle<bool> boolH;
  for(unsigned iF(0); iF != susy::nMetFilters; ++iF){
    // skip filters that are not a simple bool
    if(irregular.find(iF) != irregular.end()) continue;

    // skip filters that are not FastSim-compatible
    if(!susyEvent_->isRealData && isFastSim_ && notForFastSim.find(iF) != notForFastSim.end()) continue;

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: " << filterNames[iF];

    _event.getByLabel(edm::InputTag(filterNames[iF]), boolH);
    pass[iF] = *boolH;
  }

  // BeamHaloSummary only available in data and FullSim
  if(susyEvent_->isRealData || !isFastSim_){
    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: " << filterNames[susy::kCSCBeamHalo];

    edm::Handle<reco::BeamHaloSummary> beamHaloSummary;
    _event.getByLabel(edm::InputTag(filterNames[susy::kCSCBeamHalo]), beamHaloSummary);
    pass[susy::kCSCBeamHalo] = !(beamHaloSummary->CSCTightHaloId());
  }

  if(hcalLaser2012Filter_){
    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: " << filterNames[susy::kHcalLaser2012];
    pass[susy::kHcalLaser2012] = hcalLaser2012Filter_->filter(_event.id().run(), _event.luminosityBlock(), _event.id().event());
  }

  for(unsigned iF(0); iF != susy::nMetFilters; ++iF)
    susyEvent_->metFilterBit |= (pass[iF] ? 1 << iF : 0);
}

void
SusyNtuplizer::fillPFParticles(edm::Event const& _event, edm::EventSetup const&)
{
  if(pfCandidateCollectionTag_ == "") return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillPFParticles";

  edm::Handle<reco::PFCandidateCollection> pfH;
  _event.getByLabel(edm::InputTag(pfCandidateCollectionTag_), pfH);

  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillPFParticles: size of PFCandidateCollection = " << pfH->size();

  unsigned iPart(0);
  for(reco::PFCandidateCollection::const_iterator pItr(pfH->begin()); pItr != pfH->end(); ++pItr, ++iPart){
    if(pItr->pt() < pfParticleThreshold_) continue;
    fillPFParticle(reco::PFCandidatePtr(pfH, iPart));

    if(debugLevel_ > 2) edm::LogInfo(name()) << "e, px, py, pz = " << pItr->energy() << ", "
                                             << pItr->px() << ", " << pItr->py() << ", " << pItr->pz();
  } // for
}

void
SusyNtuplizer::fillGenInfo(edm::Event const& _event, edm::EventSetup const&)
{
  if(!storeGenInfo_ || _event.isRealData()) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillGenInfo";

  // event weighting variables, for example ptHat
  edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
  if(_event.getByLabel("generator", GenEventInfoHandle) && GenEventInfoHandle->binningValues().size() > 0)
    susyEvent_->gridParams["ptHat"] = GenEventInfoHandle->binningValues()[0];

  //get PU summary info
  edm::Handle<std::vector<PileupSummaryInfo> > pPUSummaryInfo;
  _event.getByLabel(edm::InputTag(puSummaryInfoTag_), pPUSummaryInfo);

  //fill PUSummaryInfo object
  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillGenInfo: fill PU summary information";
  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillGenInfo: size of PileupSummaryInfo collection: " << pPUSummaryInfo->size();

  susyEvent_->pu.resize(pPUSummaryInfo->size());

  for(std::vector<PileupSummaryInfo>::const_iterator iPU(pPUSummaryInfo->begin()); iPU != pPUSummaryInfo->end(); ++iPU){
    unsigned int index = iPU - pPUSummaryInfo->begin();

    if(debugLevel_ > 1)
      edm::LogInfo(name()) << "size of z position collection for BX " << index << ": " << iPU->getPU_zpositions().size() << "\n"
                           << "size of sum pT / low pT collection for BX " << index << ": " << iPU->getPU_sumpT_lowpT().size() << "\n"
                           << "size of sum pT / high pT collection for BX " << index << ": " << iPU->getPU_sumpT_highpT().size() << "\n"
                           << "size of track / low pT collection for BX " << index << ": " << iPU->getPU_ntrks_lowpT().size() << "\n"
                           << "size of track / high pT collection for BX " << index << ": " << iPU->getPU_ntrks_highpT().size() << "\n"
                           << "size of inst. lumi collection for BX " << index << ": " << iPU->getPU_instLumi().size() << "\n"
                           << "size of DataMixer event collection for BX " << index << ": " << iPU->getPU_EventID().size();

    if(debugLevel_ > 2)
      edm::LogInfo(name()) << "No. interactions for BX " << index << ": " << iPU->getPU_NumInteractions() << "\n"
                           << "BX ID for BX " << index << ": " << iPU->getBunchCrossing() << "\n"
                           << "True no. interactions for BX " << index << ": " << iPU->getTrueNumInteractions();
    //for MC samples earlier than CMSSWv4.2.8 or Fall11, comment out the last line

    susy::PUSummaryInfo& PUInfoForThisBX(susyEvent_->pu[index]);

    PUInfoForThisBX.numInteractions = iPU->getPU_NumInteractions();
    for(std::vector<float>::const_iterator i = iPU->getPU_zpositions().begin(); i != iPU->getPU_zpositions().end(); ++i){
      PUInfoForThisBX.zPositions.push_back(*i);
      if(debugLevel_ > 2) edm::LogInfo(name()) << "Z position for BX " << index
                                               << ", interaction " << (i - iPU->getPU_zpositions().begin()) << ": " << *i;
    }
    for(std::vector<float>::const_iterator i = iPU->getPU_sumpT_lowpT().begin(); i != iPU->getPU_sumpT_lowpT().end(); ++i){
      PUInfoForThisBX.sumPTLowPT.push_back(*i);
      if(debugLevel_ > 2) edm::LogInfo(name()) << "Sum pT (low pT) for BX " << index
                                               << ", interaction " << (i - iPU->getPU_sumpT_lowpT().begin()) << ": " << *i;
    }
    for(std::vector<float>::const_iterator i = iPU->getPU_sumpT_highpT().begin(); i != iPU->getPU_sumpT_highpT().end(); ++i){
      PUInfoForThisBX.sumPTHighPT.push_back(*i);
      if(debugLevel_ > 2) edm::LogInfo(name()) << "Sum pT (high pT) for BX " << index
                                               << ", interaction " << (i - iPU->getPU_sumpT_highpT().begin()) << ": " << *i;
    }
    for(std::vector<int>::const_iterator i = iPU->getPU_ntrks_lowpT().begin(); i != iPU->getPU_ntrks_lowpT().end(); ++i){
      PUInfoForThisBX.numTracksLowPT.push_back(*i);
      if(debugLevel_ > 2) edm::LogInfo(name()) << "No. tracks (low pT) for BX " << index
                                               << ", interaction " << (i - iPU->getPU_ntrks_lowpT().begin()) << ": " << *i;
    }
    for(std::vector<int>::const_iterator i = iPU->getPU_ntrks_highpT().begin(); i != iPU->getPU_ntrks_highpT().end(); ++i){
      PUInfoForThisBX.numTracksHighPT.push_back(*i);
      if(debugLevel_ > 2) edm::LogInfo(name()) << "No. tracks (high pT) for BX " << index
                                               << ", interaction " << (i - iPU->getPU_ntrks_highpT().begin()) << ": " << *i;
    }
    for(std::vector<float>::const_iterator i = iPU->getPU_instLumi().begin(); i != iPU->getPU_instLumi().end(); ++i){
      PUInfoForThisBX.instLumi.push_back(*i);
      if(debugLevel_ > 2) edm::LogInfo(name()) << "Inst. lumi for BX " << index
                                               << ", interaction " << (i - iPU->getPU_instLumi().begin()) << ": " << *i;
    }
    for(std::vector<edm::EventID>::const_iterator i = iPU->getPU_EventID().begin(); i != iPU->getPU_EventID().end(); ++i){
      PUInfoForThisBX.dataMixerRun.push_back(i->run());
      PUInfoForThisBX.dataMixerEvt.push_back(i->event());
      PUInfoForThisBX.dataMixerLumiSection.push_back(i->luminosityBlock());
      if(debugLevel_ > 2) edm::LogInfo(name()) << "Event ID for BX " << index
                                               << ", interaction " << (i - iPU->getPU_EventID().begin()) << ":\n"
                                               << " Run: " << i->run() << "\n"
                                               << " Event: " << i->event() << "\n"
                                               << " Lumi section: " << i->luminosityBlock();
    }
    PUInfoForThisBX.BX = iPU->getBunchCrossing();
    //for MC samples earlier than CMSSWv4.2.8 or Fall11, COMMENT OUT the following line
    PUInfoForThisBX.trueNumInteractions = iPU->getTrueNumInteractions();
    //for MC samples earlier than CMSSWv4.2.8 or Fall11, UNCOMMENT the following line
    //PUInfoForThisBX.trueNumInteractions = -1.0;
  }
}

void
SusyNtuplizer::fillGenParticles(edm::Event const& _event, edm::EventSetup const&)
{
  if(!storeGenInfo_ || _event.isRealData()) return;

  // Store the skimmed full decay tree
  // Only the ancestors of final state (status == 1) particles with pt > 2 GeV will be stored
  // Hadronic "blobs" are cleaned out
  // In case of loops (a particle has multiple mothers), arbitration based on the particle species and Pt is performed

  // The following algorithm will be much easier if we have a couple of auxiliary recursive functions..

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillGenParticles";

  edm::Handle<reco::GenParticleCollection> gpH;
  _event.getByLabel(edm::InputTag(genCollectionTag_), gpH);

  std::set<reco::GenParticle const*> gps;

  // Save all particles that might appear in the final decay tree
  for(reco::GenParticleCollection::const_iterator it(gpH->begin()); it != gpH->end(); it++){
    if(it->status() == 1 && it->pt() < 2.) continue;
    if(it->numberOfDaughters() == 0 && it->status() != 1) continue;
    gps.insert(&*it);
  }

  std::map<reco::GenParticle const*, reco::GenParticle const*> motherMap;
  std::map<reco::GenParticle const*, std::list<reco::GenParticle const*> > daughterMap;

  // Loop over gps, establish mother-daughter relationships + perform arbitration when necessary
  std::set<reco::GenParticle const*>::const_iterator gpend(gps.end());
  std::list<reco::GenParticle const*> rootNodes;
  for(std::set<reco::GenParticle const*>::const_iterator gItr(gps.begin()); gItr != gpend; ++gItr){
    reco::GenParticle const* part(*gItr);
    if(part->numberOfMothers() == 0){
      rootNodes.push_back(part);
      motherMap[part] = 0;
    }

    unsigned nD(part->numberOfDaughters());
    for(unsigned iD(0); iD < nD; iD++){
      reco::GenParticle const* daughter(static_cast<reco::GenParticle const*>(part->daughter(iD)));
      if(gps.find(daughter) == gpend) continue;

      if(motherMap.find(daughter) != motherMap.end()){
        // the daughter is claimed by some other mother - arbtrate
        reco::GenParticle const* existing(motherMap[daughter]);

        int thispdg(std::abs(part->pdgId()));
        bool thishad((thispdg / 100) % 10 != 0 || thispdg == 21 || (thispdg > 80 && thispdg < 101));
        int pdg(std::abs(existing->pdgId()));
        bool had((pdg / 100) % 10 != 0 || pdg == 21 || (pdg > 80 && pdg < 101));

        bool takeAway(false);
        if((thishad && had) || (!thishad && !had))
          takeAway = part->pt() > existing->pt();
        else if(!thishad && had)
          takeAway = true;

        if(!takeAway) continue;

        std::list<reco::GenParticle const*>& daughters(daughterMap[existing]);
        daughters.remove(daughter);
      }

      motherMap[daughter] = part;
      daughterMap[part].push_back(daughter);
    }
  }

  if(rootNodes.size() == 0) return;

  // flag: whether subtree with the particle as the root node is already cleaned
  std::map<reco::GenParticle const*, bool> cleanMap;
  for(std::set<reco::GenParticle const*>::const_iterator gpItr(gps.begin()); gpItr != gpend; ++gpItr)
    cleanMap[*gpItr] = false;

  std::list<reco::GenParticle const*>* sisters(&rootNodes);
  std::list<reco::GenParticle const*>::iterator pItr(sisters->begin());
  std::list<reco::GenParticle const*>::iterator pEnd(sisters->end());

  // start from the root node list, recursively clean the full tree
  // strategy: 1. clean the daughters (recursive); 2. once no daughters are left or all are cleaned, clean myself; 3. go up
  while(true){
    reco::GenParticle const* part(*pItr);
    reco::GenParticle const* mother(motherMap[part]);
    std::list<reco::GenParticle const*>* daughters(&daughterMap[part]);
    while(!cleanMap[part] && daughters->size() > 0){
      mother = part;
      sisters = daughters;
      pItr = sisters->begin();
      pEnd = sisters->end();
      part = *pItr;
      daughters = &daughterMap[part];
    }

    unsigned nD(daughters->size());
    unsigned pdg(std::abs(part->pdgId()));
    unsigned motherPdg(mother ? std::abs(mother->pdgId()) : 0);

    bool intermediateTerminal(nD == 0 && part->status() != 1);
    bool noDecay(nD == 1 && part->pdgId() == daughters->front()->pdgId());
    bool hadronicIntermediate(mother && part->status() != 1 && ((pdg / 100) % 10 != 0 || pdg == 21 || (pdg > 80 && pdg < 101)));
    bool firstHeavyHadron((motherPdg / 1000) % 10 < 4 && (motherPdg / 100) % 10 < 4 && ((pdg / 1000) % 10 >= 4 || (pdg / 100) % 10 >= 4));
    bool lightFromLight(motherPdg < 4 && pdg < 4);

    if(intermediateTerminal || noDecay || (hadronicIntermediate && !firstHeavyHadron) || lightFromLight){
      // The particle should be removed, its daughters appended to its siblings list, and their mothers set to particle's mother
      std::list<reco::GenParticle const*>::iterator dEnd(daughters->end());
      for(std::list<reco::GenParticle const*>::iterator dItr(daughters->begin()); dItr != dEnd; ++dItr){
        sisters->push_back(*dItr);
        motherMap[*dItr] = mother;
      }
      gps.erase(part);
      --pItr;
      sisters->remove(part);
      pEnd = sisters->end();
    }

    ++pItr;

    if(pItr == pEnd){
      // reached the end of the list of siblings; go up one level in hierarchy

      if(!mother) break; // we were already in the list of root nodes; cleaning is done

      part = mother;
      mother = motherMap[part];
      if(mother)
        sisters = &daughterMap[mother];
      else // I am a root node
        sisters = &rootNodes;

      pItr = std::find(sisters->begin(), sisters->end(), part);
      pEnd = sisters->end();
      cleanMap[part] = true;
    }
  }

  // serialize the particle list
  std::vector<reco::GenParticle const*> gpv;
  gpend = gps.end();
  for(std::set<reco::GenParticle const*>::const_iterator gItr(gps.begin()); gItr != gpend; ++gItr)
    gpv.push_back(*gItr);

  unsigned nP(gpv.size());
  susyEvent_->genParticles.resize(nP);

  for(unsigned iP(0); iP < nP; iP++){
    reco::GenParticle const* gp(gpv[iP]);
    susy::Particle& susyPart(susyEvent_->genParticles[iP]);

    short iM(-1);
    std::vector<reco::GenParticle const*>::const_iterator mItr(std::find(gpv.begin(), gpv.end(), motherMap[gp]));
    if(mItr != gpv.end()) iM = mItr - gpv.begin();

    susyPart.motherIndex = iM;
    susyPart.status = gp->status();
    susyPart.pdgId = gp->pdgId();
    susyPart.charge = gp->charge();
    susyPart.vertex.SetXYZ(gp->vx(), gp->vy(), gp->vz());
    susyPart.momentum.SetXYZT(gp->px(), gp->py(), gp->pz(), gp->energy());
  }
}

void
SusyNtuplizer::fillTriggerEvent(edm::Event const& _event, edm::EventSetup const&)
{
  if(!storeTriggerEvents_) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillTriggerEvent";

  triggerEvent_->initializeEvent(_event.id().run(), _event.id().event());

  try{
    edm::Handle<trigger::TriggerEvent> teH;
    _event.getByLabel(triggerEventTag_, teH);

    trigger::TriggerObjectCollection const& objects(teH->getObjects());
    unsigned nObj(objects.size());
    for(unsigned iObj(0); iObj < nObj; ++iObj){
      trigger::TriggerObject const& obj(objects.at(iObj));
      triggerEvent_->fillObject(susy::TriggerObject(obj.pt(), obj.eta(), obj.phi(), obj.mass()));
    }

    unsigned nF(teH->sizeFilters());
    for(unsigned iF(0); iF < nF; ++iF){
      trigger::Vids const& vids(teH->filterIds(iF));
      trigger::Keys const& keys(teH->filterKeys(iF));
      triggerEvent_->fillFilter(teH->filterTag(iF).label(), vids, keys);
    }
  }
  catch(...){
    triggerEvent_->fillEvent();
    throw;
  }

  triggerEvent_->fillEvent();
}

void
SusyNtuplizer::fillMet(edm::Event const& _event, edm::EventSetup const&)
{
  if(metCollectionTags_.size() == 0) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillMet";

  for(unsigned iCol(0); iCol != metCollectionTags_.size(); ++iCol){
    std::string& colName(metCollectionTags_[iCol]);

    // skip genMet if real data
    if(colName.find("gen") != std::string::npos && _event.isRealData()) continue;

    edm::Handle<edm::View<reco::MET> > metH;
    _event.getByLabel(edm::InputTag(colName), metH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMet: size of MET coll " << colName << " = " << metH->size();

    reco::MET const& recoMet(*(metH->begin()));
    susy::MET& susyMet(susyEvent_->metMap[colName]);

    susyMet.mEt.Set(recoMet.p4().px(),recoMet.p4().py());
    susyMet.sumEt = recoMet.sumEt();
    susyMet.significance = recoMet.significance();

    if(debugLevel_ > 2) edm::LogInfo(name()) << "met, metX, metY, sumEt, significance : "
                                             << susyMet.mEt.Mod() << ", " << susyMet.mEt.X() << ", " << susyMet.mEt.Y() << ", "
                                             << susyMet.sumEt << ", " << susyMet.significance;

    unsigned nCorr(recoMet.mEtCorr().size());
    susyMet.mEtCorr.resize(nCorr);

    for(unsigned i(0); i < nCorr; ++i){
      susy::CorrMETData& corr(susyMet.mEtCorr[i]);

      corr.dmEx          = recoMet.dmEx()[i];
      corr.dmEy          = recoMet.dmEy()[i];
      corr.dsumEt        = recoMet.dsumEt()[i];
      corr.dSignificance = recoMet.dSignificance()[i];

      if(debugLevel_ > 2) edm::LogInfo(name()) << "dmEx, dmEy, dsumEt, dSignificance : "
                                               << corr.dmEx << ", " << corr.dmEy << ", "
                                               << corr.dsumEt << ", " << corr.dSignificance;
    }//for
  }
}

void
SusyNtuplizer::fillPhotons(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(photonCollectionTags_.size() == 0) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillPhotons";

  edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
  edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;

  _event.getByLabel("reducedEcalRecHitsEB", barrelRecHitsHandle);
  _event.getByLabel("reducedEcalRecHitsEE", endcapRecHitsHandle);

  edm::ESHandle<CaloGeometry> cgH;
  _eventSetup.get<CaloGeometryRecord>().get(cgH);
  CaloGeometry const* caloGeometry(cgH.product());

  edm::ESHandle<CaloTopology> ctH;
  _eventSetup.get<CaloTopologyRecord>().get(ctH);
  CaloTopology const* caloTopology(ctH.product());

  edm::Handle<reco::ConversionCollection> hVetoConversions;
  _event.getByLabel("allConversions", hVetoConversions);

  edm::Handle<reco::GsfElectronCollection> hVetoElectrons;
  _event.getByLabel("gsfElectrons", hVetoElectrons);

  edm::Handle<reco::PFCandidateCollection> pfH;
  _event.getByLabel("particleFlow", pfH);

  edm::Handle<reco::VertexCollection> vtxH;
  _event.getByLabel(edm::InputTag(vtxCollectionTag_), vtxH);
  reco::VertexRef primVtxRef(vtxH, 0);

  edm::Handle<reco::VertexCollection> vtxWithBSH;
  _event.getByLabel("offlinePrimaryVerticesWithBS", vtxWithBSH);

  edm::Handle<reco::TrackCollection> trkH;
  _event.getByLabel(edm::InputTag(trackCollectionTag_), trkH);

  edm::Handle<edm::ValueMap<PFCandidatePtr> > pfTranslationH;
  _event.getByLabel("particleFlow", "photons", pfTranslationH);

  math::XYZPoint beamSpot(susyEvent_->beamSpot.X(), susyEvent_->beamSpot.Y(), susyEvent_->beamSpot.Z());

  EcalClusterLazyTools lazyTools(_event, _eventSetup, edm::InputTag("reducedEcalRecHitsEB"), edm::InputTag("reducedEcalRecHitsEE"));

  edm::Handle<reco::PhotonCollection> photonH;
  for(unsigned int iPhoC=0; iPhoC<photonCollectionTags_.size(); iPhoC++){
    std::string& collectionTag(photonCollectionTags_[iPhoC]);

    _event.getByLabel(edm::InputTag(collectionTag), photonH);
    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillPhotons: size of PhotonCollection " << collectionTag << " = " << photonH->size();

    std::vector<std::string>& idCollectionTags(photonIdCollectionTags_[collectionTag]);
    unsigned nPhoIdC(idCollectionTags.size());
    std::vector<edm::ValueMap<bool> const*> phoIds;
    for(unsigned j(0); j != nPhoIdC; ++j){
      edm::Handle<edm::ValueMap<bool> > phoIdCH;
      _event.getByLabel(edm::InputTag(idCollectionTags[j]), phoIdCH);
      phoIds.push_back(phoIdCH.product());
    }

    std::vector<std::string>& isoDepTags(photonIsoDepTags_[collectionTag]);
    edm::ValueMap<double> const* chIsoMap(0);
    edm::ValueMap<double> const* nhIsoMap(0);
    edm::ValueMap<double> const* phIsoMap(0);
    if(isoDepTags.size() != 0){
      edm::Handle<edm::ValueMap<double> > chIsoH;
      edm::Handle<edm::ValueMap<double> > nhIsoH;
      edm::Handle<edm::ValueMap<double> > phIsoH;

      _event.getByLabel(edm::InputTag(isoDepTags[susy::kChargedHadron]), chIsoH);
      _event.getByLabel(edm::InputTag(isoDepTags[susy::kNeutralHadron]), nhIsoH);
      _event.getByLabel(edm::InputTag(isoDepTags[susy::kPhoton]), phIsoH);

      chIsoMap = chIsoH.product();
      nhIsoMap = nhIsoH.product();
      phIsoMap = phIsoH.product();
    }

    susy::PhotonCollection& susyCollection(susyEvent_->photons[collectionTag]);

    int ipho = 0;
    for(reco::PhotonCollection::const_iterator it = photonH->begin(); it != photonH->end(); ++it, ++ipho){

      if(it->pt() < photonThreshold_) continue;

      reco::PhotonRef phoRef(photonH, ipho);

      susy::Photon pho;

      bool isPF(it->photonCore()->isPFlowPhoton());

      // pack fiducial bits
      pho.fidBit |= (it->isEB()          << 0);
      pho.fidBit |= (it->isEE()          << 1);
      pho.fidBit |= (it->isEBEtaGap()    << 2);
      pho.fidBit |= (it->isEBPhiGap()    << 3);
      pho.fidBit |= (it->isEERingGap()   << 4);
      pho.fidBit |= (it->isEEDeeGap()    << 5);
      pho.fidBit |= (it->isEBEEGap()     << 6);
      pho.fidBit |= (isPF                << 7);

      pho.hadronicOverEm                    = it->hadronicOverEm();
      pho.hadronicDepth1OverEm              = it->hadronicDepth1OverEm();
      pho.hadronicDepth2OverEm              = it->hadronicDepth2OverEm();
      pho.hadTowOverEm                      = it->hadTowOverEm();

      pho.e1x5                              = it->e1x5();
      pho.e2x5                              = it->e2x5();
      pho.e3x3                              = it->e3x3();
      pho.e5x5                              = it->e5x5();
      pho.maxEnergyXtal                     = it->maxEnergyXtal();
      pho.sigmaEtaEta                       = it->sigmaEtaEta();
      pho.sigmaIetaIeta                     = it->sigmaIetaIeta();
      pho.r9                                = it->r9();

      pho.ecalRecHitSumEtConeDR04           = it->ecalRecHitSumEtConeDR04();
      pho.hcalDepth1TowerSumEtConeDR04      = it->hcalDepth1TowerSumEtConeDR04();
      pho.hcalDepth2TowerSumEtConeDR04      = it->hcalDepth2TowerSumEtConeDR04();

      pho.trkSumPtSolidConeDR04             = it->trkSumPtSolidConeDR04();
      pho.trkSumPtHollowConeDR04            = it->trkSumPtHollowConeDR04();
      pho.nTrkSolidConeDR04                 = it->nTrkSolidConeDR04();
      pho.nTrkHollowConeDR04                = it->nTrkHollowConeDR04();

      pho.ecalRecHitSumEtConeDR03           = it->ecalRecHitSumEtConeDR03();
      pho.hcalDepth1TowerSumEtConeDR03      = it->hcalDepth1TowerSumEtConeDR03();
      pho.hcalDepth2TowerSumEtConeDR03      = it->hcalDepth2TowerSumEtConeDR03();
      pho.trkSumPtSolidConeDR03             = it->trkSumPtSolidConeDR03();
      pho.trkSumPtHollowConeDR03            = it->trkSumPtHollowConeDR03();
      pho.nTrkSolidConeDR03                 = it->nTrkSolidConeDR03();
      pho.nTrkHollowConeDR03                = it->nTrkHollowConeDR03();

      pho.mipChi2                           = it->mipChi2();
      pho.mipTotEnergy                      = it->mipTotEnergy();
      pho.mipSlope                          = it->mipSlope();
      pho.mipIntercept                      = it->mipIntercept();
      pho.mipNhitCone                       = it->mipNhitCone();
      pho.mipIsHalo                         = it->mipIsHalo();

      // InvalidReference implies a misconfiguration (xyIsoMaps should be zero if the IsoDeposit objects do not exist for the collection) - throw
      if(chIsoMap) pho.chargedHadronIsoDeposit = (*chIsoMap)[phoRef];
      if(nhIsoMap) pho.neutralHadronIsoDeposit = (*nhIsoMap)[phoRef];
      if(phIsoMap) pho.photonIsoDeposit        = (*phIsoMap)[phoRef];

      if(isPF){
        // InvalidReference here would imply a fatal problem of AOD - throw
        reco::PFCandidate const* pfPhoton((*pfTranslationH)[phoRef].get());
        isolator03_->fGetIsolation(pfPhoton, pfH.product(), primVtxRef, vtxH);
      }
      else
        isolator03_->fGetIsolation(&*it, pfH.product(), primVtxRef, vtxH);

      pho.chargedHadronIso                  = isolator03_->getIsolationCharged();
      pho.neutralHadronIso                  = isolator03_->getIsolationNeutral();
      pho.photonIso                         = isolator03_->getIsolationPhoton();

      // for timing
      DetId seedId(0);

      // Cluster informations
      // Photon without a SuperClusterRef implies a fatal error - throw if scRef is null
      reco::SuperClusterRef scRef(isPF ? it->pfSuperCluster() : it->superCluster());
      pho.hcalIsoConeDR04_2012 = it->hcalTowerSumEtConeDR04() + (it->hadronicOverEm() - it->hadTowOverEm()) * scRef->energy() / cosh(scRef->eta());
      pho.hcalIsoConeDR03_2012 = it->hcalTowerSumEtConeDR03() + (it->hadronicOverEm() - it->hadTowOverEm()) * scRef->energy() / cosh(scRef->eta());

      pho.superClusterPreshowerEnergy = scRef->preshowerEnergy();
      pho.superClusterPhiWidth        = scRef->phiWidth();
      pho.superClusterEtaWidth        = scRef->etaWidth();

      reco::CaloClusterPtr const& seedClusterPtr(scRef->seed());

      // Seed of a SuperCluster can be null; see e.g. RecoParticleFlow/PFProducer/plugins/PFPhotonTranslator.cc
      if(seedClusterPtr.isNonnull()){
        seedId = seedClusterPtr->seed();
        if(!seedId.null()){
          double e1000 = spr::eECALmatrix(seedId, barrelRecHitsHandle, endcapRecHitsHandle, caloGeometry, caloTopology, 1, 0, 0, 0);
          double e0100 = spr::eECALmatrix(seedId, barrelRecHitsHandle, endcapRecHitsHandle, caloGeometry, caloTopology, 0, 1, 0, 0);
          double e0010 = spr::eECALmatrix(seedId, barrelRecHitsHandle, endcapRecHitsHandle, caloGeometry, caloTopology, 0, 0, 1, 0);
          double e0001 = spr::eECALmatrix(seedId, barrelRecHitsHandle, endcapRecHitsHandle, caloGeometry, caloTopology, 0, 0, 0, 1);

          pho.e1x2 = std::max(std::max(e1000,e0100), std::max(e0010, e0001));

          std::vector<float> crysCov = EcalClusterTools::localCovariances(*seedClusterPtr, barrelRecHitsHandle.product(), caloTopology);
          pho.sigmaIphiIphi = std::sqrt(crysCov[2]);
        }
      }

      pho.superClusterIndex = fillSuperCluster(scRef);

      pho.nPixelSeeds                       = it->electronPixelSeeds().size();
      pho.passelectronveto = !ConversionTools::hasMatchedPromptElectron(it->superCluster(), hVetoElectrons, hVetoConversions, beamSpot);

      // conversion Id
      if(it->conversions().size() > 0
         && it->conversions()[0]->nTracks() == 2
         && it->conversions()[0]->conversionVertex().isValid()
         && !it->conversions()[0]->conversionVertex().isFake()) {

        pho.convInfo   = kTRUE;
        pho.convDist   = it->conversions()[0]->distOfMinimumApproach();
        pho.convDcot   = it->conversions()[0]->pairCotThetaSeparation();
        pho.convVtxChi2 = it->conversions()[0]->conversionVertex().chi2();
        pho.convVtxNdof = it->conversions()[0]->conversionVertex().ndof();
        pho.convVertex.SetXYZ(it->conversions()[0]->conversionVertex().x(),
                              it->conversions()[0]->conversionVertex().y(),
                              it->conversions()[0]->conversionVertex().z());
        pho.convDxy = it->conversions()[0]->dxy(beamSpot);
        pho.convDz  = it->conversions()[0]->dz(beamSpot);
        pho.convLxy = it->conversions()[0]->lxy(beamSpot);
        pho.convLz  = it->conversions()[0]->lz(beamSpot);
        pho.convZofPVFromTracks = it->conversions()[0]->zOfPrimaryVertexFromTracks(beamSpot);
        pho.convTrackChargeProd = (it->conversions()[0]->tracks())[0]->charge() * (it->conversions()[0]->tracks())[1]->charge();
        pho.convTrack1nHit = (it->conversions()[0]->tracks())[0]->found();
        pho.convTrack2nHit = (it->conversions()[0]->tracks())[1]->found();
        pho.convTrack1chi2 = (it->conversions()[0]->tracks())[0]->chi2();
        pho.convTrack1pT = (it->conversions()[0]->tracks())[0]->pt();
        pho.convTrack2chi2 = (it->conversions()[0]->tracks())[1]->chi2();
        pho.convTrack2pT = (it->conversions()[0]->tracks())[1]->pt();
        std::vector<math::XYZPointF> InnerPos = it->conversions()[0]->tracksInnerPosition();
        pho.convTrack1InnerZ = InnerPos[0].z();
        pho.convTrack2InnerZ = InnerPos[1].z();
        pho.convTrack1InnerX = InnerPos[0].x();
        pho.convTrack2InnerX = InnerPos[1].x();
        pho.convTrack1InnerY = InnerPos[0].y();
        pho.convTrack2InnerY = InnerPos[1].y();
        std::vector<double> signedd0 = it->conversions()[0]->tracksSigned_d0();
        pho.convTrack1Signedd0 = signedd0[0];
        pho.convTrack2Signedd0 = signedd0[1];

      }

      // Photon Id
      for(unsigned k(0); k != nPhoIdC; ++k){
        try{
          pho.idPairs[idCollectionTags[k]] = (*phoIds[k])[phoRef];
        }
        catch(cms::Exception& e){
          if(e.category() == "InvalidReference")
            edm::LogWarning("InvalidReference") << "Photon Id " << idCollectionTags[k] << " does not exist for collection " << collectionTag << " instance " << ipho;
          else
            throw;
        }
      }

      pho.caloPosition.SetXYZ(it->caloPosition().x(),it->caloPosition().y(),it->caloPosition().z());
      pho.vertex.SetXYZ(it->vx(),it->vy(),it->vz());
      pho.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());

      // store seed timing information
      if(seedId.subdetId() == EcalBarrel) {
        EcalRecHitCollection::const_iterator hit = barrelRecHitsHandle->find(seedId);
        if ((hit != barrelRecHitsHandle->end()) && hit->isTimeValid()) pho.seedTime = hit->time();
      }
      else if(seedId.subdetId() == EcalEndcap) {
        EcalRecHitCollection::const_iterator hit = endcapRecHitsHandle->find(seedId);
        if ((hit != endcapRecHitsHandle->end()) && hit->isTimeValid()) pho.seedTime = hit->time();
      }

      fillTracksAround(*it, 0.4, trkH);

      // using kt6PFJets:rho
      pho.MVAregEnergyAndErr = scEnergyCorrector_.CorrectedEnergyWithErrorV3(*it, *vtxWithBSH, susyEvent_->rho, lazyTools, _eventSetup);
      double regEnergy(pho.MVAregEnergyAndErr.first);
      pho.MVAcorrMomentum.SetVect(pho.caloPosition * (regEnergy / pho.caloPosition.Mag()));
      pho.MVAcorrMomentum.SetE(regEnergy);

      susyCollection.push_back(pho);

      if(debugLevel_ > 2) edm::LogInfo(name()) << "pt, e, hadEm : " << it->pt()
                                               << ", " << it->energy()
                                               << ", " << it->hadronicOverEm();

    }// for it
  }
}

void
SusyNtuplizer::fillElectrons(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(electronCollectionTags_.size() == 0) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillElectrons";

  math::XYZPoint beamSpot(susyEvent_->beamSpot.X(), susyEvent_->beamSpot.Y(), susyEvent_->beamSpot.Z());

  edm::Handle<reco::ConversionCollection> conversionsH;
  _event.getByLabel("allConversions", conversionsH);

  edm::Handle<reco::GsfElectronCollection> electronH;
  for(unsigned int iEleC=0; iEleC < electronCollectionTags_.size(); iEleC++){
    std::string& collectionTag(electronCollectionTags_[iEleC]);

    _event.getByLabel(edm::InputTag(collectionTag), electronH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillElectrons: size of ElectronCollection " << collectionTag << " = " << electronH->size();

    std::vector<std::string>& idCollectionTags(electronIdCollectionTags_[collectionTag]);
    unsigned nEleIdC(idCollectionTags.size());
    std::vector<edm::ValueMap<float> const*> eleIds;
    for(unsigned j(0); j != nEleIdC; ++j){
      edm::Handle<edm::ValueMap<float> > eleIdCH;
      _event.getByLabel(edm::InputTag(idCollectionTags[j]), eleIdCH);
      eleIds.push_back(eleIdCH.product());
    }

    std::vector<std::string>& isoDepTags(electronIsoDepTags_[collectionTag]);
    edm::ValueMap<double> const* chIsoMap(0);
    edm::ValueMap<double> const* nhIsoMap(0);
    edm::ValueMap<double> const* phIsoMap(0);
    if(isoDepTags.size() != 0){
      edm::Handle<edm::ValueMap<double> > chIsoH;
      edm::Handle<edm::ValueMap<double> > nhIsoH;
      edm::Handle<edm::ValueMap<double> > phIsoH;

      _event.getByLabel(edm::InputTag(isoDepTags[susy::kChargedHadron]), chIsoH);
      _event.getByLabel(edm::InputTag(isoDepTags[susy::kNeutralHadron]), nhIsoH);
      _event.getByLabel(edm::InputTag(isoDepTags[susy::kPhoton]), phIsoH);

      chIsoMap = chIsoH.product();
      nhIsoMap = nhIsoH.product();
      phIsoMap = phIsoH.product();
    }

    susy::ElectronCollection& susyCollection(susyEvent_->electrons[collectionTag]);

    int iele = 0;
    for(reco::GsfElectronCollection::const_iterator it = electronH->begin(); it != electronH->end(); ++it, ++iele){

      if(it->pt() < electronThreshold_) continue;

      reco::GsfElectronRef eleRef(electronH, iele);

      susy::Electron ele;

      bool isPF = (it->candidateP4Kind() == reco::GsfElectron::P4_PFLOW_COMBINATION);

      // fiducial bits
      ele.fidBit |= it->isEB()        << 0;
      ele.fidBit |= it->isEE()        << 1;
      ele.fidBit |= it->isEBEEGap()   << 2;
      ele.fidBit |= it->isEBEtaGap()  << 3;
      ele.fidBit |= it->isEBPhiGap()  << 4;
      ele.fidBit |= it->isEEDeeGap()  << 5;
      ele.fidBit |= it->isEERingGap() << 6;

      ele.scPixCharge = it->scPixCharge();

      ele.boolPack |= it->isGsfCtfScPixChargeConsistent() << 0;
      ele.boolPack |= it->isGsfScPixChargeConsistent()    << 1;
      ele.boolPack |= it->isGsfCtfChargeConsistent()      << 2;
      ele.boolPack |= it->ecalDrivenSeed()                << 3;
      ele.boolPack |= it->trackerDrivenSeed()             << 4;
      ele.boolPack |= it->passingCutBasedPreselection()   << 5;
      ele.boolPack |= it->passingMvaPreselection()        << 6;
      ele.boolPack |= it->ambiguous()                     << 7;
      ele.boolPack |= it->isEcalEnergyCorrected()         << 8;
      ele.boolPack |= it->isEnergyScaleCorrected()        << 9;
      ele.boolPack |= it->convFlags()                     << 10;
      ele.boolPack |= isPF                                << 11;

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

      ele.shFracInnerHits = it->shFracInnerHits();

      ele.vertex.SetXYZ(it->vx(),it->vy(),it->vz());

      if(isPF){
        ele.momentum.SetXYZT(it->p4(reco::GsfElectron::P4_PFLOW_COMBINATION).px(),
                             it->p4(reco::GsfElectron::P4_PFLOW_COMBINATION).py(),
                             it->p4(reco::GsfElectron::P4_PFLOW_COMBINATION).pz(),
                             it->p4(reco::GsfElectron::P4_PFLOW_COMBINATION).e());

        GsfElectron::ShowerShape const& ss(it->pfShowerShape());
        ele.sigmaEtaEta                       = ss.sigmaEtaEta;
        ele.sigmaIetaIeta                     = ss.sigmaIetaIeta;
        ele.sigmaIphiIphi                     = ss.sigmaIphiIphi;
        ele.e1x5                              = ss.e1x5;
        ele.e2x5Max                           = ss.e2x5Max;
        ele.e5x5                              = ss.e5x5;
        ele.r9                                = ss.r9;
        ele.hcalDepth1OverEcal                = ss.hcalDepth1OverEcal;
        ele.hcalDepth2OverEcal                = ss.hcalDepth2OverEcal;
      }
      else{
        ele.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());

        ele.sigmaEtaEta                       = it->sigmaEtaEta();
        ele.sigmaIetaIeta                     = it->sigmaIetaIeta();
        ele.sigmaIphiIphi                     = it->sigmaIphiIphi();
        ele.e1x5                              = it->e1x5();
        ele.e2x5Max                           = it->e2x5Max();
        ele.e5x5                              = it->e5x5();
        ele.r9                                = it->r9();
        ele.hcalDepth1OverEcal                = it->hcalDepth1OverEcal();
        ele.hcalDepth2OverEcal                = it->hcalDepth2OverEcal();
      }

      ele.dr03TkSumPt                = it->dr03TkSumPt();
      ele.dr03EcalRecHitSumEt        = it->dr03EcalRecHitSumEt();
      ele.dr03HcalDepth1TowerSumEt   = it->dr03HcalDepth1TowerSumEt();
      ele.dr03HcalDepth2TowerSumEt   = it->dr03HcalDepth2TowerSumEt();
      ele.dr03HcalDepth1TowerSumEtBc = it->dr03HcalDepth1TowerSumEtBc();
      ele.dr03HcalDepth2TowerSumEtBc = it->dr03HcalDepth2TowerSumEtBc();

      ele.dr04TkSumPt                = it->dr04TkSumPt();
      ele.dr04EcalRecHitSumEt        = it->dr04EcalRecHitSumEt();
      ele.dr04HcalDepth1TowerSumEt   = it->dr04HcalDepth1TowerSumEt();
      ele.dr04HcalDepth2TowerSumEt   = it->dr04HcalDepth2TowerSumEt();
      ele.dr04HcalDepth1TowerSumEtBc = it->dr04HcalDepth1TowerSumEtBc();
      ele.dr04HcalDepth2TowerSumEtBc = it->dr04HcalDepth2TowerSumEtBc();

      ele.hcalOverEcalBc             = it->hcalOverEcalBc();

      // InvalidReference implies a misconfiguration (xyIsoMaps should be zero if the IsoDeposit objects do not exist for the collection) - throw
      if(chIsoMap) ele.chargedHadronIso = (*chIsoMap)[eleRef];
      if(nhIsoMap) ele.neutralHadronIso = (*nhIsoMap)[eleRef];
      if(phIsoMap) ele.photonIso        = (*phIsoMap)[eleRef];

      ele.nMissingHits       = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      ele.passConversionVeto = !ConversionTools::hasMatchedConversion((*it), conversionsH, beamSpot);

      ele.convDist   = it->convDist();
      ele.convDcot   = it->convDcot();
      ele.convRadius = it->convRadius();

      if(isPF){
        ele.mvaStatus  = it->mvaOutput().status;
        ele.mva        = it->mvaOutput().mva;
      }

      //enum Classification { UNKNOWN=-1, GOLDEN=0, BIGBREM=1, OLDNARROW=2, SHOWERING=3, GAP=4 } ;
      ele.bremClass                         = char(it->classification());
      ele.fbrem                             = it->fbrem();

      ele.ecalEnergy                        = it->ecalEnergy();
      ele.ecalEnergyError                   = it->ecalEnergyError();
      ele.trackMomentumError                = it->trackMomentumError();

      ele.electronClusterIndex = fillCluster(it->electronCluster());
      ele.superClusterIndex    = isPF ? fillSuperCluster(it->pflowSuperCluster()) : fillSuperCluster(it->superCluster());

      ele.closestCtfTrackIndex = fillTrack(it->closestCtfTrackRef());
      ele.gsfTrackIndex        = fillGsfTrack(it->gsfTrack());

      if(ele.gsfTrackIndex != -1){
        susy::Track& track(susyEvent_->tracks[ele.gsfTrackIndex]);
        track.extrapolatedPositions["ECALInnerWall"].SetXYZ(it->trackPositionAtCalo().X(),it->trackPositionAtCalo().Y(),it->trackPositionAtCalo().Z());
      }

      for(unsigned k(0); k != nEleIdC; ++k){
        try{
          ele.idPairs[idCollectionTags[k]] = (*eleIds[k])[eleRef];
        }
        catch(cms::Exception& e){
          if(e.category() == "InvalidReference")
            edm::LogWarning("InvalidReference") << "Electron Id " << idCollectionTags[k] << " does not exist for collection " << collectionTag << " instance " << iele;
          else
            throw;
        }
      }

      susyCollection.push_back(ele);

      if(debugLevel_ > 2) edm::LogInfo(name()) << "pt, e, hadEm : " << it->pt()
                                               << ", " << it->energy()
                                               << ", " << it->hadronicOverEm();
    }// for it

  }
}

void
SusyNtuplizer::fillMuons(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(muonCollectionTags_.size() == 0) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillMuons";

  edm::Handle<reco::MuonCollection> muonH;
  for(unsigned iMuCol(0); iMuCol < muonCollectionTags_.size(); ++iMuCol){
    std::string& collectionTag(muonCollectionTags_[iMuCol]);

    _event.getByLabel(edm::InputTag(collectionTag), muonH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMuons: size of MuonCollection " << collectionTag << " = " << muonH->size();

    std::vector<std::string>& idCollectionTags(muonIdCollectionTags_[collectionTag]);
    unsigned nMuIdC(idCollectionTags.size());
    std::vector<edm::ValueMap<bool> const*> muIds;
    for(unsigned i(0); i < nMuIdC; ++i){
      edm::Handle<edm::ValueMap<bool> > muIdCH;
      _event.getByLabel(edm::InputTag(idCollectionTags[i]), muIdCH);
      muIds.push_back(muIdCH.product());
    }

    susy::MuonCollection& susyCollection(susyEvent_->muons[collectionTag]);

    int imu = 0;
    for(reco::MuonCollection::const_iterator it = muonH->begin(); it != muonH->end(); ++it, ++imu){

      if(it->pt() < muonThreshold_) continue;

      reco::MuonRef muRef(muonH, imu);

      susy::Muon mu;

      mu.type                               = it->type();
      mu.bestTrackType                      = it->muonBestTrackType();
      mu.caloCompatibility                  = it->caloCompatibility();
      mu.nMatches                           = it->numberOfMatches();
      mu.nChambers                          = it->numberOfChambers();
      mu.nMatchedStations                   = it->numberOfMatchedStations();
      if(it->combinedMuon().isNonnull()){
        reco::TrackRef combRef(it->combinedMuon());
        mu.nValidHits                         = combRef->hitPattern().numberOfValidHits();
        mu.nValidTrackerHits                  = combRef->hitPattern().numberOfValidTrackerHits();
        mu.nValidMuonHits                     = combRef->hitPattern().numberOfValidMuonHits();
        mu.nPixelLayersWithMeasurement        = combRef->hitPattern().pixelLayersWithMeasurement();
        mu.nStripLayersWithMeasurement        = combRef->hitPattern().stripLayersWithMeasurement();
      }

      reco::MuonEnergy muEnergy(it->calEnergy());
      mu.emEnergy                           = muEnergy.em;
      mu.hadEnergy                          = muEnergy.had;

      mu.trackIsoR03                        = it->isolationR03().sumPt;
      mu.ecalIsoR03                         = it->isolationR03().emEt;
      mu.hcalIsoR03                         = it->isolationR03().hadEt;
      mu.trackIsoR05                        = it->isolationR05().sumPt;
      mu.ecalIsoR05                         = it->isolationR05().emEt;
      mu.hcalIsoR05                         = it->isolationR05().hadEt;

      mu.sumChargedHadronPt03               = it->pfIsolationR03().sumChargedHadronPt;
      mu.sumChargedParticlePt03             = it->pfIsolationR03().sumChargedParticlePt;
      mu.sumNeutralHadronEt03               = it->pfIsolationR03().sumNeutralHadronEt;
      mu.sumPhotonEt03                      = it->pfIsolationR03().sumPhotonEt;
      mu.sumNeutralHadronEtHighThreshold03  = it->pfIsolationR03().sumNeutralHadronEtHighThreshold;
      mu.sumPhotonEtHighThreshold03         = it->pfIsolationR03().sumPhotonEtHighThreshold;
      mu.sumPUPt03                          = it->pfIsolationR03().sumPUPt;

      mu.sumChargedHadronPt04               = it->pfIsolationR04().sumChargedHadronPt;
      mu.sumChargedParticlePt04             = it->pfIsolationR04().sumChargedParticlePt;
      mu.sumNeutralHadronEt04               = it->pfIsolationR04().sumNeutralHadronEt;
      mu.sumPhotonEt04                      = it->pfIsolationR04().sumPhotonEt;
      mu.sumNeutralHadronEtHighThreshold04  = it->pfIsolationR04().sumNeutralHadronEtHighThreshold;
      mu.sumPhotonEtHighThreshold04         = it->pfIsolationR04().sumPhotonEtHighThreshold;
      mu.sumPUPt04                          = it->pfIsolationR04().sumPUPt;

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
      }

      Short_t* trackIndices[] = {
        &mu.trackIndex, &mu.standAloneTrackIndex, &mu.combinedTrackIndex,
        &mu.tpfmsTrackIndex, &mu.pickyTrackIndex, &mu.dytTrackIndex
      };

      reco::Muon::MuonTrackType trackTypes[] = {
        reco::Muon::InnerTrack, reco::Muon::OuterTrack, reco::Muon::CombinedTrack,
        reco::Muon::TPFMS, reco::Muon::Picky, reco::Muon::DYT
      };

      for(unsigned iT(0); iT != sizeof(trackTypes) / sizeof(reco::Muon::MuonTrackType); ++iT){
        reco::Muon::MuonTrackType trackType(trackTypes[iT]);

        if(it->isAValidMuonTrack(trackType))
          *trackIndices[iT] = fillTrack(it->muonTrack(trackType));
      }

      mu.momentum.SetXYZT(it->p4().px(),it->p4().py(),it->p4().pz(),it->p4().e());

      for(unsigned k(0); k < nMuIdC; ++k){
        try{
          mu.idPairs[idCollectionTags[k]] = (*muIds[k])[muRef];
        }
        catch(cms::Exception& e){
          if(e.category() == "InvalidReference")
            edm::LogWarning("InvalidReference") << "Muon Id " << idCollectionTags[k] << " does not exist for collection " << collectionTag << " instance " << imu;
          else
            throw;
        }
      }

      susyCollection.push_back(mu);

      if(debugLevel_ > 2) edm::LogInfo(name()) << "type, emE, hadE, pt : " << it->type()
                                               << ", " << it->calEnergy().em
                                               << ", " << it->calEnergy().had
                                               << ", " << it->pt();
    }
  }

}

void
SusyNtuplizer::fillCaloJets(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(caloJetCollectionTags_.size() == 0) return;

  // JEC naming scheme in JetMETCorrections/Configuration/python/DefaultJEC_cff.py
  //
  // For CaloJets
  //
  // ak5CaloJetsL2L3
  // ak5CaloJetsL1FastL2L3

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillCaloJets";

  edm::Handle<reco::CaloJetCollection> jetH;
  for(unsigned iJetC=0; iJetC < caloJetCollectionTags_.size(); iJetC++){
    std::string& collectionTag(caloJetCollectionTags_[iJetC]);

    _event.getByLabel(edm::InputTag(collectionTag), jetH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillCaloJets: size of " << collectionTag << " JetCollection = " << jetH->size();

    std::string key(collectionTag);
    std::string toReplace("CaloJets");
    key.replace(key.find(toReplace), toReplace.size(), "");

    susy::CaloJetCollection& susyCollection(susyEvent_->caloJets[key]);

    edm::Handle<edm::ValueMap<reco::JetID> > jetIdH;
    _event.getByLabel(key + "JetID", jetIdH);

    JetCorrector const* corrL2L3(0);
    JetCorrector const* corrL1L2L3(0);
    if(_event.isRealData()){
      corrL2L3  = JetCorrector::getJetCorrector(key + "CaloL2L3Residual", _eventSetup);
      corrL1L2L3  = JetCorrector::getJetCorrector(key + "CaloL1L2L3Residual", _eventSetup);
    }
    else{
      corrL2L3  = JetCorrector::getJetCorrector(key + "CaloL2L3", _eventSetup);
      corrL1L2L3  = JetCorrector::getJetCorrector(key + "CaloL1L2L3", _eventSetup);
    }

    int ijet(0);
    for(reco::CaloJetCollection::const_iterator it = jetH->begin(); it != jetH->end(); ++it, ++ijet){

      TLorentzVector corrP4(it->px(), it->py(), it->pz(), it->energy());
      float l1l2l3Scale(corrL1L2L3->correction(*it, _event, _eventSetup));
      corrP4 *= l1l2l3Scale;

      if(corrP4.Pt() < jetThreshold_) continue;

      reco::CaloJetRef jetRef(jetH, ijet);

      susy::CaloJet jet;

      // Basic Jet
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

      // cf. RecoJets/JetProducers/plugins/VirtualJetProducer.cc (where the jets are actually produced)
      // Vertex correction is an available option for only the CaloJets, and relies on the parameter doPVCorrection (set in RecoJets/JetProducers/python/CaloJetParameters_cfi.py)
      // When doPVCorrection = True, jet momentum is calculated wrt the primary vertex (first element of the VertexCollection) and the position is saved as jet vertex
      // i.e. this is either just the primary vertex, or the origin (if doPVCorrection = False)
      jet.vertex.SetXYZ(it->vx(),it->vy(),it->vz());

      jet.momentum.SetXYZT(it->px(),it->py(),it->pz(),it->energy());
      jet.detectorP4.SetXYZT(it->detectorP4().px(),it->detectorP4().py(),
                             it->detectorP4().pz(),it->detectorP4().energy());

      jet.jecScaleFactors["L2L3"] = corrL2L3->correction(it->p4());
      jet.jecScaleFactors["L1L2L3"] = l1l2l3Scale;

      // accessing Jet ID information
      // InvalidReference here would imply a fatal problem in AOD - throw
      reco::JetID const& jetId((*jetIdH)[jetRef]);
      jet.fHPD                          = jetId.fHPD;
      jet.fRBX                          = jetId.fRBX;
      jet.n90Hits                       = jetId.n90Hits;
      jet.fSubDetector1                 = jetId.fSubDetector1;
      jet.fSubDetector2                 = jetId.fSubDetector2;
      jet.fSubDetector3                 = jetId.fSubDetector3;
      jet.fSubDetector4                 = jetId.fSubDetector4;
      jet.restrictedEMF                 = jetId.restrictedEMF;
      jet.nHCALTowers                   = jetId.nHCALTowers;
      jet.nECALTowers                   = jetId.nECALTowers;
      jet.approximatefHPD               = jetId.approximatefHPD;
      jet.approximatefRBX               = jetId.approximatefRBX;
      jet.hitsInN90                     = jetId.hitsInN90;
      jet.numberOfHits2RPC              = jetId.numberOfHits2RPC;
      jet.numberOfHits3RPC              = jetId.numberOfHits3RPC;
      jet.numberOfHitsRPC               = jetId.numberOfHitsRPC;

      susyCollection.push_back(jet);

      if(debugLevel_ > 2) edm::LogInfo(name()) << "pt, e : " << it->pt() << ", " << it->energy();

    }// for it
  }

}

void
SusyNtuplizer::fillPFJets(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(pfJetCollectionTags_.size() == 0) return;

  // JEC naming scheme in JetMETCorrections/Configuration/python/DefaultJEC_cff.py
  //
  // For PFJets
  //
  // ak5PFJetsL2L3
  // ak5PFJetsL1FastL2L3

  // If we are to have multiple PFJetCollections in the event, the following three getter blocks
  // should be defined for each jet collection independently a la JEC.
  // i.e. we need to prepend/append the collection name to the object labels and put the blocks
  // inside the loop over the collection tags.

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillPFJets";

  // Get b-tag information
  std::vector<edm::Handle<reco::JetFloatAssociation::Container> > jetDiscriminators;
  jetDiscriminators.resize(bTagCollectionTags_.size());
  for(size_t i = 0; i < bTagCollectionTags_.size(); i++){
    _event.getByLabel(edm::InputTag(bTagCollectionTags_[i]), jetDiscriminators[i]);
  }

  // Flavour matching for MC
  reco::JetFlavourMatchingCollection const* flavMatchAlg(0);
  reco::JetFlavourMatchingCollection const* flavMatchPhy(0);
  if(!_event.isRealData()){
    edm::Handle<reco::JetFlavourMatchingCollection> flavMatchH;

    _event.getByLabel("flavourAssociationAlg", flavMatchH);
    flavMatchAlg = flavMatchH.product();

    _event.getByLabel("flavourAssociationPhy", flavMatchH);
    flavMatchPhy = flavMatchH.product();
  }

  edm::Handle<edm::View<reco::PFJet> > jetH;
  for(unsigned iJetC=0; iJetC < pfJetCollectionTags_.size(); iJetC++) {
    std::string& collectionTag(pfJetCollectionTags_[iJetC]);

    _event.getByLabel(edm::InputTag(collectionTag), jetH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillPFJets: size of " << collectionTag << " JetCollection : " << jetH->size();

    std::string key(collectionTag);
    std::string toReplace("PFJets");
    key.replace(key.find(toReplace), toReplace.size(), "");

    JetCorrector const* corrL2L3(0);
    JetCorrector const* corrL1FastL2L3(0);
    if(_event.isRealData()){
      corrL2L3  = JetCorrector::getJetCorrector(key + "PFL2L3Residual", _eventSetup);
      corrL1FastL2L3  = JetCorrector::getJetCorrector(key + "PFL1FastL2L3Residual", _eventSetup);
    }
    else{
      corrL2L3  = JetCorrector::getJetCorrector(key + "PFL2L3", _eventSetup);
      corrL1FastL2L3  = JetCorrector::getJetCorrector(key + "PFL1FastL2L3", _eventSetup);
    }

    // Get Pileup Jet ID information
    std::vector<std::string>& puJetIdTags(puJetIdCollectionTags_[collectionTag]);
    std::vector<edm::Handle<edm::ValueMap<float> > > puJetIdMVACollections(puJetIdTags.size());
    std::vector<edm::Handle<edm::ValueMap<int> > > puJetIdFlagCollections(puJetIdTags.size());
    for(size_t i = 0; i < puJetIdTags.size(); i++){
      _event.getByLabel(edm::InputTag(puJetIdTags[i] + "Discriminant"), puJetIdMVACollections[i]);
      _event.getByLabel(edm::InputTag(puJetIdTags[i] + "Id"), puJetIdFlagCollections[i]);
    }

    susy::PFJetCollection& susyCollection(susyEvent_->pfJets[key]);

    int ijet(0);
    for(edm::View<reco::PFJet>::const_iterator it = jetH->begin(); it != jetH->end(); ++it, ++ijet){

      TLorentzVector corrP4(it->px(), it->py(), it->pz(), it->energy());
      float l1fastl2l3Scale(corrL1FastL2L3->correction(*it, _event, _eventSetup));
      corrP4 *= l1fastl2l3Scale;

      if(corrP4.Pt() < jetThreshold_) continue;

      edm::RefToBase<reco::Jet> jetRef(jetH->refAt(ijet));

      susy::PFJet jet;

      // Basic Jet
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

      // cf. RecoJets/JetProducers/plugins/VirtualJetProducer.cc (where the jets are actually produced)
      // Vertex correction is an available option for only the CaloJets. For all other jet types, vertex is always 0
      // See also the description above on CaloJet vertices
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

      jet.jecScaleFactors["L2L3"] = corrL2L3->correction(it->p4());
      jet.jecScaleFactors["L1FastL2L3"] = l1fastl2l3Scale;

      // add btag for this jet
      for(size_t k = 0; k < jetDiscriminators.size(); k++){
        try{
          float value = (*(jetDiscriminators[k]))[jetRef];
          jet.bTagDiscriminators.push_back(value);
        }
        catch(cms::Exception& e){
          if(e.category() == "InvalidReference")
            edm::LogWarning("InvalidReference") << "Btag discriminator " << bTagCollectionTags_[k] << " does not exist for collection " << collectionTag << " instance " << ijet;
          else
            throw;
        }
      }

      // add tracks from this jet into track collection and tally up the local id's.
      reco::TrackRefVector trkvec(it->getTrackRefs());
      for(reco::TrackRefVector::const_iterator t_it(trkvec.begin()); t_it != trkvec.end(); ++t_it)
        jet.tracklist.push_back(fillTrack(*t_it));

      std::vector<reco::PFCandidatePtr> constituents(it->getPFConstituents());
      for(unsigned iC(0); iC != constituents.size(); ++iC)
        jet.pfParticleList.push_back(fillPFParticle(constituents[iC]));

      // if MC, add parton flavor id matches
      if(!susyEvent_->isRealData){
        try{
          jet.algDefFlavour = (*flavMatchAlg)[jetRef].getFlavour();
          jet.phyDefFlavour = (*flavMatchPhy)[jetRef].getFlavour();
        }
        catch(cms::Exception& e){
          if(e.category() == "InvalidReference")
            edm::LogWarning("InvalidReference") << "Jet flavour matching does not exist for collection " << collectionTag << " instance " << ijet;
          else
            throw;
        }
      }

      // add puJetId variables
      for(size_t k = 0; k < puJetIdTags.size(); k++){
        try{
          float mva = (*puJetIdMVACollections[k])[jetRef];
          int idFlag = (*puJetIdFlagCollections[k])[jetRef];
          jet.puJetIdDiscriminants.push_back(mva);
          jet.puJetIdFlags.push_back(idFlag);
        }
        catch(cms::Exception& e){
          if(e.category() == "InvalidReference")
            edm::LogWarning("InvalidReference") << "PU jet Id does not exist for collection " << collectionTag << " instance " << ijet;
          else
            throw;
        }
      }

      susyCollection.push_back(jet);

      if(debugLevel_ > 2) edm::LogInfo(name()) << "pt, e : " << it->pt() << ", " << it->energy();

    }// for it

  }// for PFJet
}

void
SusyNtuplizer::fillJPTJets(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(jptJetCollectionTags_.size() == 0) return;

  // JPTJetCollection was commented out in November 2011.
  // Reviving the filler just for the sake of code clarity - no guarantee this function works under the current environment.

  // JEC naming scheme in JetMETCorrections/Configuration/python/DefaultJEC_cff.py
  //
  // For JPTJets
  //
  // ak5JPTJetsL2L3
  // ak5JPTJetsL1FastL2L3
  //

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillJPTJets";

  edm::Handle<reco::JPTJetCollection> jetH;
  for(unsigned iJetC=0; iJetC < jptJetCollectionTags_.size(); iJetC++){
    std::string& collectionTag(jptJetCollectionTags_[iJetC]);

    _event.getByLabel(edm::InputTag(collectionTag), jetH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillJPTJets: size of " << collectionTag << " JetCollection : " << jetH->size();

    std::string key(collectionTag);

    susy::JPTJetCollection& susyCollection(susyEvent_->jptJets[key]);

    // "Residual" corrections were commented out and no splitting was made wrt isRealData
    // Modifying assuming the usage follows the other jet collections (2013.3.29 Y.I.)
    // This part also has issues with multi-collection usage, but JPT jets are unused anyway..
    JetCorrector const* corrL2L3(0);
    JetCorrector const* corrL1L2L3(0);
    if(_event.isRealData()){
      corrL2L3 = JetCorrector::getJetCorrector("ak5JPTL2L3Residual", _eventSetup);
      corrL1L2L3 = JetCorrector::getJetCorrector("ak5JPTL1L2L3Residual", _eventSetup);
    }
    else{
      corrL2L3 = JetCorrector::getJetCorrector("ak5JPTL2L3", _eventSetup);
      corrL1L2L3  = JetCorrector::getJetCorrector("ak5JPTL1L2L3", _eventSetup);
    }

    int ijet = 0;
    for(reco::JPTJetCollection::const_iterator it = jetH->begin(); it != jetH->end(); ++it, ++ijet){

      TLorentzVector corrP4(it->px(), it->py(), it->pz(), it->energy());
      float l2l3Scale(corrL2L3->correction(it->p4()));

      corrP4 *= l2l3Scale;

      if(corrP4.Pt() < jetThreshold_) continue;

      reco::JPTJetRef jetRef(jetH,ijet++);

      susy::JPTJet jet;

      // Basic Jet
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

      // cf. RecoJets/JetProducers/plugins/VirtualJetProducer.cc (where the jets are actually produced)
      // Vertex correction is an available option for only the CaloJets. For all other jet types, vertex is always 0
      // See also the description above on CaloJet vertices
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

      jet.jecScaleFactors["L2L3"] = l2l3Scale;
      jet.jecScaleFactors["L1L2L3"] = corrL1L2L3->correction(*it, reco::JetBaseRef(jetRef), _event, _eventSetup);

      susyCollection.push_back(jet);

      if(debugLevel_ > 2) edm::LogInfo(name()) << "pt, e : " << it->pt() << ", " << it->energy();
    }// for it
  }

}

unsigned
SusyNtuplizer::fillTrack(reco::TrackRef const& _trkRef)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillTrack";

  if(_trkRef.isNull()) return -1;

  bool existed(false);
  unsigned index(fillTrackCommon(edm::refToPtr(_trkRef), existed));

  if(!existed){
    susy::Track& track(susyEvent_->tracks[index]);

    track.charge = _trkRef->charge();
    for(int i=0; i<reco::Track::dimension; i++) track.error[i] = _trkRef->error(i);
    track.momentum.SetXYZT(_trkRef->px(),_trkRef->py(),_trkRef->pz(),_trkRef->p());
  }

  return index;
}

unsigned
SusyNtuplizer::fillGsfTrack(reco::GsfTrackRef const& _trkRef)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillGsfTrack";

  if(_trkRef.isNull()) return -1;

  bool existed(false);
  unsigned index(fillTrackCommon(edm::Ptr<reco::Track>(edm::refToPtr(_trkRef)), existed));

  if(!existed){
    susy::Track& track(susyEvent_->tracks[index]);

    track.charge = _trkRef->chargeMode();
    for(int i=0; i<reco::GsfTrack::dimensionMode; i++) track.error[i] = _trkRef->errorMode(i);
    track.momentum.SetXYZT(_trkRef->pxMode(),_trkRef->pyMode(),_trkRef->pzMode(),_trkRef->pMode());
  }

  return index;
}

unsigned
SusyNtuplizer::fillSuperCluster(reco::SuperClusterRef const& _scRef)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillSuperCluster";

  if(susyEvent_->superClusters.size() != productStore_.superClusters.size())
    throw cms::Exception("RuntimeError") << "Number of buffered SCs does not match the number of SCs in the susyEvent";

  if(_scRef.isNull()) return -1;

  std::pair<SuperClusterStore::iterator, bool> insertion(productStore_.superClusters.insert(std::pair<reco::SuperClusterRef, unsigned>(_scRef, productStore_.superClusters.size())));

  if(!insertion.second) return insertion.first->second;
  else{
    susy::SuperCluster sc;

    sc.energy = _scRef->energy();
    sc.preshowerEnergy = _scRef->preshowerEnergy();
    sc.phiWidth = _scRef->phiWidth();
    sc.etaWidth = _scRef->etaWidth();
    sc.position.SetXYZ(_scRef->x(),_scRef->y(),_scRef->z());

    for(reco::CaloClusterPtrVector::const_iterator it = _scRef->clustersBegin(); it != _scRef->clustersEnd(); it++){
      unsigned index(fillCluster(*it));
      sc.basicClusterIndices.push_back(index);
      if(_scRef->seed() == *it) sc.seedClusterIndex = index;
    }

    susyEvent_->superClusters.push_back(sc);

    return susyEvent_->superClusters.size() - 1;
  }
}

unsigned
SusyNtuplizer::fillCluster(reco::CaloClusterPtr const& _clPtr)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillCluster";

  if(susyEvent_->clusters.size() != productStore_.basicClusters.size())
    throw cms::Exception("RuntimeError") << "Number of buffered SCs does not match the number of SCs in the susyEvent";

  if(_clPtr.isNull()) return -1;

  std::pair<CaloClusterStore::iterator, bool> insertion(productStore_.basicClusters.insert(std::pair<reco::CaloClusterPtr, unsigned>(_clPtr, productStore_.basicClusters.size())));

  if(!insertion.second) return insertion.first->second;
  else{
    susy::Cluster cl;

    cl.energy = _clPtr->energy();
    cl.position.SetXYZ(_clPtr->x(),_clPtr->y(),_clPtr->z());
    cl.nCrystals = _clPtr->size();

    susyEvent_->clusters.push_back(cl);

    return susyEvent_->clusters.size() - 1;
  }
}

unsigned
SusyNtuplizer::fillPFParticle(reco::PFCandidatePtr const& _partPtr)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillPFParticle";

  // temporary measure (I believe pfParticles should be just a vector instead of map<TString, vector>)
  susy::PFParticleCollection& col(susyEvent_->pfParticles[pfCandidateCollectionTag_]);

  if(col.size() != productStore_.pfCandidates.size())
    throw cms::Exception("RuntimeError") << "Number of buffered PFCandidates does not match the number of pfParticles in the susyEvent";

  if(_partPtr.isNull()) return -1;

  std::pair<PFCandidateStore::iterator, bool> insertion(productStore_.pfCandidates.insert(std::pair<reco::PFCandidatePtr, unsigned>(_partPtr, productStore_.pfCandidates.size())));

  if(!insertion.second) return insertion.first->second;
  else{
    susy::PFParticle pf;

    pf.pdgId         = _partPtr->translateTypeToPdgId(_partPtr->particleId());
    pf.charge        = _partPtr->charge();
    pf.ecalEnergy    = _partPtr->ecalEnergy();
    pf.rawEcalEnergy = _partPtr->rawEcalEnergy();
    pf.hcalEnergy    = _partPtr->hcalEnergy();
    pf.rawHcalEnergy = _partPtr->rawHcalEnergy();
    pf.pS1Energy     = _partPtr->pS1Energy();
    pf.pS2Energy     = _partPtr->pS2Energy();

    pf.vertex.SetXYZ(_partPtr->vx(),_partPtr->vy(),_partPtr->vz());
    pf.positionAtECALEntrance.SetXYZ(_partPtr->positionAtECALEntrance().x(),_partPtr->positionAtECALEntrance().y(),_partPtr->positionAtECALEntrance().z());
    pf.momentum.SetXYZT(_partPtr->px(),_partPtr->py(),_partPtr->pz(),_partPtr->energy());

    col.push_back(pf);

    return col.size() - 1;
  }
}

unsigned
SusyNtuplizer::fillTrackCommon(edm::Ptr<reco::Track> const& _trkPtr, bool& _existed)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillTrackCommon";

  if(susyEvent_->tracks.size() != productStore_.tracks.size())
    throw cms::Exception("RuntimeError") << "Number of buffered tracks does not match the number of tracks in the susyEvent";

  std::pair<TrackStore::iterator, bool> insertion(productStore_.tracks.insert(TrackStore::value_type(_trkPtr, productStore_.tracks.size())));

  _existed = !insertion.second;

  if(_existed) return insertion.first->second;
  else{
    susy::Track track;

    track.algorithm = _trkPtr->algo();
    track.quality = _trkPtr->qualityMask();

    track.chi2 = _trkPtr->chi2();
    track.ndof = _trkPtr->ndof();

    track.numberOfValidHits        = _trkPtr->hitPattern().numberOfValidHits();
    track.numberOfValidTrackerHits = _trkPtr->hitPattern().numberOfValidTrackerHits();
    track.numberOfValidMuonHits    = _trkPtr->hitPattern().numberOfValidMuonHits();
    track.numberOfValidPixelHits   = _trkPtr->hitPattern().numberOfValidPixelHits();
    track.numberOfValidStripHits   = _trkPtr->hitPattern().numberOfValidStripHits();
    track.vertex.SetXYZ(_trkPtr->vx(),_trkPtr->vy(),_trkPtr->vz());

    if(recoMode_){
      reco::TransientTrack ttk(transientTrackBuilder_->build(*_trkPtr));
      TrajectoryStateOnSurface const tsos(ttk.innermostMeasurementState());
      if(tsos.isValid()){
        float eta(tsos.globalPosition().eta());
        std::map<TString, Surface const*> surfaces;

        if(std::abs(eta) < susy::etaGap){
          surfaces["ECALInnerWall"] = &PFGeometry::barrelBound(PFGeometry::ECALInnerWall);
          surfaces["HCALInnerWall"] = &PFGeometry::barrelBound(PFGeometry::HCALInnerWall);
        }
        else if(eta > susy::etaGap && eta < susy::etaMax){
          surfaces["ECALInnerWall"] = &PFGeometry::positiveEndcapDisk(PFGeometry::ECALInnerWall);
          surfaces["HCALInnerWall"] = &PFGeometry::positiveEndcapDisk(PFGeometry::HCALInnerWall);
          surfaces["PS1Wall"] = &PFGeometry::positiveEndcapDisk(PFGeometry::PS1Wall);
          surfaces["PS2Wall"] = &PFGeometry::positiveEndcapDisk(PFGeometry::PS2Wall);
        }
        else if(eta < -susy::etaGap && eta > -susy::etaMax){
          surfaces["ECALInnerWall"] = &PFGeometry::negativeEndcapDisk(PFGeometry::ECALInnerWall);
          surfaces["HCALInnerWall"] = &PFGeometry::negativeEndcapDisk(PFGeometry::HCALInnerWall);
          surfaces["PS1Wall"] = &PFGeometry::negativeEndcapDisk(PFGeometry::PS1Wall);
          surfaces["PS2Wall"] = &PFGeometry::negativeEndcapDisk(PFGeometry::PS2Wall);
        }

        TrajectoryStateOnSurface state;
        for(std::map<TString, Surface const*>::iterator sItr(surfaces.begin()); sItr != surfaces.end(); ++sItr){
          Plane const* plane(dynamic_cast<Plane const*>(sItr->second));
          if(plane)
            state = propagator_->propagate(tsos, *plane);
          else
            state = propagator_->propagate(tsos, *dynamic_cast<Cylinder const*>(sItr->second));

          TVector3& position(track.extrapolatedPositions[sItr->first]);
          if(state.isValid()) position.SetXYZ(state.globalPosition().x(), state.globalPosition().y(), state.globalPosition().z());
        }
      }
    }

    susyEvent_->tracks.push_back(track);

    return susyEvent_->tracks.size() - 1;
  }
}

void
SusyNtuplizer::fillTracksAround(reco::Candidate const& _cand, double _deltaR, edm::Handle<reco::TrackCollection> const& _tracks)
{
  if(debugLevel_ > 2) edm::LogInfo(name()) << "fillTracksAround";

  unsigned iTrk(0);
  reco::TrackCollection::const_iterator tEnd(_tracks->end());
  for(reco::TrackCollection::const_iterator tItr(_tracks->begin()); tItr != tEnd; ++tItr, ++iTrk){
    double dR(reco::deltaR(_cand, *tItr));
    if(dR > _deltaR) continue;
    
    reco::TrackRef ref(_tracks, iTrk);
    fillTrack(ref);
  }
}

void
SusyNtuplizer::finalize()
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "finalize";

  delete hltConfig_;
  hltConfig_ = 0;
  delete l1GtUtils_;
  l1GtUtils_ = 0;
  delete isolator03_;
  isolator03_ = 0;
  delete hcalLaser2012Filter_;
  hcalLaser2012Filter_ = 0;

  if(storeTriggerEvents_)
    triggerEvent_->write();

  delete triggerEvent_;
  triggerEvent_ = 0;

  if(susyTree_){
    TFile* outF(susyTree_->GetCurrentFile());

    outF->cd();
    susyTree_->Write();
    delete outF;
    susyTree_ = 0;
  }

  delete susyEvent_;
  susyEvent_ = 0;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SusyNtuplizer);
