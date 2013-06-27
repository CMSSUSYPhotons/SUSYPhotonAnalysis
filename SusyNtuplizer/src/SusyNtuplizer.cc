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
// $Id: SusyNtuplizer.cc,v 1.61 2013/05/19 09:34:47 yiiyama Exp $
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/Common/interface/RefToPtr.h"

#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"

#include "DataFormats/METReco/interface/MET.h"

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

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"

//for conversion safe electron veto
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "DataFormats/JetReco/interface/JPTJet.h"
#include "DataFormats/JetReco/interface/JPTJetCollection.h"

#include "DataFormats/JetReco/interface/JetID.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

// for ecal rechit related
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Calibration/IsolatedParticles/interface/eECALMatrix.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

// Jet Energy Correction
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// pileup summary info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// b-tagging info
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

// Pileup Jet Id
#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"

// PFIsolation
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"

// Photon SC energy MVA regression
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"

// MET filters
#include "PhysicsTools/Utilities/interface/EventFilterFromListStandAlone.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"

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
#include <TMVA/MsgLogger.h>

#include "SusyEvent.h"
#include "SusyTriggerEvent.h"

typedef std::vector<std::string> VString;

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
  VString muonCollectionTags_;
  VString electronCollectionTags_;
  VString photonCollectionTags_;
  VString caloJetCollectionTags_;
  VString pfJetCollectionTags_;
  VString jptJetCollectionTags_;
  VString metCollectionTags_;
  std::map<std::string, VString> photonIsoDepTags_;
  std::map<std::string, VString> electronIsoDepTags_;
  std::map<std::string, VString> bTagCollectionTags_;
  std::map<std::string, VString> qgTagCollectionTags_;
  std::map<std::string, VString> muonIdTags_;
  std::map<std::string, std::pair<std::string, std::string> > electronMVAIdTags_;
  std::map<std::string, std::pair<std::string, std::string> > jetFlavourMatchingTags_;
  std::map<std::string, VString> puJetIdCollectionTags_;
  std::string pfPUCandidatesTag_;
  std::string photonSCRegressionWeights_;
  std::map<unsigned, std::string> metFilterTags_;

  // for HLT prescales
  HLTConfigProvider* hltConfig_;

  // for L1 menu
  L1GtUtils* l1GtUtils_;

  // PFIsolator
  PFIsolationEstimator* isolator03_;

  // text-based event veto for HcalLaserEventList MET filter
  EventFilterFromListStandAlone* hcalLaserEventListFilter_;

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

  // flag for whether or not MC is fastsim
  // certain collections are not produced by FamosSequence so we skip them
  // default: false (turned on in runOverAOD.py)
  bool isFastSim_;

  // flag for recording TriggerEvent into a separate tree (susyTriggers)
  bool storeTriggerEvents_;

  // flag for storing lumiSummary info -- simple fix for ReRecos missing this
  // default : true (turned off for 24Aug2012 and 13Jul2012 datasets in runOverAOD.py)
  bool storeLumiInfo_;

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
  muonCollectionTags_(iConfig.getParameter<VString>("muonCollectionTags")),
  electronCollectionTags_(iConfig.getParameter<VString>("electronCollectionTags")),
  photonCollectionTags_(iConfig.getParameter<VString>("photonCollectionTags")),
  caloJetCollectionTags_(iConfig.getParameter<VString>("caloJetCollectionTags")),
  pfJetCollectionTags_(iConfig.getParameter<VString>("pfJetCollectionTags")),
  jptJetCollectionTags_(iConfig.getParameter<VString>("jptJetCollectionTags")),
  metCollectionTags_(iConfig.getParameter<VString>("metCollectionTags")),
  photonIsoDepTags_(),
  electronIsoDepTags_(),
  bTagCollectionTags_(),
  qgTagCollectionTags_(),
  muonIdTags_(),
  electronMVAIdTags_(),
  jetFlavourMatchingTags_(),
  puJetIdCollectionTags_(),
  pfPUCandidatesTag_(iConfig.getParameter<std::string>("pfPUCandidatesTag")),
  photonSCRegressionWeights_(iConfig.getParameter<edm::FileInPath>("photonSCRegressionWeights").fullPath()),
  metFilterTags_(),
  hltConfig_(0),
  l1GtUtils_(0),
  isolator03_(0),
  hcalLaserEventListFilter_(0),
  debugLevel_(iConfig.getParameter<int>("debugLevel")),
  storeL1Info_(iConfig.getParameter<bool>("storeL1Info")),
  storeHLTInfo_(iConfig.getParameter<bool>("storeHLTInfo")),
  storeGenInfo_(iConfig.getParameter<bool>("storeGenInfo")),
  storeGeneralTracks_(iConfig.getParameter<bool>("storeGeneralTracks")),
  storePFJetPartonMatches_(iConfig.getParameter<bool>("storePFJetPartonMatches")),
  isFastSim_(iConfig.getParameter<bool>("isFastSim")),
  storeTriggerEvents_(iConfig.getParameter<bool>("storeTriggerEvents")),
  storeLumiInfo_(iConfig.getParameter<bool>("storeLumiInfo")),
  electronThreshold_(iConfig.getParameter<double>("electronThreshold")),
  muonThreshold_(iConfig.getParameter<double>("muonThreshold")),
  photonThreshold_(iConfig.getParameter<double>("photonThreshold")),
  jetThreshold_(iConfig.getParameter<double>("jetThreshold")),
  pfParticleThreshold_(iConfig.getParameter<double>("pfParticleThreshold")),
  scEnergyCorrector_(),
  productStore_(),
  susyEvent_(0),
  susyTree_(0),
  triggerEvent_(0)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "ctor";

  /*
    Suppress TMVA messages. Only effective for MVA methods invoked from this module.
   */
  TMVA::MsgLogger::InhibitOutput();

  /*
    Check the input file to photon supercluster energy correction.
   */
  TFile* dummyFile(TFile::Open(photonSCRegressionWeights_.c_str()));
  if(!dummyFile || dummyFile->IsZombie())
    throw cms::Exception("IOError") << "Photon SC MVA regression weight file " << photonSCRegressionWeights_ << " cannot be opened";
  delete dummyFile;

  /*
    Get tags for isolation values calculated through isoDeposit, for Photon and Electron collections.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> VString of tags].
   */
  edm::ParameterSet const& photonIsoDepTags(iConfig.getParameterSet("photonIsoDepTags"));
  for(VString::iterator tItr(photonCollectionTags_.begin()); tItr != photonCollectionTags_.end(); ++tItr){
    if(!photonIsoDepTags.existsAs<edm::ParameterSet>(*tItr)) continue;
    edm::ParameterSet const& tagsPSet(photonIsoDepTags.getParameterSet(*tItr));

    VString& tags(photonIsoDepTags_[*tItr]);
    tags.resize(3);

    tags[0] = tagsPSet.getParameter<std::string>("chargedHadron");
    tags[1] = tagsPSet.getParameter<std::string>("neutralHadron");
    tags[2] = tagsPSet.getParameter<std::string>("photon");
  }
  edm::ParameterSet const& electronIsoDepTags(iConfig.getParameterSet("electronIsoDepTags"));
  for(VString::iterator tItr(electronCollectionTags_.begin()); tItr != electronCollectionTags_.end(); ++tItr){
    if(!electronIsoDepTags.existsAs<edm::ParameterSet>(*tItr)) continue;
    edm::ParameterSet const& tagsPSet(electronIsoDepTags.getParameterSet(*tItr));

    VString& tags(electronIsoDepTags_[*tItr]);
    tags.resize(3);

    tags[0] = tagsPSet.getParameter<std::string>("chargedHadron");
    tags[1] = tagsPSet.getParameter<std::string>("neutralHadron");
    tags[2] = tagsPSet.getParameter<std::string>("photon");
  }

  /*
    Get tags for muon ID.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> vector tags].
  */
  edm::ParameterSet const& muonIdTags(iConfig.getParameterSet("muonIdTags"));
  for(VString::iterator tItr(muonCollectionTags_.begin()); tItr != muonCollectionTags_.end(); ++tItr){
    if(!muonIdTags.existsAs<edm::ParameterSet>(*tItr)) continue;
    edm::ParameterSet const& tagsPSet(muonIdTags.getParameterSet(*tItr));

    VString& tags(muonIdTags_[*tItr]);
    tags.resize(6);

    tags[0] = tagsPSet.getParameter<std::string>("TMLastStationLoose");
    tags[1] = tagsPSet.getParameter<std::string>("TMLastStationTight");
    tags[2] = tagsPSet.getParameter<std::string>("TMOneStationLoose");
    tags[3] = tagsPSet.getParameter<std::string>("TMOneStationTight");
    tags[4] = tagsPSet.getParameter<std::string>("TMLastStationLowPtLoose");
    tags[5] = tagsPSet.getParameter<std::string>("TMLastStationLowPtTight");
  }

  /*
    Get tags for MVA-based electron ID.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> pair of tags].
  */
  edm::ParameterSet const& electronMVAIdTags(iConfig.getParameterSet("electronMVAIdTags"));
  for(VString::iterator tItr(electronCollectionTags_.begin()); tItr != electronCollectionTags_.end(); ++tItr){
    if(!electronMVAIdTags.existsAs<edm::ParameterSet>(*tItr)) continue;
    edm::ParameterSet const& tagsPSet(electronMVAIdTags.getParameterSet(*tItr));

    std::pair<std::string, std::string>& tagPair(electronMVAIdTags_[*tItr]);

    tagPair.first = tagsPSet.getParameter<std::string>("triggering");
    tagPair.second = tagsPSet.getParameter<std::string>("nonTriggering");
  }

  /*
    Get tags for b-tagging discriminators for PFJet collections.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> VString of tags].
  */
  edm::ParameterSet const& bTagCollectionTags(iConfig.getParameterSet("bTagCollectionTags"));
  for(VString::iterator jItr(pfJetCollectionTags_.begin()); jItr != pfJetCollectionTags_.end(); ++jItr){
    if(!bTagCollectionTags.existsAs<edm::ParameterSet>(*jItr)) continue;
    edm::ParameterSet const& tagsPSet(bTagCollectionTags.getParameterSet(*jItr));

    VString& tags(bTagCollectionTags_[*jItr]);
    tags.resize(susy::nBTagDiscriminators);

    tags[susy::kTCHE] = tagsPSet.getParameter<std::string>("TrackCountingHighEff");
    tags[susy::kTCHP] = tagsPSet.getParameter<std::string>("TrackCountingHighPur");
    tags[susy::kJP] = tagsPSet.getParameter<std::string>("JetProbability");
    tags[susy::kJBP] = tagsPSet.getParameter<std::string>("JetBProbability");
    tags[susy::kSSV] = tagsPSet.getParameter<std::string>("SimpleSecondaryVertex");
    tags[susy::kCSV] = tagsPSet.getParameter<std::string>("CombinedSecondaryVertex");
    tags[susy::kCSVMVA] = tagsPSet.getParameter<std::string>("CombinedSecondaryVertexMVA");
    tags[susy::kSE] = tagsPSet.getParameter<std::string>("SoftElectron");
    tags[susy::kSM] = tagsPSet.getParameter<std::string>("SoftMuon");
  }

  /*
    Get tags for quark/gluon tagging discriminators for PFJet collections.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> VString of tags].
  */
  edm::ParameterSet const& qgTagCollectionTags(iConfig.getParameterSet("qgTagCollectionTags"));
  for(VString::iterator jItr(pfJetCollectionTags_.begin()); jItr != pfJetCollectionTags_.end(); ++jItr){
    if(!qgTagCollectionTags.existsAs<std::string>(*jItr)) continue;
    std::string qgTagModule(qgTagCollectionTags.getParameter<std::string>(*jItr));

    VString& tags(qgTagCollectionTags_[*jItr]);
    tags.resize(susy::nQGDiscriminators);

    tags[susy::kQuarkLikelihood] = qgTagModule + ":qgLikelihood";
    tags[susy::kGluonMLP] = qgTagModule + ":qgMLP";
  }

  /*
    Get tags for MC jet-parton matching for PFJet collections.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> pair of tags].
  */
  edm::ParameterSet const& jetFlavourMatchingTags(iConfig.getParameterSet("jetFlavourMatchingTags"));
  for(VString::iterator jItr(pfJetCollectionTags_.begin()); jItr != pfJetCollectionTags_.end(); ++jItr){
    if(!jetFlavourMatchingTags.existsAs<edm::ParameterSet>(*jItr)) continue;
    edm::ParameterSet const& tagsPSet(jetFlavourMatchingTags.getParameterSet(*jItr));

    std::pair<std::string, std::string>& tagPair(jetFlavourMatchingTags_[*jItr]);

    tagPair.first = tagsPSet.getParameter<std::string>("alg");
    tagPair.second = tagsPSet.getParameter<std::string>("phy");
  }

  /*
    Get tags for pileup jet ID for PFJet collections.
    Since the values are only defined with regard to specific collections,
    the configuration is also organized into sub-configs for each collection.
    The tags are stored as a map [collection name -> VString of tags].
  */
  edm::ParameterSet const& puJetIdConfigs(iConfig.getParameterSet("puJetId"));
  for(VString::iterator jItr(pfJetCollectionTags_.begin()); jItr != pfJetCollectionTags_.end(); ++jItr){
    if(!puJetIdConfigs.existsAs<edm::ParameterSet>(*jItr)) continue;
    edm::ParameterSet const& idConfig(puJetIdConfigs.getParameterSet(*jItr));
    std::string tag(idConfig.getParameter<std::string>("tag"));
    VString algos(idConfig.getParameter<VString>("algorithms"));

    VString& collectionTags(puJetIdCollectionTags_[*jItr]);
    collectionTags.resize(susy::nPUJetIdAlgorithms);

    for(unsigned iA(0); iA != algos.size(); ++iA){
      if(algos[iA] == "full") collectionTags[susy::kPUJetIdFull] = tag + ":full";
      else if(algos[iA] == "cutbased") collectionTags[susy::kPUJetIdCutBased] = tag + ":cutbased";
      else if(algos[iA] == "simple") collectionTags[susy::kPUJetIdSimple] = tag + ":simple";
    }
  }

  /*
    Open the output file.
  */
  TString outputFileName(iConfig.getParameter<std::string>("outputFileName"));

  TFile* outF(TFile::Open(outputFileName, "RECREATE"));
  if(!outF || outF->IsZombie())
    throw cms::Exception("IOError") << "Cannot create file " << outputFileName;

  outF->cd();    
  susyTree_ = new TTree("susyTree", "SUSY Event");
  susyTree_->SetAutoSave(10000000); // 10M events

  /*
    Construct an Event object and tell it what collections it should expect.
    The character ':' in the collection tags need to be replaced with '_'.
    Also for jet collections, we historically have the convention to strip off the
    "CaloJets" and "PFJets" from the collection name.
  */
  susyEvent_ = new susy::Event;

  for(VString::iterator tItr(metCollectionTags_.begin()); tItr != metCollectionTags_.end(); ++tItr)
    susyEvent_->metMap[TString(*tItr).ReplaceAll(":", "_")] = susy::MET();
  for(VString::iterator tItr(muonCollectionTags_.begin()); tItr != muonCollectionTags_.end(); ++tItr)
    susyEvent_->muons[TString(*tItr).ReplaceAll(":", "_")] = susy::MuonCollection();
  for(VString::iterator tItr(electronCollectionTags_.begin()); tItr != electronCollectionTags_.end(); ++tItr)
    susyEvent_->electrons[TString(*tItr).ReplaceAll(":", "_")] = susy::ElectronCollection();
  for(VString::iterator tItr(photonCollectionTags_.begin()); tItr != photonCollectionTags_.end(); ++tItr)
    susyEvent_->photons[TString(*tItr).ReplaceAll(":", "_")] = susy::PhotonCollection();
  for(VString::iterator tItr(caloJetCollectionTags_.begin()); tItr != caloJetCollectionTags_.end(); ++tItr){
    TString key(*tItr);
    key.ReplaceAll("CaloJets", "");
    susyEvent_->caloJets[key] = susy::CaloJetCollection();
  }
  for(VString::iterator tItr(pfJetCollectionTags_.begin()); tItr != pfJetCollectionTags_.end(); ++tItr){
    TString key(*tItr);
    key.ReplaceAll("PF", "");
    key.ReplaceAll("Jets", "");
    susyEvent_->pfJets[key] = susy::PFJetCollection();
  }
  for(VString::iterator tItr(jptJetCollectionTags_.begin()); tItr != jptJetCollectionTags_.end(); ++tItr)
    susyEvent_->jptJets[TString(*tItr).ReplaceAll(":", "_")] = susy::JPTJetCollection();

  /*
    Get the list of MET filter results to store.
    Only the fiters with "run" parameter set to True are stored.
    At the same time set the MET filter bitmask to the correct setting for the dataset.
    This mask will be used in Event::passMetFilters() function.
  */
  susyEvent_->metFilterMask = 0;

  edm::ParameterSet const& metFilters(iConfig.getParameterSet("metFilters"));
  VString metFilterNames(metFilters.getParameterNames());
  for(unsigned iF(0); iF != metFilterNames.size(); ++iF){
    edm::ParameterSet const& filterConfig(metFilters.getParameterSet(metFilterNames[iF]));
    if(!filterConfig.getParameter<bool>("run")) continue;

    unsigned bitPosition(-1);
    if(metFilterNames[iF] == "CSCBeamHalo") bitPosition = susy::kCSCBeamHalo;
    else if(metFilterNames[iF] == "HcalNoise") bitPosition = susy::kHcalNoise;
    else if(metFilterNames[iF] == "EcalDeadCellTP") bitPosition = susy::kEcalDeadCellTP;
    else if(metFilterNames[iF] == "EcalDeadCellBE") bitPosition = susy::kEcalDeadCellBE;
    else if(metFilterNames[iF] == "HcalLaserOccupancy") bitPosition = susy::kHcalLaserOccupancy;
    else if(metFilterNames[iF] == "TrackingFailure") bitPosition = susy::kTrackingFailure;
    else if(metFilterNames[iF] == "EEBadSC") bitPosition = susy::kEEBadSC;
    else if(metFilterNames[iF] == "HcalLaserEventList") bitPosition = susy::kHcalLaserEventList;
    else if(metFilterNames[iF] == "EcalLaserCorr") bitPosition = susy::kEcalLaserCorr;
    else if(metFilterNames[iF] == "ManyStripClus53X") bitPosition = susy::kManyStripClus53X;
    else if(metFilterNames[iF] == "TooManyStripClus53X") bitPosition = susy::kTooManyStripClus53X;
    else if(metFilterNames[iF] == "LogErrorTooManyClusters") bitPosition = susy::kLogErrorTooManyClusters;
    else if(metFilterNames[iF] == "EERingOfFire") bitPosition = susy::kEERingOfFire;
    else if(metFilterNames[iF] == "InconsistentMuon") bitPosition = susy::kInconsistentMuon;
    else if(metFilterNames[iF] == "GreedyMuon") bitPosition = susy::kGreedyMuon;
    else if(metFilterNames[iF] == "HcalLaserRECOUserStep") bitPosition = susy::kHcalLaserRECOUserStep;

    metFilterTags_[bitPosition] = filterConfig.getParameter<std::string>("tag");
    if(filterConfig.getParameter<bool>("default")) susyEvent_->metFilterMask |= (1 << bitPosition);
  }

  /*
    Get the list of gridParams to be stored.
    By its nature, any addition of the gridParams will require modifying the function fillGenInfo and
    recompiling the code. In addition, the use must specify, as a VString in the python configuration,
    what the names of the parameters are.
  */
  VString gridParamsList(iConfig.getParameter<VString>("gridParams"));
  for(VString::iterator pItr(gridParamsList.begin()); pItr != gridParamsList.end(); ++pItr)
    susyEvent_->gridParams[*pItr] = 0.;

  /*
    Pass the output tree to the Event object for branch bookings.
  */
  susyEvent_->addOutput(*susyTree_);

  /*
    Create a trigger tree if configured.
  */
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

  if(metFilterTags_.find(susy::kHcalLaserEventList) != metFilterTags_.end()){
    edm::FileInPath listPath("EventFilter/HcalRawToDigi/data/HCALLaser2012AllDatasets.txt.gz");
    gzFile dummyFP(gzopen(listPath.fullPath().c_str(), "r"));
    if(dummyFP != 0){
      gzclose(dummyFP);
      hcalLaserEventListFilter_ = new EventFilterFromListStandAlone(listPath.fullPath());
    }
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

  if(storeL1Info_){
    l1GtUtils_->getL1GtRunCache(iRun, _eventSetup, true, false); // use event setup, do not use L1TriggerMenuLite

    TString currentConfig(l1GtUtils_->l1TriggerMenu());
    currentConfig.ReplaceAll("/", "_");
    currentConfig.ReplaceAll(".", "p");
    if(currentConfig != susyEvent_->l1Map.getMenuName()){
      int err(0);
      L1GtTriggerMenu const* menu(l1GtUtils_->ptrL1TriggerMenuEventSetup(err));
      if(err != 0)
	throw cms::Exception("RuntimeError") << "L1GtUtils failed to return the trigger menu";

      std::vector<std::string> paths;
      AlgorithmMap const& l1GtAlgorithms(menu->gtAlgorithmMap());
      for(CItAlgo iAlgo(l1GtAlgorithms.begin()); iAlgo != l1GtAlgorithms.end(); ++iAlgo)
	paths.push_back(iAlgo->second.algoName());

      susyEvent_->l1Map.setMenu(currentConfig, paths);
    }
  }

  if(storeHLTInfo_){
    //intialize HLTConfigProvider
    bool menuChanged(false);
    if(!hltConfig_->init(iRun, _eventSetup, "HLT", menuChanged))
      throw cms::Exception("RuntimeError") << "HLTConfigProvider::init() returned non 0";

    if(menuChanged){
      edm::LogInfo(name()) << "beginRun: HLT configuration changed to " << hltConfig_->tableName();

      TString currentConfig(hltConfig_->tableName());
      currentConfig.ReplaceAll("/online/collisions/", "");
      currentConfig.ReplaceAll("/cdaq/physics/", "");
      currentConfig.ReplaceAll("/", "_");
      currentConfig.ReplaceAll(".", "p");

      susyEvent_->hltMap.setMenu(currentConfig, hltConfig_->triggerNames());
    }
  }

  if(photonCollectionTags_.size() != 0){
    // initialize Photon SC energy MVA regression
    // EventSetup is not used if the third argument is false
    scEnergyCorrector_.Initialize(_eventSetup, photonSCRegressionWeights_, false);
  }
}

// ------------ method called to for each event  ------------
// fill the tree variables
void
SusyNtuplizer::analyze(edm::Event const& _event, edm::EventSetup const& _eventSetup)
{
  if(debugLevel_ > 0) edm::LogInfo(name()) << "analyze";

  /*
    Framework will guarantee that the destructor will be called in case of exceptions.
    There is no need to try-catch in the individual filler functions.
  */

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

  if(storeLumiInfo_) fillLumiSummary(_event, _eventSetup);

  fillBeamSpot(_event, _eventSetup);

  fillRhos(_event, _eventSetup);

  fillTriggerMaps(_event, _eventSetup);

  fillGeneralTracks(_event, _eventSetup);

  fillMetFilters(_event, _eventSetup);

  fillMet(_event, _eventSetup);

  fillPhotons(_event, _eventSetup);

  fillElectrons(_event, _eventSetup);

  fillMuons(_event, _eventSetup);

  fillCaloJets(_event, _eventSetup);

  fillPFJets(_event, _eventSetup);

  //  fillJPTJets(_event, _eventSetup);

  fillPFParticles(_event, _eventSetup);

  fillVertices(_event, _eventSetup);

  fillGenInfo(_event, _eventSetup);

  fillGenParticles(_event, _eventSetup);

  fillTriggerEvent(_event, _eventSetup);

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

    int err(0);
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

      susyEvent_->l1Map.set(algoName, prescale, decisionBeforeMask);
    }

  } // L1 only if data or fullsim

  if(storeHLTInfo_){

    if(debugLevel_ > 0) edm::LogInfo(name()) << "fillTriggerMaps: HLT map";

    edm::Handle<edm::TriggerResults> hltH;
    _event.getByLabel(edm::InputTag("TriggerResults", "", "HLT"), hltH);

    unsigned nHlt(hltH->size());

    edm::TriggerNames const& hltTriggerNames(_event.triggerNames(*hltH));
    if(nHlt != hltTriggerNames.size())
      throw cms::Exception("RuntimeError") << "TriggerPathName size mismatches !!!";

    // loop over hlt paths
    for(unsigned i(0); i != nHlt; ++i){
      std::string const& path(hltTriggerNames.triggerName(i));

      int prescale(hltConfig_->prescaleValue(_event, _eventSetup, path));

      susyEvent_->hltMap.set(path, prescale, hltH->accept(i));

      if(debugLevel_ > 2) edm::LogInfo(name()) << path << " : " << hltH->accept(i);
    }

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

  VertexHigherPtSquared sumPt2Calculator;

  susyEvent_->vertices.resize(vertices.size());

  for(unsigned iV(0); iV != vertices.size(); ++iV){
    reco::Vertex const& recoVtx(vertices[iV]);
    susy::Vertex& susyVtx(susyEvent_->vertices[iV]);

    susyVtx.tracksSize = recoVtx.tracksSize();
    susyVtx.sumPt2     = sumPt2Calculator.sumPtSquared(recoVtx);
    susyVtx.chi2       = recoVtx.chi2();
    susyVtx.ndof       = recoVtx.ndof();
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

  std::set<unsigned> irregular;
  irregular.insert(susy::kCSCBeamHalo);
  irregular.insert(susy::kHcalLaserEventList);
  irregular.insert(susy::kHcalLaserRECOUserStep);

  std::vector<bool> pass(susy::nMetFilters, true);

  std::map<unsigned, std::string>::iterator fItr;

  edm::Handle<bool> boolH;
  for(unsigned iF(0); iF != susy::nMetFilters; ++iF){
    // skip filters that are not a simple bool in the edm::Event
    if(irregular.find(iF) != irregular.end()) continue;
    // skip filters disabled in the job configuration
    fItr = metFilterTags_.find(iF);
    if(fItr == metFilterTags_.end()) continue;

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: " << fItr->first;

    _event.getByLabel(edm::InputTag(fItr->second), boolH);

    switch(iF){
    case susy::kManyStripClus53X:
    case susy::kTooManyStripClus53X:
    case susy::kLogErrorTooManyClusters:
    case susy::kLogErrorTooManyTripletsPairs:
    case susy::kLogErrorTooManySeeds:
      pass[iF] = !*boolH;
      break;
    default:
      pass[iF] = *boolH;
    }
  }

  if((fItr = metFilterTags_.find(susy::kCSCBeamHalo)) != metFilterTags_.end()){
    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: " << fItr->first;

    edm::Handle<reco::BeamHaloSummary> beamHaloSummary;
    _event.getByLabel(edm::InputTag(fItr->second), beamHaloSummary);

    pass[susy::kCSCBeamHalo] = !(beamHaloSummary->CSCTightHaloId());
  }

  if(metFilterTags_.find(susy::kHcalLaserEventList) != metFilterTags_.end() && hcalLaserEventListFilter_){
    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: Hcal laser filter based on an event list (Run2012ABC)";

    pass[susy::kHcalLaserEventList] = hcalLaserEventListFilter_->filter(_event.id().run(), _event.luminosityBlock(), _event.id().event());
  }

  if(metFilterTags_.find(susy::kHcalLaserRECOUserStep) != metFilterTags_.end()){
    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMetFilters: Hcal laser filter run at RECO (Run2012 2013 rerecos)";

    edm::Handle<edm::TriggerResults> recoPathResultsH;
    _event.getByLabel(edm::InputTag("TriggerResults", "", "RECO"), recoPathResultsH);
    edm::TriggerNames const& recoPathNames(_event.triggerNames(*recoPathResultsH));
    unsigned iUserStep(recoPathNames.triggerIndex("user_step"));
    if(iUserStep == recoPathNames.size())
      throw cms::Exception("RuntimeError") << "user_step not found in edmTriggerResults_TriggerResults__RECO";

    pass[susy::kHcalLaserRECOUserStep] = recoPathResultsH->accept(iUserStep);
  }

  for(unsigned iF(0); iF != susy::nMetFilters; ++iF)
    if(!pass[iF]) susyEvent_->metFilterBit &= ~(1 << iF);
}

void
SusyNtuplizer::fillPFParticles(edm::Event const& _event, edm::EventSetup const&)
{
  if(pfCandidateCollectionTag_ == "") return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillPFParticles";

  edm::Handle<reco::PFCandidateCollection> pfH;
  _event.getByLabel(edm::InputTag(pfCandidateCollectionTag_), pfH);

  edm::Handle<reco::PFCandidateCollection> pfPUH;
  _event.getByLabel(edm::InputTag(pfPUCandidatesTag_), pfPUH);

  if(debugLevel_ > 1) edm::LogInfo(name()) << "fillPFParticles: size of PFCandidateCollection = " << pfH->size();

  unsigned long iPart(0);
  for(reco::PFCandidateCollection::const_iterator pItr(pfH->begin()); pItr != pfH->end(); ++pItr, ++iPart){
    if(pItr->pt() < pfParticleThreshold_) continue;
    fillPFParticle(reco::PFCandidatePtr(pfH, iPart));

    if(debugLevel_ > 2) edm::LogInfo(name()) << "e, px, py, pz = " << pItr->energy() << ", "
                                             << pItr->px() << ", " << pItr->py() << ", " << pItr->pz();
  }

  std::set<unsigned long> puKeys;

  for(reco::PFCandidateCollection::const_iterator puItr(pfPUH->begin()); puItr != pfPUH->end(); ++puItr)
    puKeys.insert(puItr->sourceCandidatePtr(0).key());

  for(PFCandidateStore::iterator sItr(productStore_.pfCandidates.begin()); sItr != productStore_.pfCandidates.end(); ++sItr)
    if(puKeys.find(sItr->first.key()) != puKeys.end()) susyEvent_->pfParticles[sItr->second].isPU = kTRUE;
}

void
SusyNtuplizer::fillGenInfo(edm::Event const& _event, edm::EventSetup const&)
{
  if(!storeGenInfo_ || _event.isRealData()) return;

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillGenInfo";

  /*
    Filling gridParams. The user has to modify this part for any additional parameter,
    and give its name in the VString of gridParams in the python configuration.
   */
  /*--------------*/
  edm::Handle<GenEventInfoProduct> GenEventInfoHandle;
  if(_event.getByLabel("generator", GenEventInfoHandle) && GenEventInfoHandle->binningValues().size() > 0)
    susyEvent_->gridParams["ptHat"] = GenEventInfoHandle->binningValues()[0];

  /*--------------*/

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

    susy::MET& susyMet(susyEvent_->metMap.find(TString(colName).ReplaceAll(":", "_"))->second);

    edm::Handle<edm::View<reco::MET> > metH;
    _event.getByLabel(edm::InputTag(colName), metH);

    if(metH.isValid()){

      if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMet: size of MET coll " << colName << " = " << metH->size();

      reco::MET const& recoMet(*(metH->begin()));

      susyMet.mEt.Set(recoMet.p4().px(),recoMet.p4().py());
      susyMet.sumEt = recoMet.sumEt();
      susyMet.significance = recoMet.significance();

      if(debugLevel_ > 2) edm::LogInfo(name()) << "met, metX, metY, sumEt, significance : "
                                               << susyMet.mEt.Mod() << ", " << susyMet.mEt.X() << ", " << susyMet.mEt.Y() << ", "
                                               << susyMet.sumEt << ", " << susyMet.significance;

    }
    else{
      // FIX FOR NoPUMet

      if(debugLevel_ > 1) edm::LogInfo(name()) << "fillMet: MET " << colName << " invalid";

      susyMet.mEt.Set(0., 0.);
      susyMet.sumEt = -1.;
      susyMet.significance = -1.;
    }
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
  _event.getByLabel(pfCandidateCollectionTag_, pfH);

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

    VString& isoDepTags(photonIsoDepTags_[collectionTag]);
    edm::ValueMap<double> const* chIsoMap(0);
    edm::ValueMap<double> const* nhIsoMap(0);
    edm::ValueMap<double> const* phIsoMap(0);
    if(isoDepTags.size() != 0){
      edm::Handle<edm::ValueMap<double> > chIsoH;
      edm::Handle<edm::ValueMap<double> > nhIsoH;
      edm::Handle<edm::ValueMap<double> > phIsoH;

      _event.getByLabel(edm::InputTag(isoDepTags[0]), chIsoH);
      _event.getByLabel(edm::InputTag(isoDepTags[1]), nhIsoH);
      _event.getByLabel(edm::InputTag(isoDepTags[2]), phIsoH);

      chIsoMap = chIsoH.product();
      nhIsoMap = nhIsoH.product();
      phIsoMap = phIsoH.product();
    }

    susy::PhotonCollection& susyCollection(susyEvent_->photons.find(TString(collectionTag).ReplaceAll(":", "_"))->second);

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

      // Now minor gymnastics for grabbing the "worst isolation from any other
      // vertex" -- excluding the PV since we already have that one...

      // This is only relevant for the charged had iso -- the others don't
      // care what vx they come from...

      // Assign things at first to the PV -- this is the default if there is
      // only one vx

      pho.worstOtherVtxChargedHadronIso=pho.chargedHadronIso;
      pho.worstOtherVtxChargedHadronIsoVtxIdx=0;


      if (vtxH->size()>1) {

        // dummy vars to hold the maxima as we loop...
        Float_t dumchIso(0);
        Int_t dumchIdx(0);

        for (unsigned int ivx=1;ivx<vtxH->size();++ivx) {
          reco::VertexRef tmpVxRef(vtxH,ivx);
          // this is going to look very familiar -- some degree of laziness here...
          if (isPF) {
            // InvalidReference here would imply a fatal problem of AOD - throw
            reco::PFCandidate const* pfPhoton((*pfTranslationH)[phoRef].get());
            isolator03_->fGetIsolation(pfPhoton, pfH.product(), tmpVxRef, vtxH);
          }  
          else {
            isolator03_->fGetIsolation(&*it, pfH.product(), tmpVxRef, vtxH);
          }

          Float_t thischIso=isolator03_->getIsolationCharged();


          if (thischIso>dumchIso) {
            dumchIso=thischIso;
            dumchIdx=ivx;
          }

        } // vx for loop

        // now we assign what we got...

        pho.worstOtherVtxChargedHadronIso=dumchIso;
        pho.worstOtherVtxChargedHadronIsoVtxIdx=dumchIdx;

      } // vx size > 1

      // done with other vx stuff...


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

      pho.nPixelSeeds = it->electronPixelSeeds().size();
      pho.passelectronveto = !ConversionTools::hasMatchedPromptElectron(it->superCluster(), hVetoElectrons, hVetoConversions, beamSpot);

      // conversion Id
      reco::ConversionRefVector conversions(it->conversions());
      if(conversions.size() != 0){
        reco::Conversion const& conversion(*conversions[0]);

        if(conversion.nTracks() == 2
           && conversion.conversionVertex().isValid()
           && !conversion.conversionVertex().isFake()) {
          reco::Track const& track0(*conversion.tracks().at(0));
          reco::Track const& track1(*conversion.tracks().at(1));

          pho.convInfo   = kTRUE;
          pho.convDist   = conversion.distOfMinimumApproach();
          pho.convDcot   = conversion.pairCotThetaSeparation();
          pho.convVtxChi2 = conversion.conversionVertex().chi2();
          pho.convVtxNdof = conversion.conversionVertex().ndof();
          pho.convVertex.SetXYZ(conversion.conversionVertex().x(),
                                conversion.conversionVertex().y(),
                                conversion.conversionVertex().z());
          pho.convDxy = conversion.dxy(beamSpot);
          pho.convDz  = conversion.dz(beamSpot);
          pho.convLxy = conversion.lxy(beamSpot);
          pho.convLz  = conversion.lz(beamSpot);
          pho.convZofPVFromTracks = conversion.zOfPrimaryVertexFromTracks(beamSpot);
          pho.convTrackChargeProd = track0.charge() * track1.charge();
          pho.convTrack1nHit = track0.found();
          pho.convTrack2nHit = track1.found();
          pho.convTrack1chi2 = track0.chi2();
          pho.convTrack1pT = track0.pt();
          pho.convTrack2chi2 = track1.chi2();
          pho.convTrack2pT = track1.pt();
          std::vector<math::XYZPointF> InnerPos = conversion.tracksInnerPosition();
          pho.convTrack1InnerZ = InnerPos[0].z();
          pho.convTrack2InnerZ = InnerPos[1].z();
          pho.convTrack1InnerX = InnerPos[0].x();
          pho.convTrack2InnerX = InnerPos[1].x();
          pho.convTrack1InnerY = InnerPos[0].y();
          pho.convTrack2InnerY = InnerPos[1].y();
          std::vector<double> signedd0 = conversion.tracksSigned_d0();
          pho.convTrack1Signedd0 = signedd0[0];
          pho.convTrack2Signedd0 = signedd0[1];
        }
      }

      pho.caloPosition.SetXYZ(it->caloPosition().x(),it->caloPosition().y(),it->caloPosition().z());
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
      //MVA supercluster regression corrected energy and error.
      //Photons use (0,0,0) as vertex, so for the best mass resolution they need to be corrected to the primary vertex as follows, where PhoOne and PhoTwo are susy::Photon* :
      /*
      TVector3 vPos = (event.vertices[0]).position;
      TVector3 dirPhoOne=PhoOne->caloPosition - vPos,dirPhoTwo=PhoTwo->caloPosition - vPos;
      TVector3 pOne=dirPhoOne.Unit()*PhoOne->MVAregEnergy,pTwo=dirPhoTwo.Unit()*PhoTwo->MVAregEnergy;
      TLorentzVector p4PhoOneVcorr(pOne.x(),pOne.y(),pOne.z(),PhoOne->MVAregEnergy),p4PhoTwoVcorr(pTwo.x(),pTwo.y(),pTwo.z(),PhoTwo->MVAregEnergy);
      float InvarMassMVAcorrVertexCorr=(p4PhoOneVcorr+p4PhoTwoVcorr).M();
      */
      std::pair<double,double> scEcor = scEnergyCorrector_.CorrectedEnergyWithErrorV3(*it, *vtxWithBSH, susyEvent_->rho, lazyTools, _eventSetup);
      pho.MVAregEnergy = scEcor.first;
      pho.MVAregErr = scEcor.second;

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

    VString& isoDepTags(electronIsoDepTags_[collectionTag]);
    edm::ValueMap<double> const* chIsoMap(0);
    edm::ValueMap<double> const* nhIsoMap(0);
    edm::ValueMap<double> const* phIsoMap(0);
    if(isoDepTags.size() != 0){
      edm::Handle<edm::ValueMap<double> > vmH;

      _event.getByLabel(edm::InputTag(isoDepTags[0]), vmH);
      chIsoMap = vmH.product();
      _event.getByLabel(edm::InputTag(isoDepTags[1]), vmH);
      nhIsoMap = vmH.product();
      _event.getByLabel(edm::InputTag(isoDepTags[2]), vmH);
      phIsoMap = vmH.product();
    }

    edm::ValueMap<float> const* mvaTrigMap(0);
    edm::ValueMap<float> const* mvaNonTrigMap(0);
    if(electronMVAIdTags_.find(collectionTag) != electronMVAIdTags_.end()){
      edm::Handle<edm::ValueMap<float> > vmH;
      std::pair<std::string, std::string>& tags(electronMVAIdTags_[collectionTag]);

      _event.getByLabel(edm::InputTag(tags.first), vmH);
      mvaTrigMap = vmH.product();
      _event.getByLabel(edm::InputTag(tags.second), vmH);
      mvaNonTrigMap = vmH.product();
    }

    susy::ElectronCollection& susyCollection(susyEvent_->electrons.find(TString(collectionTag).ReplaceAll(":", "_"))->second);

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
      ele.boolPack |= isPF                                << 10;

      ele.scPixCharge = it->scPixCharge();

      ele.convFlag = it->convFlags();

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

      ele.trackPositionAtVtx.SetXYZ(it->trackPositionAtVtx().X(),it->trackPositionAtVtx().Y(),it->trackPositionAtVtx().Z());
      ele.trackPositionAtCalo.SetXYZ(it->trackPositionAtCalo().X(),it->trackPositionAtCalo().Y(),it->trackPositionAtCalo().Z());
      ele.trackMomentumAtVtx.SetXYZT(it->trackMomentumAtVtx().X(),it->trackMomentumAtVtx().Y(),
                                     it->trackMomentumAtVtx().Z(),it->trackMomentumAtVtx().R());
      ele.trackMomentumAtCalo.SetXYZT(it->trackMomentumAtCalo().X(),it->trackMomentumAtCalo().Y(),
                                      it->trackMomentumAtCalo().Z(),it->trackMomentumAtCalo().R());
      ele.trackMomentumOut.SetXYZT(it->trackMomentumOut().X(),it->trackMomentumOut().Y(),
                                   it->trackMomentumOut().Z(),it->trackMomentumOut().R());
      ele.trackMomentumAtEleClus.SetXYZT(it->trackMomentumAtEleClus().X(),it->trackMomentumAtEleClus().Y(),
                                         it->trackMomentumAtEleClus().Z(),it->trackMomentumAtEleClus().R());
      ele.trackMomentumAtVtxWithConstraint.SetXYZT(it->trackMomentumAtVtxWithConstraint().X(),it->trackMomentumAtVtxWithConstraint().Y(),
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

      // InvalidReference implies a misconfiguration - throw
      if(mvaTrigMap) ele.mvaTrig    = (*mvaTrigMap)[eleRef];
      if(mvaNonTrigMap) ele.mvaNonTrig = (*mvaNonTrigMap)[eleRef];

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

    VString& idTags(muonIdTags_[collectionTag]);
    std::vector<edm::ValueMap<bool> const*> idMaps(6, NULL);
    if(idTags.size() != 0){
      edm::Handle<edm::ValueMap<bool> > vmH;

      for(unsigned iId(0); iId != 6; ++iId){
        _event.getByLabel(edm::InputTag(idTags[iId]), vmH);
        idMaps[iId] = vmH.product();
      }
    }

    susy::MuonCollection& susyCollection(susyEvent_->muons.find(TString(collectionTag).ReplaceAll(":", "_"))->second);

    int imu = 0;
    for(reco::MuonCollection::const_iterator it = muonH->begin(); it != muonH->end(); ++it, ++imu){

      if(it->pt() < muonThreshold_) continue;

      reco::MuonRef muRef(muonH, imu);

      susy::Muon mu;

      Short_t* trackIndices[] = {
        &mu.trackIndex, &mu.standAloneTrackIndex, &mu.combinedTrackIndex,
        &mu.tpfmsTrackIndex, &mu.pickyTrackIndex, &mu.dytTrackIndex
      };

      reco::Muon::MuonTrackType trackTypes[] = {
        reco::Muon::InnerTrack, reco::Muon::OuterTrack, reco::Muon::CombinedTrack,
        reco::Muon::TPFMS, reco::Muon::Picky, reco::Muon::DYT
      };

      mu.type                               = it->type();
      mu.bestTrackType                      = it->muonBestTrackType();
      // SWGuideMuonId#HighPT_Muon
      mu.highPtBestTrackType                = it->innerTrack().isNonnull() ? muon::tevOptimized(*it, 200., 17., 40., 0.25).second : 0;

      for(unsigned iId(0); iId != 6; ++iId)
        if(idMaps[iId] && (*idMaps[iId])[muRef]) mu.qualityFlags |= (0x1 << iId);

      mu.nMatches                           = it->numberOfMatches();
      mu.stationMask                        = it->stationMask();
      mu.nMatchedStations                   = it->numberOfMatchedStations();
      mu.nChambers                          = it->numberOfChambers();

      if(it->combinedMuon().isNonnull()){
        reco::TrackRef combRef(it->combinedMuon());
        mu.nValidHits                         = combRef->hitPattern().numberOfValidHits();
        mu.nValidTrackerHits                  = combRef->hitPattern().numberOfValidTrackerHits();
        mu.nValidMuonHits                     = combRef->hitPattern().numberOfValidMuonHits();
        mu.nPixelLayersWithMeasurement        = combRef->hitPattern().pixelLayersWithMeasurement();
        mu.nStripLayersWithMeasurement        = combRef->hitPattern().stripLayersWithMeasurement();
      }

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

      mu.caloCompatibility                  = it->caloCompatibility();
      mu.segmentCompatibility               = muon::segmentCompatibility(*it);

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

      for(unsigned iT(0); iT != sizeof(trackTypes) / sizeof(reco::Muon::MuonTrackType); ++iT){
        reco::Muon::MuonTrackType trackType(trackTypes[iT]);

        if(it->isAValidMuonTrack(trackType))
          *trackIndices[iT] = fillTrack(it->muonTrack(trackType));
      }

      if(it->isPFMuon()){
        reco::Candidate::LorentzVector pfP4(it->pfP4());
        mu.momentum.SetXYZT(pfP4.px(), pfP4.py(), pfP4.pz(), pfP4.e());
      }
      else
        mu.momentum.SetXYZT(it->p4().px(), it->p4().py(), it->p4().pz(), it->p4().e());

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

    TString jetAlgo(collectionTag);
    jetAlgo.ReplaceAll("CaloJets", "");

    std::string jetIdTag((jetAlgo + "JetID").Data());
    std::string jetCorrName((jetAlgo + "Calo").Data());
    TString jetAlgoUC(jetAlgo);
    jetAlgoUC.ToUpper();
    std::string jecParamsName((jetAlgoUC + "Calo").Data());

    edm::Handle<edm::ValueMap<reco::JetID> > jetIdH;
    _event.getByLabel(edm::InputTag(jetIdTag), jetIdH);

    JetCorrector const* corrL2L3(0);
    JetCorrector const* corrL1L2L3(0);
    if(_event.isRealData()){
      corrL2L3  = JetCorrector::getJetCorrector(jetCorrName + "L2L3Residual", _eventSetup);
      corrL1L2L3  = JetCorrector::getJetCorrector(jetCorrName + "L1L2L3Residual", _eventSetup);
    }
    else{
      corrL2L3  = JetCorrector::getJetCorrector(jetCorrName + "L2L3", _eventSetup);
      corrL1L2L3  = JetCorrector::getJetCorrector(jetCorrName + "L1L2L3", _eventSetup);
    }

    edm::ESHandle<JetCorrectorParametersCollection> jecParamsHndl;
    _eventSetup.get<JetCorrectionsRecord>().get(jecParamsName, jecParamsHndl);
    JetCorrectionUncertainty jecUncert((*jecParamsHndl)["Uncertainty"]);

    susy::CaloJetCollection& susyCollection(susyEvent_->caloJets.find(jetAlgo)->second);

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

      jecUncert.setJetEta(corrP4.Eta());
      jecUncert.setJetPt(corrP4.Pt());
      try{
        jet.jecUncertainty = jecUncert.getUncertainty(true); // true => error high, false => error low. Only symmetric errors are provided so far
      }
      catch(cms::Exception& e){
        // JetCorrectionUncertainty throws when the arguments are out-of-bounds (JetCorrector silently fails by returning unity)
        if(e.category() == "SimpleJetCorrectionUncertainty"){
          edm::LogWarning("SimpleJetCorrectionUncertainty") << "CaloJet Eta = " << corrP4.Eta() << " Pt = " << corrP4.Pt() << " is out of range for JEC uncertainty determination";
          jet.jecUncertainty = -1.;
        }
        else
          throw;
      }

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

  if(debugLevel_ > 0) edm::LogInfo(name()) << "fillPFJets";

  edm::Handle<edm::View<reco::PFJet> > jetH;
  for(unsigned iJetC=0; iJetC < pfJetCollectionTags_.size(); iJetC++) {
    std::string& collectionTag(pfJetCollectionTags_[iJetC]);

    _event.getByLabel(edm::InputTag(collectionTag), jetH);

    if(debugLevel_ > 1) edm::LogInfo(name()) << "fillPFJets: size of " << collectionTag << " JetCollection : " << jetH->size();

    TString jetAlgo(collectionTag);
    jetAlgo.Remove(jetAlgo.Index("PF")); // e.g. "ak5"

    std::string jetCorrName((jetAlgo + "PF").Data());
    TString jetAlgoUC(jetAlgo);
    jetAlgoUC.ToUpper();
    std::string jecParamsName((jetAlgoUC + "PF").Data());
    TString susyColName(jetAlgo);
    if(collectionTag.find("chs") != std::string::npos){
      jetCorrName += "chs";
      jecParamsName += "chs";
      susyColName += "chs";
    }

    JetCorrector const* corrL2L3(0);
    JetCorrector const* corrL1FastL2L3(0);
    if(_event.isRealData()){
      corrL2L3  = JetCorrector::getJetCorrector(jetCorrName + "L2L3Residual", _eventSetup);
      corrL1FastL2L3  = JetCorrector::getJetCorrector(jetCorrName + "L1FastL2L3Residual", _eventSetup);
    }
    else{
      corrL2L3  = JetCorrector::getJetCorrector(jetCorrName + "L2L3", _eventSetup);
      corrL1FastL2L3  = JetCorrector::getJetCorrector(jetCorrName + "L1FastL2L3", _eventSetup);
    }

    edm::ESHandle<JetCorrectorParametersCollection> jecParamsHndl;
    _eventSetup.get<JetCorrectionsRecord>().get(jecParamsName, jecParamsHndl);
    JetCorrectionUncertainty jecUncert((*jecParamsHndl)["Uncertainty"]);

    // Get b-tag information
    VString& bTagTags(bTagCollectionTags_[collectionTag]);
    reco::JetFloatAssociation::Container const* bTagDiscriminators[susy::nBTagDiscriminators];
    for(unsigned i(0); i != susy::nBTagDiscriminators; ++i) bTagDiscriminators[i] = 0;
    if(bTagTags.size() != 0){
      for(unsigned iT(0); iT != susy::nBTagDiscriminators; ++iT){
        edm::Handle<reco::JetFloatAssociation::Container> hndl;
        _event.getByLabel(edm::InputTag(bTagTags[iT]), hndl);
        bTagDiscriminators[iT] = hndl.product();
      }
    }

    // Get quark-gluon discrimination information
    VString& qgTagTags(qgTagCollectionTags_[collectionTag]);
    edm::ValueMap<float> const* qgDiscriminators[susy::nQGDiscriminators];
    for(unsigned i(0); i != susy::nQGDiscriminators; ++i) qgDiscriminators[i] = 0;
    if(qgTagTags.size() != 0){
      for(unsigned iT(0); iT != susy::nQGDiscriminators; ++iT){
        edm::Handle<edm::ValueMap<float> > hndl;
        _event.getByLabel(edm::InputTag(qgTagTags[iT]), hndl);
        if(hndl.isValid()) qgDiscriminators[iT] = hndl.product();
      }
    }

    // Flavour matching for MC
    reco::JetFlavourMatchingCollection const* flavMatchAlg(0);
    reco::JetFlavourMatchingCollection const* flavMatchPhy(0);
    if(!_event.isRealData() && jetFlavourMatchingTags_.find(collectionTag) != jetFlavourMatchingTags_.end()){
      edm::Handle<reco::JetFlavourMatchingCollection> flavMatchH;
      std::pair<std::string, std::string>& tagPair(jetFlavourMatchingTags_[collectionTag]);

      _event.getByLabel(edm::InputTag(tagPair.first), flavMatchH);
      flavMatchAlg = flavMatchH.product();
      _event.getByLabel(edm::InputTag(tagPair.second), flavMatchH);
      flavMatchPhy = flavMatchH.product();
    }

    // Get Pileup Jet ID information
    VString& puJetIdAlgoTags(puJetIdCollectionTags_[collectionTag]);
    edm::ValueMap<float> const* puJetIdMVAValues[susy::nPUJetIdAlgorithms];
    edm::ValueMap<int> const* puJetIdFlags[susy::nPUJetIdAlgorithms];
    if(puJetIdAlgoTags.size() != 0){
      for(unsigned iA(0); iA != susy::nPUJetIdAlgorithms; ++iA){
        if(puJetIdAlgoTags[iA] == "") continue;
        edm::Handle<edm::ValueMap<float> > discHndl;
        edm::Handle<edm::ValueMap<int> > idHndl;
        _event.getByLabel(edm::InputTag(puJetIdAlgoTags[iA] + "Discriminant"), discHndl);
        _event.getByLabel(edm::InputTag(puJetIdAlgoTags[iA] + "Id"), idHndl);
        puJetIdMVAValues[iA] = discHndl.product();
        puJetIdFlags[iA] = idHndl.product();
      }
    }

    susy::PFJetCollection& susyCollection(susyEvent_->pfJets.find(susyColName)->second);

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

      jecUncert.setJetEta(corrP4.Eta());
      jecUncert.setJetPt(corrP4.Pt());
      try{
        jet.jecUncertainty = jecUncert.getUncertainty(true); // true => error high, false => error low. Only symmetric errors are provided so far
      }
      catch(cms::Exception& e){
        // JetCorrectionUncertainty throws when the arguments are out-of-bounds (JetCorrector silently fails by returning unity)
        if(e.category() == "SimpleJetCorrectionUncertainty"){
          edm::LogWarning("SimpleJetCorrectionUncertainty") << "PFJet Eta = " << corrP4.Eta() << " Pt = " << corrP4.Pt() << " is out of range for JEC uncertainty determination";
          jet.jecUncertainty = -1.;
        }
        else
          throw;
      }

      // add btag for this jet
      if(bTagTags.size() != 0){
        for(unsigned iT(0); iT != susy::nBTagDiscriminators; ++iT){
          try{
            jet.bTagDiscriminators[iT] = (*(bTagDiscriminators[iT]))[jetRef];
          }
          catch(cms::Exception& e){
            if(e.category() == "InvalidReference")
              edm::LogWarning("InvalidReference") << "Btag discriminator " << bTagTags[iT] << " does not exist for collection " << collectionTag << " instance " << ijet;
            else
              throw;
          }
        }
      }

      if(qgTagTags.size() != 0){
        for(unsigned iT(0); iT != susy::nQGDiscriminators; ++iT){
          if(!qgDiscriminators[iT]) continue;
          try{
            jet.qgDiscriminators[iT] = (*(qgDiscriminators[iT]))[jetRef];
          }
          catch(cms::Exception& e){
            if(e.category() == "InvalidReference")
              edm::LogWarning("InvalidReference") << "Quark-gluon discriminator " << qgTagTags[iT] << " does not exist for collection " << collectionTag << " instance " << ijet;
            else
              throw;
          }
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
      if(!susyEvent_->isRealData && flavMatchAlg && flavMatchPhy){
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
      if(puJetIdAlgoTags.size() != 0){
        for(unsigned iA(0); iA != susy::nPUJetIdAlgorithms; ++iA){
          if(puJetIdAlgoTags[iA] == "") continue;
          try{
            jet.puJetIdDiscriminants[iA] = (*puJetIdMVAValues[iA])[jetRef];
            jet.puJetIdFlags[iA] = (*puJetIdFlags[iA])[jetRef];
          }
          catch(cms::Exception& e){
            if(e.category() == "InvalidReference")
              edm::LogWarning("InvalidReference") << "PU jet Id does not exist for collection " << collectionTag << " instance " << ijet;
            else
              throw;
          }
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

    TString jetAlgo(collectionTag);

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

    edm::ESHandle<JetCorrectorParametersCollection> jecParamsHndl;
    _eventSetup.get<JetCorrectionsRecord>().get("AK5JPT", jecParamsHndl);
    JetCorrectionUncertainty jecUncert((*jecParamsHndl)["Uncertainty"]);

    susy::JPTJetCollection& susyCollection(susyEvent_->jptJets.find(jetAlgo)->second);

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

      jecUncert.setJetEta(corrP4.Eta());
      jecUncert.setJetPt(corrP4.Pt());
      try{
        jet.jecUncertainty = jecUncert.getUncertainty(true); // true => error high, false => error low. Only symmetric errors are provided so far
      }
      catch(cms::Exception& e){
        // JetCorrectionUncertainty throws when the arguments are out-of-bounds (JetCorrector silently fails by returning unity)
        if(e.category() == "SimpleJetCorrectionUncertainty"){
          edm::LogWarning("SimpleJetCorrectionUncertainty") << "JPTJet Eta = " << corrP4.Eta() << " Pt = " << corrP4.Pt() << " is out of range for JEC uncertainty determination";
          jet.jecUncertainty = -1.;
        }
        else
          throw;
      }

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
    track.ptError = _trkRef->ptError();
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
    track.ptError = _trkRef->ptModeError();
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

  if(susyEvent_->pfParticles.size() != productStore_.pfCandidates.size())
    throw cms::Exception("RuntimeError") << "Number of buffered PFCandidates does not match the number of pfParticles in the susyEvent";

  if(_partPtr.isNull()) return -1;

  std::pair<PFCandidateStore::iterator, bool> insertion(productStore_.pfCandidates.insert(std::pair<reco::PFCandidatePtr, unsigned>(_partPtr, productStore_.pfCandidates.size())));

  if(!insertion.second) return insertion.first->second;
  else{
    susy::PFParticle pf;

    pf.pdgId         = _partPtr->translateTypeToPdgId(_partPtr->particleId());
    pf.charge        = _partPtr->charge();
    pf.ecalEnergy    = _partPtr->ecalEnergy();
    pf.hcalEnergy    = _partPtr->hcalEnergy();

    pf.vertex.SetXYZ(_partPtr->vx(),_partPtr->vy(),_partPtr->vz());
    pf.momentum.SetXYZT(_partPtr->px(),_partPtr->py(),_partPtr->pz(),_partPtr->energy());

    susyEvent_->pfParticles.push_back(pf);

    return susyEvent_->pfParticles.size() - 1;
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
  delete hcalLaserEventListFilter_;
  hcalLaserEventListFilter_ = 0;

  if(storeTriggerEvents_)
    triggerEvent_->write();

  delete triggerEvent_;
  triggerEvent_ = 0;

  delete susyEvent_;
  susyEvent_ = 0;

  if(susyTree_){
    TFile* outF(susyTree_->GetCurrentFile());

    outF->cd();
    susyTree_->Write();
    delete outF;
    susyTree_ = 0;
  }

}

//define this as a plug-in
DEFINE_FWK_MODULE(SusyNtuplizer);
