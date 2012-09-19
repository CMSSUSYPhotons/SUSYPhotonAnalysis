import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

# change this to 0 if you run on FullSim MC
isFastSim = 1

# These are fixes for the JetProbability b-tagger calibrations as recommended by BTV.
# See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration#Calibration_in_52x_Data_and_MC
# Leaving these as all 0s results in no changes to the GlobalTag used.
is52xData = 0
is53xData = 0
is53xMC = 0

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.suppressWarning = cms.untracked.vstring('newSecondaryVertexTagInfos')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring( )
                            )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections/Configuration/JetCorrectionServices_cff")
process.load("JetMETCorrections/Type1MET/caloMETCorrections_cff")

# SusyNtuplizer options
process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)
if isFastSim:
    process.susyNtuplizer.isFastSim = cms.bool(True)
# For FNAL users
#process.susyNtuplizer.photonSCRegressionWeights = "/eos/uscms/store/user/lpcpjm/NtuplizerData/gbrv3ph_52x.root"

process.metAnalysisSequence = cms.Sequence(process.producePFMETCorrections*
                                           process.produceCaloMETCorrections)

#Calculate rho restricted to barrel for photon pileup subtraction
process.kt6PFJetsRhoBarrelOnly = process.kt4PFJets.clone(
    src = cms.InputTag('particleFlow'),
    rParam = cms.double(0.6),
    #Eta range of jets to be considered for Rho calculation
    #Should be at most (jet acceptance - jet radius)
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax=cms.double(1.4442),
    #Eta range of ghost jets to be considered for Rho calculation - must be greater than Rho_EtaMax
    Ghost_EtaMax=cms.double(2.5)
    )

#Calculate rho restricted to barrel for photon pileup subtraction
process.kt6PFJetsRho25 = process.kt4PFJets.clone(
    src = cms.InputTag('particleFlow'),
    rParam = cms.double(0.6),
    #Eta range of jets to be considered for Rho calculation
    #Should be at most (jet acceptance - jet radius)
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax=cms.double(2.5)
    )

# HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

# HCAL laser events filter
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
process.hcalLaserEventFilter.taggingMode = cms.bool(True)

# EcalDeadCellTriggerPrimitiveFilter
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

# EcalDeadCellBoundaryEnergyFilter
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)

# Tracking failure filter
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3Residual')
process.trackingFailureFilter.taggingMode = cms.bool(True)

# EE Bad SC Filter
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.eeBadScFilter.taggingMode = cms.bool(True)

# EE ring of fire
process.load('RecoMET.METFilters.eeNoiseFilter_cfi')
process.eeNoiseFilter.taggingMode = cms.bool(True)

# Inconsistent muon pf candidate filter
process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')
process.inconsistentMuonPFCandidateFilter.taggingMode = cms.bool(True)

# Greedy muon pf candidate filter
process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')
process.greedyMuonPFCandidateFilter.taggingMode = cms.bool(True)

#Add up all MET filters
if realData or not isFastSim:
    process.metFiltersSequence = cms.Sequence(
    	process.HBHENoiseFilterResultProducer *
    	process.hcalLaserEventFilter *
    	process.EcalDeadCellTriggerPrimitiveFilter *
    	process.EcalDeadCellBoundaryEnergyFilter *
    	process.goodVertices *
   	process.trackingFailureFilter *
	process.eeBadScFilter *
	process.eeNoiseFilter *
	process.inconsistentMuonPFCandidateFilter *
	process.greedyMuonPFCandidateFilter
	)
else:
    process.metFiltersSequence = cms.Sequence(
    	process.hcalLaserEventFilter *
    	process.EcalDeadCellTriggerPrimitiveFilter *
    	process.EcalDeadCellBoundaryEnergyFilter *
    	process.goodVertices *
   	process.trackingFailureFilter *
	process.eeBadScFilter *
	process.eeNoiseFilter *
	process.inconsistentMuonPFCandidateFilter *
	process.greedyMuonPFCandidateFilter
	)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.myPartons = cms.EDProducer("PartonSelector",   
                                   withLeptons = cms.bool(False),
                                   src = cms.InputTag("genParticles")
                                   )

process.flavourByRef = cms.EDProducer("JetPartonMatcher",
                                      jets = cms.InputTag("ak5PFJets"),
                                      coneSizeToAssociate = cms.double(0.3),
                                      partons = cms.InputTag("myPartons")
                                      )

process.flavourAssociationAlg = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("flavourByRef"),
    physicsDefinition = cms.bool(False)
    )

process.flavourAssociationPhy = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("flavourByRef"),
    physicsDefinition = cms.bool(True)
    )

process.JetFlavourMatching = cms.Sequence(
    process.myPartons *
    process.flavourByRef *
    process.flavourAssociationAlg *
    process.flavourAssociationPhy
)

# PileupJetId
from CMGTools.External.pujetidsequence_cff import puJetId, puJetMva

process.recoPuJetId = puJetId.clone(
   jets = cms.InputTag("ak5PFJets"),
   applyJec = cms.bool(True),
   inputIsCorrected = cms.bool(False)
)

process.recoPuJetMva = puJetMva.clone(
   jets = cms.InputTag("ak5PFJets"),
   applyJec = cms.bool(True),
   inputIsCorrected = cms.bool(False),                
   jetids = cms.InputTag("recoPuJetId")
)

process.recoPuJetIdSqeuence = cms.Sequence(
    process.recoPuJetId * 
    process.recoPuJetMva
)

if realData:
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2012C/DoublePhoton/AOD/PromptReco-v2/000/202/016/D6785FDC-23F5-E111-9DA4-0030486780B4.root'
        )
    process.GlobalTag.globaltag = 'GR_R_53_V14::All'

    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3Residual")
    # JEC for data
    process.jet = cms.Sequence(
	process.ak5PFJets *
        # CaloJets
        process.ak5CaloJetsL2L3Residual * process.ak5CaloJetsL1L2L3Residual *
        # PFJets
        process.ak5PFJetsL2L3Residual * process.ak5PFJetsL1FastL2L3Residual *
        # Barrel only Rho calculation
        process.kt6PFJetsRhoBarrelOnly *
	# Rho25 calculation
	process.kt6PFJetsRho25
        )

else:
    process.source.fileNames = cms.untracked.vstring(
	'dcap:///pnfs/cms/WAX/11/store/mc/Summer12/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v1/0001/FE8F81B3-C494-E111-B50D-003048D476BC.root'
        )
    process.GlobalTag.globaltag = 'START53_V11::All'
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3")
    # JEC for MC
    process.jet = cms.Sequence(
	process.ak5PFJets *
        # CaloJets
        process.ak5CaloJetsL2L3 * process.ak5CaloJetsL1L2L3 *
        # PFJets
        process.ak5PFJetsL2L3 * process.ak5PFJetsL1FastL2L3 *
        # Barrel only Rho calculation
        process.kt6PFJetsRhoBarrelOnly *
	# Rho25 calculation
	process.kt6PFJetsRho25 *
	# Jet-flavour matching for b-tagging
	process.JetFlavourMatching
        )
    process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3')
    if isFastSim:
	process.susyNtuplizer.muonIDCollectionTags = cms.vstring()
	process.susyNtuplizer.muonCollectionTags = cms.vstring("muons")

# IsoDeposit
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.isoDeposit = cms.Sequence( process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)


# b-tagging stuff ===============================================================
# Re-run b-tagging with PFJets as input

# b-tagging general configuration
process.load("RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")

# create a new jets and tracks associaiton
process.newJetTracksAssociatorAtVertex = process.ic5PFJetTracksAssociatorAtVertex.clone()
process.newJetTracksAssociatorAtVertex.jets = cms.InputTag("ak5PFJets")
process.newJetTracksAssociatorAtVertex.tracks = cms.InputTag("generalTracks")

# impact parameter b-tag
process.newImpactParameterTagInfos = process.impactParameterTagInfos.clone()
process.newImpactParameterTagInfos.jetTracks = "newJetTracksAssociatorAtVertex"
process.newTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone()
process.newTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
process.newTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone()
process.newTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
process.newJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone()
process.newJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
process.newJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone()
process.newJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )

# secondary vertex b-tag
process.newSecondaryVertexTagInfos = process.secondaryVertexTagInfos.clone()
process.newSecondaryVertexTagInfos.trackIPTagInfos = "newImpactParameterTagInfos"
process.newSimpleSecondaryVertexBJetTags = process.simpleSecondaryVertexBJetTags.clone()
process.newSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
process.newCombinedSecondaryVertexBJetTags = process.combinedSecondaryVertexBJetTags.clone()
process.newCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), 
cms.InputTag("newSecondaryVertexTagInfos") )
process.newCombinedSecondaryVertexMVABJetTags = process.combinedSecondaryVertexMVABJetTags.clone()
process.newCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), 
cms.InputTag("newSecondaryVertexTagInfos") )

# soft electron b-tag
process.newSoftElectronTagInfos = process.softElectronTagInfos.clone()
process.newSoftElectronTagInfos.jets = "ak5PFJets"
process.newSoftElectronBJetTags = process.softElectronBJetTags.clone()
process.newSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftElectronTagInfos") )

# soft muon b-tag
process.newSoftMuonTagInfos = process.softMuonTagInfos.clone()
process.newSoftMuonTagInfos.jets = "ak5PFJets"
process.newSoftMuonBJetTags = process.softMuonBJetTags.clone()
process.newSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
process.newSoftMuonByIP3dBJetTags = process.softMuonByIP3dBJetTags.clone()
process.newSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
process.newSoftMuonByPtBJetTags = process.softMuonByPtBJetTags.clone()
process.newSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )

process.newJetTracksAssociator = cms.Sequence(
    process.newJetTracksAssociatorAtVertex
)

process.newJetBtaggingIP = cms.Sequence(
    process.newImpactParameterTagInfos * (
        process.newTrackCountingHighEffBJetTags +
        process.newTrackCountingHighPurBJetTags +
        process.newJetProbabilityBJetTags +
        process.newJetBProbabilityBJetTags
    )
)

process.newJetBtaggingSV = cms.Sequence(
    process.newImpactParameterTagInfos *
    process.newSecondaryVertexTagInfos * (
        process.newSimpleSecondaryVertexBJetTags +
        process.newCombinedSecondaryVertexBJetTags +
        process.newCombinedSecondaryVertexMVABJetTags
    )
)

process.newJetBtaggingEle = cms.Sequence(
    process.softElectronCands * 
    process.newSoftElectronTagInfos *
    process.newSoftElectronBJetTags
)

process.newJetBtaggingMu = cms.Sequence(
    process.newSoftMuonTagInfos * (
        process.newSoftMuonBJetTags +
        process.newSoftMuonByIP3dBJetTags +
        process.newSoftMuonByPtBJetTags
    )
)

process.newJetBtagging = cms.Sequence(
    process.newJetBtaggingIP +
    process.newJetBtaggingSV +
    process.newJetBtaggingEle +
    process.newJetBtaggingMu
)

if is52xData and realData:
    process.GlobalTag.toGet = cms.VPSet(
        cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
        tag = cms.string("TrackProbabilityCalibration_2D_2012DataTOT_v1_offline"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
        tag = cms.string("TrackProbabilityCalibration_3D_2012DataTOT_v1_offline"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )

if is53xData and realData:
    process.GlobalTag.toGet = cms.VPSet(
  	cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
       	tag = cms.string("TrackProbabilityCalibration_2D_Data53X_v2"),
       	connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       	tag = cms.string("TrackProbabilityCalibration_3D_Data53X_v2"),
       	connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )

if is53xMC and not realData:
    process.GlobalTag.toGet = cms.VPSet(
  	cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
        tag = cms.string("TrackProbabilityCalibration_2D_MC53X_v2"),
        connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")),
    cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
       	tag = cms.string("TrackProbabilityCalibration_3D_MC53X_v2"),
       	connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU"))
    )

process.p = cms.Path(
    process.jet *
    process.newJetTracksAssociator *
    process.newJetBtagging *
    process.metAnalysisSequence *
    process.metFiltersSequence *
    process.recoPuJetIdSqeuence *
    process.isoDeposit *
    process.susyNtuplizer
    )
