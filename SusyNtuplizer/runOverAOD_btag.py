import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring( )
                            )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
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

# HBHENoiseFilterResultProducer
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')

# HCAL laser events filter
process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")

# EcalDeadCellTriggerPrimitiveFilter
process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')

# EcalDeadCellBoundaryEnergyFilter
process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')

# Tracking failure filter
# this is not recommended at the moment, but let's keep it for later use
process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3Residual')

#Add up all MET filters
process.metFiltersSequence = cms.Sequence(
    process.HBHENoiseFilterResultProducer *
    process.hcalLaserEventFilter *
    process.EcalDeadCellTriggerPrimitiveFilter *
    process.EcalDeadCellBoundaryEnergyFilter *
    process.goodVertices *
    process.trackingFailureFilter
)

if realData:
    process.source.fileNames = cms.untracked.vstring(
	'/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/706/DA8B61A9-BE83-E111-8BCB-001D09F2906A.root'
        )
    process.GlobalTag.globaltag = 'GR_R_52_V7::All'

    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3Residual")
    # JEC for data
    process.jet = cms.Sequence(
        # CaloJets
        process.ak5CaloJetsL2L3Residual * process.ak5CaloJetsL1L2L3Residual *
        # PFJets
        process.ak5PFJetsL2L3Residual * process.ak5PFJetsL1FastL2L3Residual *
        # Barrel only Rho calculation
        process.kt6PFJetsRhoBarrelOnly
        )

else:
    process.source.fileNames = cms.untracked.vstring(
	'dcap:///pnfs/cms/WAX/11/store/mc/Summer12/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S7_START52_V9-v1/0001/FE8F81B3-C494-E111-B50D-003048D476BC.root'
        )
    process.GlobalTag.globaltag = 'START52_V9::All'
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3")
    # JEC for MC
    process.jet = cms.Sequence(
        # CaloJets
        process.ak5CaloJetsL2L3 * process.ak5CaloJetsL1L2L3 *
        # PFJets
        process.ak5PFJetsL2L3 * process.ak5PFJetsL1FastL2L3 *
        # Barrel only Rho calculation
        process.kt6PFJetsRhoBarrelOnly
        )
    process.trackingFailureFilter.JetSource = cms.InputTag('ak5PFJetsL2L3')

# IsoDeposit
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

process.isoDeposit = cms.Sequence( process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)


# b-tagging stuff ===============================================================
# Re-run b-tagging with PFJets as input

# b-tagging general configuration
process.load("RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")

# create a new jets and tracks associaiton
process.newJetTracksAssociatorAtVertex = process.ic5JetTracksAssociatorAtVertex.clone()
process.newJetTracksAssociatorAtVertex.jets = "ak5PFJets"
#process.newJetTracksAssociatorAtVertex.tracks = "generalTracks"

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

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.myPartons = cms.EDProducer("PartonSelector",   
                                   withLeptons = cms.bool(False),
                                   src = cms.InputTag("genParticles")
                                   )

process.flavourByRef = cms.EDProducer("JetPartonMatcher",
                                      #jets = cms.InputTag("iterativeCone5CaloJets"),
                                      #jets = cms.InputTag("ak5CaloJetsL2L3"),
                                      jets = cms.InputTag("ak5PFJetsL2L3"),
                                      coneSizeToAssociate = cms.double(0.3),
                                      partons = cms.InputTag("myPartons")
                                      )

process.JetFlavourMatching = cms.Sequence(
    process.myPartons *
    process.flavourByRef
)

if realData:
    process.p = cms.Path(
	process.newJetTracksAssociator *
	process.newJetBtagging *
	process.metAnalysisSequence *
	process.jet *
	process.metFiltersSequence *
	process.isoDeposit *
	process.susyNtuplizer
	)

else:
    process.p = cms.Path(
	process.newJetTracksAssociator *
	process.newJetBtagging *
	process.metAnalysisSequence *
	process.jet *
	process.JetFlavourMatching *
	process.metFiltersSequence *
	process.isoDeposit *
	process.susyNtuplizer
	)

