import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

#Flag HCAL Noise events as such
process.HBHENoiseFilterResultProducer = cms.EDProducer(
    'HBHENoiseFilterResultProducer',
    label = cms.InputTag('hcalnoise','','RECO'),
    minRatio = cms.double(-999),
    maxRatio = cms.double(999),
    minHPDHits = cms.int32(17),
    minRBXHits = cms.int32(999),
    minHPDNoOtherHits = cms.int32(10),
    minZeros = cms.int32(10),
    minHighEHitTime = cms.double(-9999.0),
    maxHighEHitTime = cms.double(9999.0),
    maxRBXEMF = cms.double(-999.0),
    minNumIsolatedNoiseChannels = cms.int32(99999),
    minIsolatedNoiseSumE = cms.double(99999.),
    minIsolatedNoiseSumEt = cms.double(99999.),
    useTS4TS5 = cms.bool(True)
    )

#Flag ECAL dead cell events as such
from JetMETAnalysis.ecalDeadCellTools.EcalDeadCellEventFilter_cfi import *
process.ecalDeadCellTPfilter = EcalDeadCellEventFilter.clone(
    tpDigiCollection = cms.InputTag("ecalTPSkim"),
    # when activated, the filter does not filter event.
    # the filter is however storing a bool in the event, that can be used to take the
    # filtering decision a posteriori
    taggingMode = cms.bool( True ),
    debug = cms.untracked.bool( False ),
    etValToBeFlagged = cms.double(63.75),
    doEEfilter = cms.untracked.bool( False ), # turn it on by default
    makeProfileRoot = cms.untracked.bool( False )
    )
#Add up all MET filters
process.metFiltersSequence = cms.Sequence( process.HBHENoiseFilterResultProducer * process.ecalDeadCellTPfilter)


if realData:
    process.source.fileNames = cms.untracked.vstring(
        '/store/data/Run2011A/DoubleElectron/AOD/May10ReReco-v1/0000/003D325C-547B-E011-81D4-001A928116C2.root'
        )
    process.GlobalTag.globaltag = 'GR_R_42_V21A::All'

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
        #'dcap:///pnfs/cms/WAX/resilient/lpcpjm/PrivateMC/BinoSignalPoints_5_7_11/reco/1250_1200_225/reco_1250_1200_225_1.root'
        #'dcap:////pnfs/cms/WAX/11/store/user/lpcpjm/PrivateMC/FastSim/binolikegrid2/gen/fastsim_400_400_375.root'
        'file:/eos/uscms/store/user/lpcpjm/PrivateMC/FastSim/428p7/Spectra_gB/gen/fastsim_2000_1015_55.root'

        )
    process.GlobalTag.globaltag = 'START42_V15B::All'
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


    


# b-tagging stuff =====================================================================
# It seems that the default b-tagging info coresponds to CaloJets, while we want PFJets
# Re-run b-tagging with PFJets as input


# b-tagging general configuration
process.load("RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi")
process.load("RecoBTag.Configuration.RecoBTag_cff")

# create a new jets and tracks association
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
process.newCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )
process.newCombinedSecondaryVertexMVABJetTags = process.combinedSecondaryVertexMVABJetTags.clone()
process.newCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )

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


    

if realData:
    print "Running on data ......................................................."
    process.p = cms.Path(
        process.newJetTracksAssociator *
        process.newJetBtagging *
        process.metFiltersSequence *
        process.metAnalysisSequence *
        process.jet *
        process.susyNtuplizer )
else:
    print "Running on MC ........................................................."
    process.p = cms.Path(
        process.newJetTracksAssociator *
        process.newJetBtagging *
        #process.metFiltersSequence * # doesn not work on MC; comment out according to David Morse
        process.metAnalysisSequence *
        process.jet *
        process.susyNtuplizer )
