import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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
process.susyNtuplizer.storeBtagging = cms.bool(False)
    
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
    noiselabel = cms.InputTag('hcalnoise','','RECO'),
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

# CSCBeamHaloFilter
process.muonsFromCosmics = process.muontiming.clone(
    MuonCollection = cms.InputTag("muonsFromCosmics")
    )
process.MessageLogger.suppressWarning = cms.untracked.vstring("CSCHaloData")
process.CSCBeamHaloFilterResultProducer = cms.Sequence( process.muonsFromCosmics * process.BeamHaloId)

# HCAL laser events filter
from RecoMET.METFilters.hcalLaserEventFilter_cfi import *

# Tracking failure filter
process.goodVerticesForFilter = cms.EDFilter(
    "VertexSelector",
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)
from RecoMET.METFilters.trackingFailureFilter_cfi import *
process.trackingFailureFilterProducer = trackingFailureFilter.clone(
    JetSource = cms.InputTag('ak5PFJetsL2L3Residual'),
    VertexSource = cms.InputTag('goodVerticesForFilter'),
    taggingMode = cms.bool(True)
)
process.trackingFailureFilterSequence = cms.Sequence(
    process.goodVerticesForFilter *
    process.trackingFailureFilterProducer
)

#Add up all MET filters
process.metFiltersSequence = cms.Sequence(
    process.HBHENoiseFilterResultProducer *
    process.ecalDeadCellTPfilter *
#    process.CSCBeamHaloFilterResultProducer *
    process.trackingFailureFilterSequence
)

if realData:
    process.source.fileNames = cms.untracked.vstring(
        #'/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/450/3CC2524E-ED80-E111-BCFB-BCAEC518FF52.root',
        #'/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/995/7493F473-1286-E111-B648-003048F118E0.root',
	#'/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/192/257/4A3E3DCF-BC8E-E111-BA82-5404A640A639.root',
	'/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/190/706/DA8B61A9-BE83-E111-8BCB-001D09F2906A.root'

        )
    process.GlobalTag.globaltag = 'GR_R_52_V7::All'
    #process.GlobalTag.globaltag = 'GR_P_V32::All'

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
        'dcap:///pnfs/cms/WAX/resilient/lpcpjm/PrivateMC/BinoSignalPoints_5_7_11/reco/1250_1200_225/reco_1250_1200_225_1.root'
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


process.p = cms.Path(
    process.metAnalysisSequence *
    process.jet *
    process.metFiltersSequence *
    process.susyNtuplizer
    )
   
outfile = open('config.py','w')
print >> outfile,process.dumpPython()
outfile.close()
