import FWCore.ParameterSet.Config as cms

# change this to 0 if you run on MC files
realData = 1

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

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

if realData:
    process.source.fileNames = cms.untracked.vstring('/store/data/Run2011A/DoubleElectron/AOD/May10ReReco-v1/0000/003D325C-547B-E011-81D4-001A928116C2.root')
    process.GlobalTag.globaltag = 'GR_R_42_V14::All'
else:
    process.source.fileNames = cms.untracked.vstring('/store/relval/CMSSW_4_2_3/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START42_V12-v2/0062/AA819E58-077B-E011-8C9C-0018F3D095FC.root')
    process.GlobalTag.globaltag = 'START42_V12::All'

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.ak5PFL1Fastjet.useCondDB = False

# to avoid crashes due to looking for kt6PFJets:rho
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.kt6PFJets = process.kt4PFJets.clone(
    rParam = cms.double(0.6),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    voronoiRfact = cms.double(0.9)
    )

# SusyNtuplizer options
process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)

# JEC
process.jet = cms.Sequence(
    process.kt6PFJets *
    # CaloJets
    process.ak5CaloJetsL2L3 * process.ak5CaloJetsL1FastL2L3 *
    # PFJets
    process.ak5PFJetsL2L3 * process.ak5PFJetsL1FastL2L3 *
    # JPTJets
    process.ak5JPTJetsL2L3 * process.ak5JPTJetsL1FastL2L3
    )

process.p = cms.Path( process.jet * process.susyNtuplizer )
