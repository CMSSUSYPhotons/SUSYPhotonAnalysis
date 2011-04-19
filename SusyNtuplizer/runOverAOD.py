import FWCore.ParameterSet.Config as cms

process = cms.Process("SNT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(
# MC
    '/store/relval/CMSSW_4_2_0/RelValZTT/GEN-SIM-RECO/START42_V9-v1/0054/107DB9B4-7D5E-E011-91E9-001A92810AEA.root'
# Data
#    '/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/449/CEC7419A-4C50-E011-B232-0030486733B4.root'
                                                              )
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V10::All'
#process.GlobalTag.globaltag = 'GR_R_311_V2::All'

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.ak5CaloL1Offset.useCondDB = False
process.ak5PFL1Offset.useCondDB = False
process.ak5JPTL1Offset.useCondDB = False

process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)
process.susyNtuplizer.recoMode = cms.bool(True)
process.susyNtuplizer.outputFileName = cms.string("susyEvents.root")
process.susyNtuplizer.storeGenInfos = cms.bool(True)
process.susyNtuplizer.electronCollectionTags = cms.vstring("gsfElectrons","pfElectronTranslator:pf")
process.susyNtuplizer.photonCollectionTags = cms.vstring("photons","pfPhotonTranslator:pfphot")

process.jec = cms.Sequence(process.ak5CaloJetsL2L3 * process.ak7CaloJetsL2L3 * process.kt4CaloJetsL2L3 * process.kt6CaloJetsL2L3
                           * process.ak5CaloJetsL2L3Residual * process.ak7CaloJetsL2L3Residual * process.kt4CaloJetsL2L3Residual * process.kt6CaloJetsL2L3Residual
                           * process.ak5PFJetsL2L3 * process.ak7PFJetsL2L3 * process.kt4PFJetsL2L3 * process.kt6PFJetsL2L3
                           * process.ak5PFJetsL2L3Residual * process.ak7PFJetsL2L3Residual * process.kt4PFJetsL2L3Residual * process.kt6PFJetsL2L3Residual
                           * process.ak5JPTJetsL2L3 * process.ak5JPTJetsL2L3Residual)

process.p = cms.Path( process.jec * process.susyNtuplizer )
