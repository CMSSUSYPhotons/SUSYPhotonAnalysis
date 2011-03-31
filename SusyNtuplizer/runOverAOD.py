import FWCore.ParameterSet.Config as cms

process = cms.Process("SNT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_4_1_3/RelValTTbar/GEN-SIM-RECO/START311_V2-v1/0037/648B6AA5-C751-E011-8208-001A928116C6.root')
#                            fileNames = cms.untracked.vstring(
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/442/AC2B6850-5B50-E011-A75A-0030487C5CFA.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/443/E4041B08-5850-E011-BF52-0030487CD704.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/444/92F211A3-AC50-E011-9850-0030487CD906.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/445/3465DB0B-5750-E011-912F-003048F118DE.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/446/C80AC79A-5C50-E011-B2B9-003048F1BF68.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/447/EE2523F0-5550-E011-9FDE-001617C3B76E.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/449/CEC7419A-4C50-E011-B232-0030486733B4.root')

)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_311_V2::All'
#process.GlobalTag.globaltag = 'GR_P_V16::All'

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.ak5CaloL1Offset.useCondDB = False
process.ak5PFL1Offset.useCondDB = False
process.ak5JPTL1Offset.useCondDB = False

process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)
#process.susyNtuplizer.recoMode = cms.bool(True)
process.susyNtuplizer.outputFileName = cms.string("susyEvents.root")
process.susyNtuplizer.storeGenInfos = cms.bool(True)

process.jec = cms.Sequence(process.ak5CaloJetsL2L3 * process.ak7CaloJetsL2L3 * process.kt4CaloJetsL2L3 * process.kt6CaloJetsL2L3 *
                           process.ak5CaloJetsL2L3Residual * process.ak7CaloJetsL2L3Residual * process.kt4CaloJetsL2L3Residual * process.kt6CaloJetsL2L3Residual *
                           process.ak5PFJetsL2L3 * process.ak7PFJetsL2L3 * process.kt4PFJetsL2L3 * process.kt6PFJetsL2L3 *
                           process.ak5PFJetsL2L3Residual * process.ak7PFJetsL2L3Residual * process.kt4PFJetsL2L3Residual * process.kt6PFJetsL2L3Residual *
                           process.ak5JPTJetsL2L3 * process.ak5JPTJetsL2L3Residual)

process.p = cms.Path( process.jec * process.susyNtuplizer )

