import FWCore.ParameterSet.Config as cms

process = cms.Process("SNT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring('/store/relval/CMSSW_4_1_3/RelValTTbar/GEN-SIM-RECO/START311_V2-v1/0037/648B6AA5-C751-E011-8208-001A928116C6.root')
#                            fileNames = cms.untracked.vstring(
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/377/C4B65E07-0C4F-E011-80D6-0030487C7E18.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/379/108678E6-0B4F-E011-B68D-0030487CD178.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/383/609C0CE6-0C4F-E011-9D31-0030487CD812.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/384/BE02D7AB-084F-E011-BBDF-0030487C7392.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/386/64EC0AC8-0A4F-E011-847D-00304879BAB2.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/403/8696468F-364F-E011-BBC1-003048D2BED6.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/404/36E04BE4-374F-E011-ABD1-001D09F23D1D.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/405/CAC914AE-D64F-E011-A7E2-0030487CD718.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/406/D0214F28-604F-E011-8F8A-003048F11C5C.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/410/866EF731-604F-E011-A5D1-00304879BAB2.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/413/A2F31523-EB4F-E011-BE17-0030487CBD0A.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/421/C6C9A2FB-774F-E011-9A7A-003048D2BD66.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/423/70D787D7-7C4F-E011-A7D5-0030487C8E00.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/425/B8B05C57-844F-E011-866B-003048F117EA.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/427/CC551DB3-814F-E011-898F-0030487CD718.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/428/E82B74C9-B34F-E011-9E8C-001617C3B5D8.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/431/B4D4627D-2D50-E011-B7B7-0030487A3C92.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/432/C0F4946D-1150-E011-8B81-001D09F2546F.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/433/24F4D05C-B54F-E011-82A3-001617C3B706.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/439/FEE9F0FA-B34F-E011-A364-003048F1C836.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/442/AC2B6850-5B50-E011-A75A-0030487C5CFA.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/443/E4041B08-5850-E011-BF52-0030487CD704.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/444/92F211A3-AC50-E011-9850-0030487CD906.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/445/3465DB0B-5750-E011-912F-003048F118DE.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/446/C80AC79A-5C50-E011-B2B9-003048F1BF68.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/447/EE2523F0-5550-E011-9FDE-001617C3B76E.root',
#'/store/data/Run2011A/PhotonHad/AOD/PromptReco-v1/000/160/449/CEC7419A-4C50-E011-B232-0030486733B4.root')

)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = cms.int32(0)
#process.susyNtuplizer.printTreeVariables = cms.int32(2)
#process.susyNtuplizer.recoPatMode = cms.bool(True)
process.susyNtuplizer.outputFileName = cms.string("susyEvents.root")
process.susyNtuplizer.storeGenInfos = cms.bool(True)

process.p = cms.Path( process.susyNtuplizer )

