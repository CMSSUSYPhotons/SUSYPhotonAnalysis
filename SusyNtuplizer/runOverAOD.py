import FWCore.ParameterSet.Config as cms

process = cms.Process("SNT")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring(
# MC
#    '/store/relval/CMSSW_4_2_3/RelValPhotonJets_Pt_10/GEN-SIM-RECO/START42_V12-v2/0062/AA819E58-077B-E011-8C9C-0018F3D095FC.root'
    '/store/data/Run2011A/DoubleElectron/AOD/May10ReReco-v1/0004/003F6C25-7C7B-E011-8C25-0026189438EF.root'
    )
)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# for DATA 2011 with 4_2_X
#process.GlobalTag.globaltag = 'GR_R_42_V12::All'
# for MC with 4_2_X
#process.GlobalTag.globaltag = 'START42_V12::All'

# for 4_1_3
# 4_2_X jet corrections are not implemented for L2L3Residuals.
# So we are still using the old tags for the time being...
process.GlobalTag.globaltag = 'GR_R_311_V2::All'

process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.ak5CaloL1Offset.useCondDB = False
process.ak5PFL1Offset.useCondDB = False
process.ak5JPTL1Offset.useCondDB = False


# SusyNtuplizer options
process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")

process.jec = cms.Sequence(process.ak5CaloJetsL2L3 * process.ak7CaloJetsL2L3 * process.kt4CaloJetsL2L3 * process.kt6CaloJetsL2L3
                           * process.ak5CaloJetsL2L3Residual * process.ak7CaloJetsL2L3Residual * process.kt4CaloJetsL2L3Residual * process.kt6CaloJetsL2L3Residual
                           * process.ak5PFJetsL2L3 * process.ak7PFJetsL2L3 * process.kt4PFJetsL2L3 * process.kt6PFJetsL2L3
                           * process.ak5PFJetsL2L3Residual * process.ak7PFJetsL2L3Residual * process.kt4PFJetsL2L3Residual * process.kt6PFJetsL2L3Residual
                           * process.ak5JPTJetsL2L3 * process.ak5JPTJetsL2L3Residual)

process.p = cms.Path( process.jec * process.susyNtuplizer )
