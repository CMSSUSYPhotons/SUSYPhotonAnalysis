import FWCore.ParameterSet.Config as cms

globalTag = 'FT_53_V21_AN3'

# PF, Calo, or JPT
jetType = 'PF'

process = cms.Process("RA3")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = globalTag + '::All'

process.saveCorrections = cms.EDAnalyzer(
    'JetCorrectorDBReader', 
    payloadName    = cms.untracked.string('AK5' + jetType),
    printScreen    = cms.untracked.bool(True),
    createTextFile = cms.untracked.bool(True),
    globalTag      = cms.untracked.string(globalTag)
)

process.p = cms.Path(
    process.saveCorrections
)
