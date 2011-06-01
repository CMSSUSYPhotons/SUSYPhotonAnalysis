import FWCore.ParameterSet.Config as cms

realData = 0

process = cms.Process("RA3")

process.load('FWCore/MessageService/MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if realData:
    process.GlobalTag.globaltag = 'GR_R_42_V14::All'
else:
    process.GlobalTag.globaltag = 'START42_V12::All'


process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

process.readAK5Calo = cms.EDAnalyzer(
    'JetCorrectorDBReader', 
    payloadName    = cms.untracked.string('AK5Calo'),
    printScreen    = cms.untracked.bool(True),
    createTextFile = cms.untracked.bool(True),
    globalTag      = cms.untracked.string('Jec11_V1')
    )

process.readAK5PF = cms.EDAnalyzer(
    'JetCorrectorDBReader', 
    payloadName    = cms.untracked.string('AK5PF'),
    printScreen    = cms.untracked.bool(True),
    createTextFile = cms.untracked.bool(True),
    globalTag      = cms.untracked.string('Jec11_V1')
    )

process.readAK5JPT = cms.EDAnalyzer(
    'JetCorrectorDBReader', 
    payloadName    = cms.untracked.string('AK5JPT'),
    printScreen    = cms.untracked.bool(True),
    createTextFile = cms.untracked.bool(True),
    globalTag      = cms.untracked.string('Jec11_V1')
    )


process.p = cms.Path( process.readAK5Calo * process.readAK5PF * process.readAK5JPT )
