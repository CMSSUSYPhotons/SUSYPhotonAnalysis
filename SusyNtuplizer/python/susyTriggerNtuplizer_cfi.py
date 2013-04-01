import FWCore.ParameterSet.Config as cms

susyTriggerNtuplizer = cms.EDAnalyzer("SusyTriggerNtuplizer",
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD"),
    fileName = cms.untracked.string("susyTriggers.root")
)
