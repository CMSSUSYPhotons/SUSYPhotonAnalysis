import FWCore.ParameterSet.Config as cms

susyTriggerNtuplizer = cms.EDAnalyzer("SusyTriggerNtuplizer",
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD"),
    filterTags = cms.untracked.vstring(),
    fileName = cms.untracked.string("susyTriggers.root")
)
