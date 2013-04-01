import FWCore.ParameterSet.Config as cms

susyNtuplizer = cms.EDAnalyzer('SusyNtuplizer',
    lumiSummaryTag = cms.string("lumiProducer"),
    vtxCollectionTag = cms.string("offlinePrimaryVertices"),
    trackCollectionTag = cms.string("generalTracks"),
    pfCandidateCollectionTag = cms.string("particleFlow"),
    muonCollectionTags = cms.vstring(
        "muons",
        "muonsFromCosmics"
    ),
    electronCollectionTags = cms.vstring(
        "gsfElectrons"
    ),
    photonCollectionTags = cms.vstring(
        "photons",
        "pfPhotonTranslator:pfphot"
    ),
    caloJetCollectionTags = cms.vstring(
        "ak5CaloJets"
    ),
    pfJetCollectionTags = cms.vstring(
        "ak5PFJets"
    ),
    jptJetCollectionTags = cms.vstring(),
    muonIdTags = cms.PSet(
        muons = cms.vstring(
            "muons:muidAllArbitrated",
            "muons:muidGMStaChiCompatibility",
            "muons:muidGMTkChiCompatibility",
            "muons:muidGMTkKinkTight",
            "muons:muidGlobalMuonPromptTight",
            "muons:muidTM2DCompatibilityLoose",
            "muons:muidTM2DCompatibilityTight",
            "muons:muidTMLastStationAngLoose",
            "muons:muidTMLastStationAngTight",
            "muons:muidTMLastStationLoose",
            "muons:muidTMLastStationOptimizedLowPtLoose",
            "muons:muidTMLastStationOptimizedLowPtTight",
            "muons:muidTMLastStationTight",
            "muons:muidTMOneStationAngLoose",
            "muons:muidTMOneStationAngTight",
            "muons:muidTMOneStationLoose",
            "muons:muidTMOneStationTight",
            "muons:muidTrackerMuonArbitrated"
        )
    ),
    electronIdTags = cms.PSet(
        gsfElectrons = cms.vstring(
            "eidLoose",
            "eidRobustHighEnergy",
            "eidRobustLoose",
            "eidRobustTight",
            "eidTight",
            "mvaTrigV0",
            "mvaNonTrigV0"
        )
    ),
    photonIdTags = cms.PSet(
        photons = cms.vstring(
            "PhotonIDProd:PhotonCutBasedIDLoose",
            "PhotonIDProd:PhotonCutBasedIDLooseEM",
            "PhotonIDProd:PhotonCutBasedIDTight"
        )
    ),
    photonIsoDepTags = cms.PSet(
        photons = cms.PSet(
            chargedHadron = cms.string('phPFIsoValueCharged03PFIdPFIso'),
            neutralHadron = cms.string('phPFIsoValueNeutral03PFIdPFIso'),
            photon = cms.string('phPFIsoValueGamma03PFIdPFIso')
        )
    ),
    electronIsoDepTags = cms.PSet(
        gsfElectrons = cms.PSet(
            chargedHadron = cms.string('elPFIsoValueCharged03PFIdPFIso'),
            neutralHadron = cms.string('elPFIsoValueNeutral03PFIdPFIso'),
            photon = cms.string('elPFIsoValueGamma03PFIdPFIso')
        )
    ),
    metCollectionTags = cms.vstring(
        "pfMet",
        "genMetTrue",
        "pfType1CorrectedMet",
        "pfType1p2CorrectedMet",
        "pfType1CorrectedMetsysShiftCorr",
        "pfType1p2CorrectedMetsysShiftCorr"
    ),
    bTagCollectionTags = cms.vstring(
        "newTrackCountingHighEffBJetTags",
        "newTrackCountingHighPurBJetTags",
        "newJetProbabilityBJetTags",
        "newJetBProbabilityBJetTags",
        "newSimpleSecondaryVertexBJetTags",
        "newCombinedSecondaryVertexBJetTags",
        "newCombinedSecondaryVertexMVABJetTags",
        "newSoftElectronBJetTags",
        "newSoftMuonBJetTags"
    ),
    puJetIdTags = cms.PSet(
        ak5PFJets = cms.vstring(
            "recoPuJetMva:full",
            "recoPuJetMva:cutbased",
            "recoPuJetMva:simple"
        )
    ),
    genCollectionTag = cms.string("genParticles"),
    puSummaryInfoTag = cms.string("addPileupInfo"),
    triggerEventTag = cms.string("hltTriggerSummaryAOD"),
    photonSCRegressionWeights = cms.string("/afs/cern.ch/user/b/bendavid/cmspublic/regweights52xV3/gbrv3ph_52x.root"),
    muonThreshold = cms.double(2.0),
    electronThreshold = cms.double(2.0),
    photonThreshold = cms.double(10.0),
    jetThreshold = cms.double(20.0),
    pfParticleThreshold = cms.double(1.0),
    debugLevel = cms.int32(0),
    storeL1Info = cms.bool(True),
    storeHLTInfo = cms.bool(True),
    storeGenInfo = cms.bool(True),
    storeGeneralTracks = cms.bool(False),
    storePFJetPartonMatches = cms.bool(True),
    storeTriggerEvents = cms.bool(False),
    recoMode = cms.bool(False),
    isFastSim = cms.bool(False),
    outputFileName = cms.string("susyEvents.root"),
    triggerFileName = cms.string("susyTriggers.root")
)
