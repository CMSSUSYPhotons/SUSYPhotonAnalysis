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
        "ak5PFJets",
        "ak5PFchsJets"
    ),
    jptJetCollectionTags = cms.vstring(),
    electronMVAIdTags = cms.PSet(
        gsfElectrons = cms.PSet(
            triggering = cms.string("mvaTrigV0"),
            nonTriggering = cms.string("mvaNonTrigV0")
        )
    ),
    muonIdTags = cms.PSet(
        muons = cms.PSet(
            TMLastStationLoose = cms.string("muons:muidTMLastStationLoose"),
            TMLastStationTight = cms.string("muons:muidTMLastStationTight"),
            TMOneStationLoose = cms.string("muons:muidTMOneStationLoose"),
            TMOneStationTight = cms.string("muons:muidTMOneStationTight"),
            TMLastStationLowPtLoose = cms.string("muons:muidTMLastStationOptimizedLowPtLoose"),
            TMLastStationLowPtTight = cms.string("muons:muidTMLastStationOptimizedLowPtTight")
        )
    ),
    photonIsoDepTags = cms.PSet(
        photons = cms.PSet(
            chargedHadron = cms.string('photonPFIsoValueCharged03'),
            neutralHadron = cms.string('photonPFIsoValueNeutral03'),
            photon = cms.string('photonPFIsoValueGamma03')
        )
    ),
    electronIsoDepTags = cms.PSet(
        gsfElectrons = cms.PSet(
            chargedHadron = cms.string('gsfElectronPFIsoValueCharged03'),
            neutralHadron = cms.string('gsfElectronPFIsoValueNeutral03'),
            photon = cms.string('gsfElectronPFIsoValueGamma03')
        )
    ),
    metCollectionTags = cms.vstring(
        "genMetTrue",
        "met",
        "corMetGlobalMuons",
        "caloType1CorrectedMet",
        "caloType1p2CorrectedMet",
        "tcMet",
        "pfMet",
        "pfType1CorrectedMet",
        "pfType1p2CorrectedMet",
        "pfType0CorrectedMet",
        "pfSysShiftCorrectedMet",
        "pfType01SysShiftCorrectedMet",
        "pfType01p2SysShiftCorrectedMet",
        "pfNoPileUpMet",
        "pfMVAMet"
    ),
    bTagCollectionTags = cms.PSet(
        ak5PFJets = cms.PSet(
            TrackCountingHighEff = cms.string("newTrackCountingHighEffBJetTags"),
            TrackCountingHighPur = cms.string("newTrackCountingHighPurBJetTags"),
            JetProbability = cms.string("newJetProbabilityBJetTags"),
            JetBProbability = cms.string("newJetBProbabilityBJetTags"),
            SimpleSecondaryVertex = cms.string("newSimpleSecondaryVertexBJetTags"),
            CombinedSecondaryVertex = cms.string("newCombinedSecondaryVertexBJetTags"),
            CombinedSecondaryVertexMVA = cms.string("newCombinedSecondaryVertexMVABJetTags"),
            SoftElectron = cms.string("newSoftElectronBJetTags"),
            SoftMuon = cms.string("newSoftMuonBJetTags")
        )
    ),
    jetFlavourMatchingTags = cms.PSet(
        ak5PFJets = cms.PSet(
            alg = cms.string("flavourAssociationAlg"),
            phy = cms.string("flavourAssociationPhy")
        ),
        ak5PFchsJets = cms.PSet(
            alg = cms.string("flavourAssociationCHSAlg"),
            phy = cms.string("flavourAssociationCHSPhy")
        )
    ),
    puJetId = cms.PSet(
        ak5PFJets = cms.PSet(
            tag = cms.string("recoPuJetId"),
            algorithms = cms.vstring("full", "cutbased")
        ),
        ak5PFchsJets = cms.PSet(
            tag = cms.string("recoPuJetIdChs"),
            algorithms = cms.vstring("full", "cutbased")
        )
    ),
    pfPUCandidatesTag = cms.string("pfPileUp"),
    genCollectionTag = cms.string("genParticles"),
    puSummaryInfoTag = cms.string("addPileupInfo"),
    triggerEventTag = cms.string("hltTriggerSummaryAOD"),
    gridParams = cms.vstring("ptHat"),
    photonSCRegressionWeights = cms.string("/afs/cern.ch/user/b/bendavid/cmspublic/regweights52xV3/gbrv3ph_52x.root"),
    muonThreshold = cms.double(2.0),
    electronThreshold = cms.double(2.0),
    photonThreshold = cms.double(10.0),
    jetThreshold = cms.double(20.0),
    pfParticleThreshold = cms.double(3.0),
    debugLevel = cms.int32(0),
    storeL1Info = cms.bool(True),
    storeHLTInfo = cms.bool(True),
    storeGenInfo = cms.bool(True),
    storeGeneralTracks = cms.bool(False),
    storePFJetPartonMatches = cms.bool(True),
    storeTriggerEvents = cms.bool(False),
    isFastSim = cms.bool(False),
    outputFileName = cms.string("susyEvents.root"),
    triggerFileName = cms.string("susyTriggers.root")
)
