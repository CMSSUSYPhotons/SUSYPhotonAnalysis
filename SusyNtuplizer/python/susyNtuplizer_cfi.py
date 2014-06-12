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
        "pfType01CorrectedMet",
        "pfType01p2CorrectedMet",
        "pfSysShiftCorrectedMet",
        "pfType01SysShiftCorrectedMet",
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
        ),
        ak5PFchsJets = cms.PSet(
            TrackCountingHighEff = cms.string("chsTrackCountingHighEffBJetTags"),
            TrackCountingHighPur = cms.string("chsTrackCountingHighPurBJetTags"),
            JetProbability = cms.string("chsJetProbabilityBJetTags"),
            JetBProbability = cms.string("chsJetBProbabilityBJetTags"),
            SimpleSecondaryVertex = cms.string("chsSimpleSecondaryVertexBJetTags"),
            CombinedSecondaryVertex = cms.string("chsCombinedSecondaryVertexBJetTags"),
            CombinedSecondaryVertexMVA = cms.string("chsCombinedSecondaryVertexMVABJetTags"),
            SoftElectron = cms.string("chsSoftElectronBJetTags"),
            SoftMuon = cms.string("chsSoftMuonBJetTags")
        )
    ),
    qgTagCollectionTags = cms.PSet(
        ak5PFJets = cms.string("QGTaggerAK5"),
        ak5PFchsJets = cms.string("QGTaggerAK5chs")
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
    metFilters = cms.PSet(
        CSCBeamHalo = cms.PSet(
            tag = cms.string("BeamHaloSummary"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        HcalNoise = cms.PSet(
            tag = cms.string("HBHENoiseFilterResultProducer:HBHENoiseFilterResult"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        EcalDeadCellTP = cms.PSet(
            tag = cms.string("EcalDeadCellTriggerPrimitiveFilter"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        EcalDeadCellBE = cms.PSet(
            tag = cms.string("EcalDeadCellBoundaryEnergyFilter"),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        TrackingFailure = cms.PSet(
            tag = cms.string("trackingFailureFilter"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        EEBadSC = cms.PSet(
            tag = cms.string("eeBadScFilter"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        HcalLaserOccupancy = cms.PSet(
            tag = cms.string("hcalLaserEventFilter"),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        HcalLaserEventList = cms.PSet(
            tag = cms.string(""),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        HcalLaserRECOUserStep = cms.PSet(
            tag = cms.string(""),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        EcalLaserCorr = cms.PSet(
            tag = cms.string("ecalLaserCorrFilter"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        ManyStripClus53X = cms.PSet(
            tag = cms.string("manystripclus53X"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        TooManyStripClus53X = cms.PSet(
            tag = cms.string("toomanystripclus53X"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        LogErrorTooManyClusters = cms.PSet(
            tag = cms.string("logErrorTooManyClusters"),
            run = cms.bool(True),
            default = cms.bool(True)
        ),
        LogErrorTooManyTripletsPairs = cms.PSet(
            tag = cms.string("logErrorTooManyTripletsPairs"),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        LogErrorTooManySeeds = cms.PSet(
            tag = cms.string("logErrorTooManySeeds"),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        EERingOfFire = cms.PSet(
            tag = cms.string("eeNoiseFilter"),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        InconsistentMuon = cms.PSet(
            tag = cms.string("inconsistentMuonPFCandidateFilter"),
            run = cms.bool(True),
            default = cms.bool(False)
        ),
        GreedyMuon = cms.PSet(
            tag = cms.string("greedyMuonPFCandidateFilter"),
            run = cms.bool(True),
            default = cms.bool(False)
        )
    ),
    hltFilterTags = cms.vstring(),
    photonSCRegressionWeights = cms.FileInPath("SUSYPhotonAnalysis/SusyNtuplizer/data/gbrv3ph_52x.root"),
    muonThreshold = cms.double(2.0),
    electronThreshold = cms.double(2.0),
    photonThreshold = cms.double(10.0),
    jetThreshold = cms.double(20.0),
    pfParticleThreshold = cms.double(3.0),
    genParticleThreshold = cms.double(2.0),
    debugLevel = cms.int32(0),
    selectEvents = cms.vstring(),
    storeL1Info = cms.bool(True),
    storeHLTInfo = cms.bool(True),
    storeGenInfo = cms.bool(True),
    storeGeneralTracks = cms.bool(False),
    storePFJetPartonMatches = cms.bool(True),
    storeTriggerEvents = cms.bool(False),
    storeLumiInfo = cms.bool(True),
    isFastSim = cms.bool(False),
    outputFileName = cms.string("susyEvents.root"),
    triggerFileName = cms.string("susyTriggers.root")
)
