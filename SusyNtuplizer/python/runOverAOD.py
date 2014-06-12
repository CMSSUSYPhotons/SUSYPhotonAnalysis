def configure( dataset, sourceNames=[], hltPaths=[], maxEvents = -1, outputName = '', runNoPUMVAMetSequence=False ):
    """ This is the main routine used to configure the nTuplizer. Please check
    the correctness of all relevant default parameters.

    dataset: Please choose one of the entries in the lists.
    sourceNames: list of AOD filenames
    hltPaths: list of trigger hlt paths (will be ignored for MC)
    maxEvents: maximal number of events
    runNoPUMVAMetSequence: please read README for more details
    """
    collisionDatasets = [
        '52xPrompt',     # Run2012[AB]{C}-PromptReco{_v1}
        '52x23May2012',  # Run2012A-23May2012
        '53xPromptC',    # Run2012C-PromptReco-v2
        '53xPromptD',    # Run2012D-PromptReco
        '53x13July2012', # Run2012[AB]-13Jul2012
        '53x06Aug2012',  # Run2012A-06Aug2012
        '53x24Aug2012',  # Run2012C-24Aug2012
        '53x11Dec2012',  # Run2012C-EcalRecover_11Dec2012
        '53x16Jan2013',  # Run2012D-16Jan2013
        '53x22Jan2013'   # Run2012[ABCD]-22Jan2013
    ]
    mcDatasets = [
        '52xFullSim',
        '52xFastSim',
        '53xFullSim',
        '53xFastSim',
        '70xFullSim'
    ]

    isRealData = dataset in collisionDatasets
    isMC = dataset in mcDatasets
    isFastSim = isMC and 'FastSim' in dataset
    is53x = '53x' in dataset
    is52x = '52x' in dataset
    is70x = '70x' in dataset

    if not isRealData and not isMC:
        raise RuntimeError("Dataset " + dataset + " not defined")

    import FWCore.ParameterSet.Config as cms

    ##########################
    ### Initialize process ###
    ##########################
    process = cms.Process("RA3")
    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False),
        FailPath = cms.untracked.vstring(
            'FatalRootError'
        )
    )

    process.source = cms.Source("PoolSource",
        noEventSort = cms.untracked.bool(True),
        duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
        fileNames = cms.untracked.vstring(sourceNames)
    )

    #####################
    ### MessageLogger ###
    #####################
    process.load('FWCore.MessageService.MessageLogger_cfi')
    process.MessageLogger.cerr.FwkReport.reportEvery = 1000
    process.MessageLogger.suppressWarning = cms.untracked.vstring(
        'newSecondaryVertexTagInfos',
        'chsSecondaryVertexTagInfos',
        'pfCandidateToVertexAssoc',
        'manystripclus53X',
        'toomanystripclus53X'
    )
    process.MessageLogger.suppressError = cms.untracked.vstring('ecalLaserCorrFilter')
    process.MessageLogger.categories.append('SusyNtuplizer')
    process.MessageLogger.cerr.SusyNtuplizer = cms.untracked.PSet( limit = cms.untracked.int32(100) )

    #############################
    ### Conditions & Services ###
    #############################
    process.load('Configuration.Geometry.GeometryIdeal_cff')
    process.load('Configuration.StandardSequences.MagneticField_38T_cff')
    process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

    # The following global tags are the latest as of April 8th 2013. The user is strongly recommended to check the
    # SWGuideFrontierConditions twiki before running the configuration.
    if dataset == '52xPrompt':
        process.GlobalTag.globaltag = 'GR_P_V39_AN3::All'
    elif dataset == '52x23May2012':
        process.GlobalTag.globaltag = 'FT_P_V32B_AN4::All'
    elif dataset == '53xPromptC' or dataset == '53xPromptD':
        process.GlobalTag.globaltag = 'GR_P_V42_AN4::All'
    elif dataset == '53x13July2012':
        process.GlobalTag.globaltag = 'FT_53_V6C_AN4::All'
    elif dataset == '53x06Aug2012':
        process.GlobalTag.globaltag = 'FT_53_V6C_AN4::All'
    elif dataset == '53x24Aug2012':
        process.GlobalTag.globaltag = 'FT53_V10A_AN4::All'
    elif dataset == '53x11Dec2012':
        process.GlobalTag.globaltag = 'FT_P_V42C_AN4::All'
    elif dataset == '53x16Jan2013':
        process.GlobalTag.globaltag = 'FT_P_V43E_AN4::All'
    elif dataset == '53x22Jan2013':
        process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'
    elif isMC and is52x:
        process.GlobalTag.globaltag = 'START52_V16::All'
    elif isMC and is53x:
        process.GlobalTag.globaltag = 'START53_V25::All'
    elif isMC and is70x:
        process.GlobalTag.globaltag = 'POSTLS170_V6::All'

    #####################
    ### SusyNtuplizer ###
    #####################
    process.load("SUSYPhotonAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
    process.susyNtuplizer.debugLevel = 0
    process.susyNtuplizer.isFastSim = isFastSim
    process.susyNtuplizer.caloJetCollectionTags = []
    if outputName:
        process.susyNtuplizer.outputFileName = outputName

    #########################
    ### HLT result filter ###
    #########################
    if len(hltPaths) != 0:
        process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
        process.hltHighLevel.HLTPaths = hltPaths
    else:
        process.hltHighLevel = cms.Sequence() # placeholder

    ###################################
    ### HCAL laser pollution filter ###
    ###################################
    # A fraction of 2012 data had HCAL lasers firing coincident with bunch crossing
    # Unlike other MET filters (see below), the filter provided by the HCAL DPG does
    # not allow to run on "tagging mode", i.e. saving the result of the filter but
    # not actually filtering out the event.
    process.load('EventFilter.HcalRawToDigi.hcallasereventfilter2012_cfi')

    ##########################################
    ### Good vertex collection (transient) ###
    ##########################################
    process.goodVertices = cms.EDFilter("VertexSelector",
        filter = cms.bool(False),
        src = cms.InputTag("offlinePrimaryVertices"),
        cut = cms.string("!isFake && ndof >= 4 && abs(z) <= 24 && position.rho < 2")
    )
    process.primaryVertex = cms.EDFilter("PATSingleVertexSelector",
        mode = cms.string('firstVertex'),
        vertices = cms.InputTag('goodVertices'),
        filter = cms.bool(False)
    )

    process.vertexSelectionSequence = cms.Sequence(
        process.goodVertices +
        process.primaryVertex
    )

    ###########################
    ### PU Rho calculations ###
    ###########################
    from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets

    # 70X does not come with kt PF jets any more - rho correction will not be used in the future..
    process.kt6PFJets = kt4PFJets.clone(
        rParam = cms.double(0.6),
        doRhoFastjet = cms.bool(True),
        doAreaFastjet = cms.bool(True),
        voronoiRfact = cms.double(0.9)
    )

    #Calculate rho restricted to barrel for photon pileup subtraction
    process.kt6PFJetsRhoBarrelOnly = kt4PFJets.clone(
        src = cms.InputTag('particleFlow'),
        rParam = cms.double(0.6),
        #Eta range of jets to be considered for Rho calculation
        #Should be at most (jet acceptance - jet radius)
        doRhoFastjet = cms.bool(True),
        Rho_EtaMax=cms.double(1.4442),
        #Eta range of ghost jets to be considered for Rho calculation - must be greater than Rho_EtaMax
        Ghost_EtaMax=cms.double(2.5)
    )

    #Calculate rho restricted to barrel for photon pileup subtraction
    process.kt6PFJetsRho25 = kt4PFJets.clone(
        src = cms.InputTag('particleFlow'),
        rParam = cms.double(0.6),
        #Eta range of jets to be considered for Rho calculation
        #Should be at most (jet acceptance - jet radius)
        doRhoFastjet = cms.bool(True),
        Rho_EtaMax=cms.double(2.5)
    )

    process.puRhoSequence = cms.Sequence(
        process.kt6PFJets +
        process.kt6PFJetsRhoBarrelOnly +
        process.kt6PFJetsRho25
    )

    ###############################
    ### PF-based reconstruction ###
    ###############################
    # Run the full event reconstruction starting from the full list of PFCandidates (particleFlow).
    # Serves as inputs to isoDeposit, CHS jets, and NoPUMET.
    process.load("CommonTools.ParticleFlow.PFBRECO_cff")

    process.pfPileUp.Vertices = cms.InputTag("goodVertices")
    process.pfPileUp.checkClosestZVertex = False

    process.pfBasedRecoSequence = cms.Sequence(
        process.particleFlowPtrs +
        process.pfNoPileUpSequence +
        process.pfParticleSelectionSequence +
        process.pfPhotonSequence +
        process.pfMuonSequence +
        process.pfNoMuon +
        process.pfElectronSequence +
        process.pfNoElectron
    )

    ####################################
    ### Photon & electron isoDeposit ###
    ####################################
    # pf isolation tools seems broken in 70X (Jun 12 2014, Y.I.)
#    import CommonTools.ParticleFlow.Tools.pfIsolation as pfIsolation
#    process.load('CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff')
#    process.load('CommonTools.ParticleFlow.Isolation.pfPhotonIsolation_cff')
#
#    process.eleIsoSequence = pfIsolation.setupPFElectronIso(process, 'gedGsfElectrons', newpostfix = 'GedEle')
#    process.phoIsoSequence = pfIsolation.setupPFPhotonIso(process, 'photons', newpostfix = 'OldPh')
#    process.gedPhoIsoSequence = pfIsolation.setupPFPhotonIso(process, 'gedPhotons', newpostfix = 'GedPh')

    process.pfIsolationSequence = cms.Sequence(
#        process.eleIsoSequence +
#        process.phoIsoSequence +
#        process.gedPhoIsoSequence
    )

    ###############################################
    ### PFchs (charged hadron subtraction) jets ###
    ###############################################
    # pfJets comes from PFBRECO.

    process.ak5PFJetsCHS = process.pfJets.clone()
    process.ak5PFJetsCHS.doAreaFastjet = True

    process.pfCHSJetSequence = cms.Sequence(
        process.ak5PFJetsCHS
    )

    ###########
    ### JEC ###
    ###########
    process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
    process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")

    # Setting up JEC ESProducers for ak5PFchs. This block will be included in the JetCorrectionServices_cff
    # in a future tag by JetMET (April 8, 2013)
    process.ak5PFchsL1Fastjet = process.ak5PFL1Fastjet.clone(algorithm = cms.string('AK5PFchs'))
    process.ak5PFchsL2Relative = process.ak5PFL2Relative.clone(algorithm = cms.string('AK5PFchs'))
    process.ak5PFchsL3Absolute = process.ak5PFL3Absolute.clone(algorithm = cms.string('AK5PFchs'))
    process.ak5PFchsResidual = process.ak5PFResidual.clone(algorithm = cms.string('AK5PFchs'))
    process.ak5PFchsL2L3 = cms.ESProducer('JetCorrectionESChain',
        correctors = cms.vstring('ak5PFchsL2Relative', 'ak5PFchsL3Absolute')
    )
    process.ak5PFchsL2L3Residual = process.ak5PFchsL2L3.clone()
    process.ak5PFchsL2L3Residual.correctors.append('ak5PFchsResidual')
    process.ak5PFchsL1FastL2L3 = process.ak5PFchsL2L3.clone()
    process.ak5PFchsL1FastL2L3.correctors.insert(0, 'ak5PFchsL1Fastjet')
    process.ak5PFchsL1FastL2L3Residual = process.ak5PFchsL1FastL2L3.clone()
    process.ak5PFchsL1FastL2L3Residual.correctors.append('ak5PFchsResidual')

    ###############################
    ### Corrected MET producers ###
    ###############################
    process.load("JetMETCorrections.Type1MET.caloMETCorrections_cff")
    process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

    if isRealData:
        process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3Residual")
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
    else:
        process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3")
        process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")

    from JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi import pfMEtSysShiftCorr
    process.pfMEtSysShiftCorr = pfMEtSysShiftCorr.clone(
        srcVertices = cms.InputTag('goodVertices')
    )

    # http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/JetMETCorrections/Type1MET/python/pfMETsysShiftCorrections_cfi.py?revision=1.6&view=markup
    # pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data (_mc)
    if isRealData:
        process.pfMEtSysShiftCorr.parameter = cms.PSet(
            px = cms.string("+0.2661 + 0.3217*Nvtx"),
            py = cms.string("-0.2251 - 0.1747*Nvtx")
        )
    else:
        process.pfMEtSysShiftCorr.parameter = cms.PSet(
            px = cms.string("+0.1166 + 0.0200*Nvtx"),
            py = cms.string("+0.2764 - 0.1280*Nvtx")
        )

    # Remove the track-based type 0 correction
    process.producePFMETCorrections.remove(process.pfchsMETcorr)

    # Use PF-based type 0 correction. Reference: AN-2012/333
    # Note that type 0 and type 1 are not truly orthogonal; the implicit assumption is that the PU
    # contributions are fully captured in the high-pt jet offset, low-pt jet, and non-clustering
    # energy terms in the type 1 MET correction formula.
    from CommonTools.RecoUtils.pfcand_assomap_cfi import PFCandAssoMap
    from JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi import pfMETcorrType0

    # Full implementation requires the reconstruction of displaced vertex, which we are skipping here.
    # (Displaced vertex reconstruction requires TrackExtra information and is impossible from AOD)
    process.pfCandidateToVertexAssoc = PFCandAssoMap.clone()

    # Attention: Type0PFMETcorrInputProducer does not produce sumEt of the correction.
    process.pfType0MetCorrection = pfMETcorrType0.clone(
        srcPFCandidateToVertexAssociations = cms.InputTag('pfCandidateToVertexAssoc'),
        srcHardScatterVertex = cms.InputTag('primaryVertex')
    )

    process.pfType01CorrectedMet = process.pfType1CorrectedMet.clone(
        srcType1Corrections = cms.VInputTag(
            cms.InputTag('pfType0MetCorrection'),
            cms.InputTag('pfJetMETcorr', 'type1')
        )
    )
    process.pfType01p2CorrectedMet = process.pfType1p2CorrectedMet.clone(
        srcType1Corrections = cms.VInputTag(
            cms.InputTag('pfType0MetCorrection'),
            cms.InputTag('pfJetMETcorr', 'type1')
        )
    )
    process.pfSysShiftCorrectedMet = process.pfType1CorrectedMet.clone(
        srcType1Corrections = cms.VInputTag(
            cms.InputTag('pfMEtSysShiftCorr')
        )
    )
    process.pfType01SysShiftCorrectedMet = process.pfType1CorrectedMet.clone(
        srcType1Corrections = cms.VInputTag(
            cms.InputTag('pfType0MetCorrection'),
            cms.InputTag('pfJetMETcorr', 'type1'),
            cms.InputTag('pfMEtSysShiftCorr')
        )
    )

    process.correctedMetSequence = cms.Sequence(
        process.produceCaloMETCorrections +
        process.producePFMETCorrections +
        process.pfCandidateToVertexAssoc +
        process.pfType0MetCorrection +
        process.pfType01CorrectedMet +
        process.pfType01p2CorrectedMet +
        process.pfMEtSysShiftCorr +
        process.pfSysShiftCorrectedMet +
        process.pfType01SysShiftCorrectedMet
    )

    #############################
    ### MVA-based electron ID ###
    #############################
    # ELECTRON MVA NOT MEANINGFUL FOR A WHILE
#    process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(
#        process.mvaTrigV0 *
#        process.mvaNonTrigV0
    )

    #################
    ### PU jet ID ###
    #################
#    from RecoJets.JetProducers.PileupJetIDParams_cfi import full_5x_chs, cutbased
#    from RecoJets.JetProducers.PileupJetID_cfi import pileupJetIdProducer, pileupJetIdProducerChs
#
#    process.recoPuJetIdChs = pileupJetIdProducerChs.clone(
#        jets = cms.InputTag("ak5PFJetsCHS"),
#        algos = cms.VPSet(full_5x_chs, cutbased),
#        applyJec = cms.bool(True),
#        inputIsCorrected = cms.bool(False),
#        residualsTxt = cms.FileInPath("SUSYPhotonAnalysis/SusyNtuplizer/python/runOverAOD.py") # not used, but has to point to an existing file
#    )

    process.recoPuJetIdSequence = cms.Sequence(
#        process.recoPuJetIdChs
    )

    ##################################
    ### Jet-parton matching for MC ###
    ##################################
    if isMC:
        process.myPartons = cms.EDProducer("PartonSelector",
            withLeptons = cms.bool(False),
            src = cms.InputTag("genParticles")
        )
        process.flavourByRef = cms.EDProducer("JetPartonMatcher",
            jets = cms.InputTag("ak5PFJets"),
            coneSizeToAssociate = cms.double(0.3),
            partons = cms.InputTag("myPartons")
        )
        process.flavourAssociationAlg = cms.EDProducer("JetFlavourIdentifier",
            srcByReference = cms.InputTag("flavourByRef"),
            physicsDefinition = cms.bool(False)
        )
        process.flavourAssociationPhy = cms.EDProducer("JetFlavourIdentifier",
            srcByReference = cms.InputTag("flavourByRef"),
            physicsDefinition = cms.bool(True)
        )
        process.flavourByRefCHS = cms.EDProducer("JetPartonMatcher",
            jets = cms.InputTag("ak5PFJetsCHS"),
            coneSizeToAssociate = cms.double(0.3),
            partons = cms.InputTag("myPartons")
        )
        process.flavourAssociationCHSAlg = cms.EDProducer("JetFlavourIdentifier",
            srcByReference = cms.InputTag("flavourByRefCHS"),
            physicsDefinition = cms.bool(False)
        )
        process.flavourAssociationCHSPhy = cms.EDProducer("JetFlavourIdentifier",
            srcByReference = cms.InputTag("flavourByRefCHS"),
            physicsDefinition = cms.bool(True)
        )
        process.JetFlavourMatchingSequence = cms.Sequence(
            process.myPartons *
            process.flavourByRef *
            process.flavourByRefCHS *
            process.flavourAssociationAlg *
            process.flavourAssociationPhy *
            process.flavourAssociationCHSAlg *
            process.flavourAssociationCHSPhy
        )
    else:
        process.JetFlavourMatchingSequence = cms.Sequence() # placeholder


    #################
    ### b-tagging ###
    #################
    # Re-run b-tagging with PFJets as input

    # b-tagging general configuration
    process.load("RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi")
    process.load("RecoBTag.Configuration.RecoBTag_cff")

    # create a new jets and tracks associaiton
    process.newJetTracksAssociatorAtVertex = process.ic5PFJetTracksAssociatorAtVertex.clone(
        jets = cms.InputTag("ak5PFJets"),
        tracks = cms.InputTag("generalTracks")
        )
    process.chsJetTracksAssociatorAtVertex = process.ic5PFJetTracksAssociatorAtVertex.clone(
        jets = cms.InputTag("ak5PFJetsCHS"),
        tracks = cms.InputTag("generalTracks")
        )

    # impact parameter b-tag
    process.newImpactParameterTagInfos = process.impactParameterTagInfos.clone(
        jetTracks = cms.InputTag("newJetTracksAssociatorAtVertex")
        )
    process.chsImpactParameterTagInfos = process.impactParameterTagInfos.clone(
        jetTracks = cms.InputTag("chsJetTracksAssociatorAtVertex")
        )
    # TCHE
    process.newTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
        )
    process.chsTrackCountingHighEffBJetTags = process.trackCountingHighEffBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsImpactParameterTagInfos") )
        )
    # TCHP
    process.newTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
        )
    process.chsTrackCountingHighPurBJetTags = process.trackCountingHighPurBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsImpactParameterTagInfos") )
        )
    # JP
    process.newJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
        )
    process.chsJetProbabilityBJetTags = process.jetProbabilityBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsImpactParameterTagInfos") )
        )
    # JBP
    process.newJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
        )
    process.chsJetBProbabilityBJetTags = process.jetBProbabilityBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsImpactParameterTagInfos") )
        )

    # secondary vertex b-tag
    process.newSecondaryVertexTagInfos = process.secondaryVertexTagInfos.clone(
        trackIPTagInfos = cms.InputTag("newImpactParameterTagInfos")
        )
    process.chsSecondaryVertexTagInfos = process.secondaryVertexTagInfos.clone(
        trackIPTagInfos = cms.InputTag("chsImpactParameterTagInfos")
        )
    # SSV
    process.newSimpleSecondaryVertexBJetTags = process.simpleSecondaryVertexBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
        )
    process.chsSimpleSecondaryVertexBJetTags = process.simpleSecondaryVertexBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsSecondaryVertexTagInfos") )
        )
    # CSV
    process.newCombinedSecondaryVertexBJetTags = process.combinedSecondaryVertexBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"),
                                  cms.InputTag("newSecondaryVertexTagInfos")
            )
        )
    process.chsCombinedSecondaryVertexBJetTags = process.combinedSecondaryVertexBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsImpactParameterTagInfos"),
                                  cms.InputTag("chsSecondaryVertexTagInfos")
            )
        )
    # CSVMVA
    process.newCombinedSecondaryVertexMVABJetTags = process.combinedSecondaryVertexMVABJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"),
                                  cms.InputTag("newSecondaryVertexTagInfos")
            )
        )
    process.chsCombinedSecondaryVertexMVABJetTags = process.combinedSecondaryVertexMVABJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsImpactParameterTagInfos"),
                                  cms.InputTag("chsSecondaryVertexTagInfos")
            )
        )
    # soft electron b-tag
    process.newSoftPFElectronsTagInfos = process.softPFElectronsTagInfos.clone(
        jets = "ak5PFJets"
        )
    process.newSoftPFElectronBJetTags = process.softPFElectronBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newSoftPFElectronsTagInfos") )
        )
    process.chsSoftPFElectronsTagInfos = process.softPFElectronsTagInfos.clone(
        jets = "ak5PFJetsCHS"
        )
    process.chsSoftPFElectronBJetTags = process.softPFElectronBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsSoftPFElectronsTagInfos") )
        )
    # soft muon b-tag
    process.newSoftPFMuonsTagInfos = process.softPFMuonsTagInfos.clone(
        jets = "ak5PFJets"
        )
    process.newSoftPFMuonBJetTags = process.softPFMuonBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newSoftPFMuonsTagInfos") )
        )
    process.newSoftPFMuonByIP3dBJetTags = process.softPFMuonByIP3dBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newSoftPFMuonsTagInfos") )
        )
    process.newSoftPFMuonByPtBJetTags = process.softPFMuonByPtBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("newSoftPFMuonsTagInfos") )
        )
    process.chsSoftPFMuonsTagInfos = process.softPFMuonsTagInfos.clone(
        jets = "ak5PFJetsCHS"
        )
    process.chsSoftPFMuonBJetTags = process.softPFMuonBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsSoftPFMuonsTagInfos") )
        )
    process.chsSoftPFMuonByIP3dBJetTags = process.softPFMuonByIP3dBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsSoftPFMuonsTagInfos") )
        )
    process.chsSoftPFMuonByPtBJetTags = process.softPFMuonByPtBJetTags.clone(
        tagInfos = cms.VInputTag( cms.InputTag("chsSoftPFMuonsTagInfos") )
        )

    process.newJetTracksAssociator = cms.Sequence(
        process.newJetTracksAssociatorAtVertex
    )
    process.chsJetTracksAssociator = cms.Sequence(
        process.chsJetTracksAssociatorAtVertex
    )

    process.newJetBtaggingIP = cms.Sequence(
        process.newImpactParameterTagInfos * (
            process.newTrackCountingHighEffBJetTags +
            process.newTrackCountingHighPurBJetTags +
            process.newJetProbabilityBJetTags +
            process.newJetBProbabilityBJetTags
        )
    )
    process.chsJetBtaggingIP = cms.Sequence(
        process.chsImpactParameterTagInfos * (
            process.chsTrackCountingHighEffBJetTags +
            process.chsTrackCountingHighPurBJetTags +
            process.chsJetProbabilityBJetTags +
            process.chsJetBProbabilityBJetTags
        )
    )

    process.newJetBtaggingSV = cms.Sequence(
        process.newImpactParameterTagInfos *
        process.newSecondaryVertexTagInfos * (
            process.newSimpleSecondaryVertexBJetTags +
            process.newCombinedSecondaryVertexBJetTags +
            process.newCombinedSecondaryVertexMVABJetTags
        )
    )
    process.chsJetBtaggingSV = cms.Sequence(
        process.chsImpactParameterTagInfos *
        process.chsSecondaryVertexTagInfos * (
            process.chsSimpleSecondaryVertexBJetTags +
            process.chsCombinedSecondaryVertexBJetTags +
            process.chsCombinedSecondaryVertexMVABJetTags
        )
    )

    process.newJetBtaggingEle = cms.Sequence(
        process.newSoftPFElectronsTagInfos *
        process.newSoftPFElectronBJetTags
    )
    process.chsJetBtaggingEle = cms.Sequence(
        process.chsSoftPFElectronsTagInfos *
        process.chsSoftPFElectronBJetTags
    )

    process.newJetBtaggingMu = cms.Sequence(
        process.newSoftPFMuonsTagInfos * (
            process.newSoftPFMuonBJetTags +
            process.newSoftPFMuonByIP3dBJetTags +
            process.newSoftPFMuonByPtBJetTags
        )
    )
    process.chsJetBtaggingMu = cms.Sequence(
        process.chsSoftPFMuonsTagInfos * (
            process.chsSoftPFMuonBJetTags +
            process.chsSoftPFMuonByIP3dBJetTags +
            process.chsSoftPFMuonByPtBJetTags
        )
    )

    process.newJetBtagging = cms.Sequence(
        process.newJetBtaggingIP +
        process.newJetBtaggingSV +
        process.newJetBtaggingEle +
        process.newJetBtaggingMu
    )
    process.chsJetBtagging = cms.Sequence(
        process.chsJetBtaggingIP +
        process.chsJetBtaggingSV +
        process.chsJetBtaggingEle +
        process.chsJetBtaggingMu
    )

    # These are fixes for the JetProbability b-tagger calibrations as recommended by BTV.
    # See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagJetProbabilityCalibration#Calibration_in_52x_Data_and_MC
    useSpecialBTagCalibration = True
    if isRealData and is52x:
        btagTrackProbability2DTag = "TrackProbabilityCalibration_2D_2012DataTOT_v1_offline"
        btagTrackProbability3DTag = "TrackProbabilityCalibration_3D_2012DataTOT_v1_offline"
    elif isRealData and is53x:
        btagTrackProbability2DTag = "TrackProbabilityCalibration_2D_Data53X_v2"
        btagTrackProbability3DTag = "TrackProbabilityCalibration_3D_Data53X_v2"
    elif isMC and is53x:
        btagTrackProbability2DTag = "TrackProbabilityCalibration_2D_MC53X_v2"
        btagTrackProbability3DTag = "TrackProbabilityCalibration_3D_MC53X_v2"
    else:
        useSpecialBTagCalibration = False

    if useSpecialBTagCalibration:
        process.GlobalTag.toGet = cms.VPSet(
            cms.PSet(record = cms.string("BTagTrackProbability2DRcd"),
                tag = cms.string(btagTrackProbability2DTag),
                connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")
            ),
            cms.PSet(record = cms.string("BTagTrackProbability3DRcd"),
                tag = cms.string(btagTrackProbability3DTag),
                connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_BTAU")
            )
        )

    ###########################
    ### Quark-gluon tagging ###
    ###########################
#    from QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff import goodOfflinePrimaryVerticesQG, QGTagger
#
#    process.goodOfflinePrimaryVerticesQG = goodOfflinePrimaryVerticesQG
#
#    process.QGTaggerAK5 = QGTagger.clone(
#        srcJets = cms.InputTag('ak5PFJets'),
#        srcRho = cms.InputTag('kt6PFJets', 'rho'),
#        srcRhoIso = cms.InputTag('kt6PFJetsRho25', 'rho')
#    )
#    process.QGTaggerAK5chs = QGTagger.clone(
#        srcJets = cms.InputTag('ak5PFJetsCHS'),
#        useCHS = cms.untracked.bool(True),
#        srcRho = cms.InputTag('kt6PFJets', 'rho'),
#        srcRhoIso = cms.InputTag('kt6PFJetsRho25', 'rho')
#    )
#
#    if isRealData:
#        process.QGTaggerAK5.jec = cms.untracked.string('ak5PFL1FastL2L3Residual')
#        process.QGTaggerAK5chs.jec = cms.untracked.string('ak5PFchsL1FastL2L3Residual')
#    else:
#        process.QGTaggerAK5.jec = cms.untracked.string('ak5PFL1FastL2L3')
#        process.QGTaggerAK5chs.jec = cms.untracked.string('ak5PFchsL1FastL2L3')

    process.QGTaggingSequence = cms.Sequence(
#        process.goodOfflinePrimaryVerticesQG +
#        process.QGTaggerAK5 +
#        process.QGTaggerAK5chs
    )

    ###################
    ### MET filters ###
    ###################
    # HBHENoiseFilterResultProducer
    process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')

    # HCAL laser events filter (for 2011 data and 2012D prompt)
    process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
    process.hcalLaserEventFilter.taggingMode = cms.bool(True)

    # EcalDeadCellTriggerPrimitiveFilter
    process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
    process.EcalDeadCellTriggerPrimitiveFilter.taggingMode = cms.bool(True)

    # EcalDeadCellBoundaryEnergyFilter
    process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
    process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(True)

    # Tracking failure filter
    process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
    process.trackingFailureFilter.taggingMode = cms.bool(True)

    # EE Bad SC Filter
    process.load('RecoMET.METFilters.eeBadScFilter_cfi')
    process.eeBadScFilter.taggingMode = cms.bool(True)

    # EE ring of fire
    process.load('RecoMET.METFilters.eeNoiseFilter_cfi')
    process.eeNoiseFilter.taggingMode = cms.bool(True)

    # Inconsistent muon pf candidate filter
    process.load('RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi')
    process.inconsistentMuonPFCandidateFilter.taggingMode = cms.bool(True)

    # Greedy muon pf candidate filter
    process.load('RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi')
    process.greedyMuonPFCandidateFilter.taggingMode = cms.bool(True)

    # The ECAL laser correction filter (needs correct GT to work)
    process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
    process.ecalLaserCorrFilter.taggingMode = cms.bool(True)

    # The tracking POG filters
    process.load('RecoMET.METFilters.trackingPOGFilters_cff')
    process.manystripclus53X.taggedMode = cms.untracked.bool(True)
    process.toomanystripclus53X.taggedMode = cms.untracked.bool(True)
    process.logErrorTooManyClusters.taggedMode = cms.untracked.bool(True)
    process.logErrorTooManyTripletsPairs.taggedMode = cms.untracked.bool(True)
    process.logErrorTooManySeeds.taggedMode = cms.untracked.bool(True)

    #Add up all MET filters
    process.metFiltersSequence = cms.Sequence(
        process.HBHENoiseFilterResultProducer +
        process.hcalLaserEventFilter +
        process.EcalDeadCellTriggerPrimitiveFilter +
        process.EcalDeadCellBoundaryEnergyFilter +
        process.trackingFailureFilter +
        process.eeBadScFilter +
        process.eeNoiseFilter +
        process.inconsistentMuonPFCandidateFilter +
        process.greedyMuonPFCandidateFilter +
        process.ecalLaserCorrFilter +
        process.manystripclus53X +
        process.toomanystripclus53X +
        process.logErrorTooManyClusters +
        process.logErrorTooManyTripletsPairs +
        process.logErrorTooManySeeds
    )

    ######################
    ### NoPU & MVA MET ###
    ######################

    if runNoPUMVAMetSequence:
        from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams
        from JetMETCorrections.Configuration.JetCorrectionProducers_cff import ak5PFJetsL1
        if isRealData:
            process.ak5PFJetsL123Corrected = ak5PFJetsL1.clone(correctors = ['ak5PFL1FastL2L3Residual'])
        else:
            process.ak5PFJetsL123Corrected = ak5PFJetsL1.clone(correctors = ['ak5PFL1FastL2L3'])

        process.recoPuJetIdCorrected = pileupJetIdProducer.clone(
            jets = cms.InputTag("ak5PFJetsL123Corrected"),
            algos = cms.VPSet(full_5x, cutbased),
            residualsTxt = cms.FileInPath("RecoJets/JetProducers/data/mva_JetID_v1.weights.xml") # not used, but has to point to an existing file
        )

        # NoPU MET
        from JetMETCorrections.METPUSubtraction.noPileUpPFMET_cff import noPileUpPFMEtData, noPileUpPFMEt
        process.pfNoPileUpMetData = noPileUpPFMEtData.clone(
            srcJets = cms.InputTag('ak5PFJetsL123Corrected'),
            srcJetIds = cms.InputTag('recoPuJetIdCorrected', 'fullId'),
            srcPFCandToVertexAssociations = cms.InputTag('pfCandidateToVertexAssoc'),
            srcHardScatterVertex = cms.InputTag('primaryVertex')
        )

        process.pfNoPileUpMet = noPileUpPFMEt.clone(
            srcMVAMEtData = cms.InputTag('pfNoPileUpMetData'),
            srcLeptons = cms.VInputTag(
                'pfIsolatedPhotons',
                'pfIsolatedMuons',
                'pfIsolatedElectrons'
            ),
            srcType0Correction = cms.InputTag('pfType0MetCorrection'),
            saveInputs = cms.bool(False)
        )

        # MVA MET
        from JetMETCorrections.METPUSubtraction.mvaPFMET_cff import pfMEtMVA
        process.pfMVAMet = pfMEtMVA.clone(
            srcCorrJets = cms.InputTag('ak5PFJetsL123Corrected'),
            srcLeptons = cms.VInputTag(
                'pfIsolatedPhotons',
                'pfIsolatedMuons',
                'pfIsolatedElectrons'
            )
        )

        process.noPUMVAMetSequence = cms.Sequence(
            process.ak5PFJetsL123Corrected +
            process.recoPuJetIdCorrected +
            process.pfNoPileUpMetData +
            process.pfNoPileUpMet +
            process.pfMVAMet
        )
    else:
        process.noPUMVAMetSequence = cms.Sequence() # placeholder

    #####################
    ### Finalize path ###
    #####################
    process.standard_step = cms.Path(
        process.hltHighLevel +
        process.hcallasereventfilter2012 +
        process.vertexSelectionSequence +
        process.puRhoSequence +
        process.pfBasedRecoSequence +
        process.pfIsolationSequence +
        process.pfCHSJetSequence +
        process.correctedMetSequence +
        process.eidMVASequence +
        process.recoPuJetIdSequence +
        process.JetFlavourMatchingSequence +
        process.newJetTracksAssociator +
        process.chsJetTracksAssociator +
        process.newJetBtagging +
        process.chsJetBtagging +
        process.QGTaggingSequence +
        process.metFiltersSequence
    )

    process.optional_step = cms.Path(
        process.hltHighLevel +
        process.hcallasereventfilter2012 +
        process.noPUMVAMetSequence
    )

    process.ntuplizer_step = cms.Path(
        process.hltHighLevel +
        process.hcallasereventfilter2012 +
        process.susyNtuplizer
    )

    process.schedule = cms.Schedule(process.standard_step,process.optional_step,process.ntuplizer_step)


    ###################################################################
    ### Dataset-dependent sequence and event content configurations ###
    ###################################################################

    if isRealData:
        process.susyNtuplizer.metCollectionTags.remove('genMetTrue')
        process.susyNtuplizer.jetFlavourMatchingTags = cms.PSet()
        process.susyNtuplizer.gridParams = cms.vstring()

    if isFastSim:
        process.susyNtuplizer.muonCollectionTags = cms.vstring("muons")
        process.susyNtuplizer.muonIdTags = cms.PSet()

        process.metFiltersSequence.remove(process.HBHENoiseFilterResultProducer)
        process.metFiltersSequence.remove(process.logErrorTooManyClusters)
        process.metFiltersSequence.remove(process.logErrorTooManyTripletsPairs)
        process.metFiltersSequence.remove(process.logErrorTooManySeeds)
        # met filters with run = False will automatically be default = False
        process.susyNtuplizer.metFilters.CSCBeamHalo.run = False
        process.susyNtuplizer.metFilters.HcalNoise.run = False
        process.susyNtuplizer.metFilters.LogErrorTooManyClusters.run = False
        process.susyNtuplizer.metFilters.LogErrorTooManyTripletsPairs.run = False
        process.susyNtuplizer.metFilters.LogErrorTooManySeeds.run = False

    if is52x or isFastSim:
        process.metFiltersSequence.remove(process.manystripclus53X)
        process.metFiltersSequence.remove(process.toomanystripclus53X)
        process.susyNtuplizer.metFilters.ManyStripClus53X.run = False
        process.susyNtuplizer.metFilters.TooManyStripClus53X.run = False

    if not runNoPUMVAMetSequence:
        process.susyNtuplizer.metCollectionTags.remove('pfNoPileUpMet')
        process.susyNtuplizer.metCollectionTags.remove('pfMVAMet')

    if dataset in ['52xPrompt', '52x23May2012', '53xPromptC', '53x13July2012', '53x06Aug2012', '53x24Aug2012', '53x11Dec2013']:
        process.susyNtuplizer.metFilters.HcalLaserEventList.default = True
        # The filter is actually applied at the sequence level; all events saved in the ntuples would have passed it

    if dataset in ['53xPromptD', '53x16Jan2013']:
        process.susyNtuplizer.metFilters.HcalLaserOccupancy.default = True

    if dataset == '53x22Jan2013':
        process.susyNtuplizer.metFilters.HcalLaserRECOUserStep.default = True
    else:
        process.susyNtuplizer.metFilters.HcalLaserRECOUserStep.run = False

    if dataset in ['53x13July2012', '53x24Aug2012']:
        process.susyNtuplizer.storeLumiInfo = cms.bool(False)

    if dataset == '53x22Jan2013':
        process.correctedMetSequence.remove(process.pfMEtSysShiftCorr)
        process.correctedMetSequence.remove(process.pfSysShiftCorrectedMet)
        process.correctedMetSequence.remove(process.pfType01SysShiftCorrectedMet)
        process.susyNtuplizer.metCollectionTags.remove('pfSysShiftCorrectedMet')
        process.susyNtuplizer.metCollectionTags.remove('pfType01SysShiftCorrectedMet')

    if is70x:
        process.standard_step.remove(process.pfCHSJetSequence)
    else:
        process.puRhoSequence.remove(process.kt6PFJets)

    return process
