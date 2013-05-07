###### Set the dataset type to use proper settings ######

dataset = '53x22Jan2013'

#########################################################

###### Set the input file name if running locally  ######

sourceNames = [
#    'root://xrootd.unl.edu//store/data/Run2012A/Photon/AOD/22Jan2013-v1/20000/FEF664CA-ED68-E211-A3FA-003048678B08.root'
#    'root://xrootd.unl.edu//store/mc/Summer12_DR53X/GVJets_Incl_8TeV-madgraph/AODSIM/PU_S10_START53_V7C-v1/00001/FEE4A970-1333-E211-A759-00261894398B.root'    
]

#########################################################

### HLT results filter (leave empty if not filtering) ###

hltPaths = [
    'HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*',
    'HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*',
    'HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*',
    'HLT_Photon36_R9Id85_Photon22_R9Id85_v*',
    'HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*'
]

#########################################################

# Maximum number of events to process (ignored in CRAB) #

maxEvents = -1

#########################################################

collisionDatasets = [
    '52xPrompt',     # Run2012[AB]-PromptReco
    '52x23May2012',  # Run2012A-23May2012
    '53xPrompt',     # Run2012C-PromptReco-v2, Run2012D-PromptReco
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
    '53xFastSim'
]

isRealData = dataset in collisionDatasets
isMC = dataset in mcDatasets
isFastSim = isMC and 'FastSim' in dataset
is53x = '53x' in dataset
is52x = '52x' in dataset

if not isRealData and not isMC:
    raise RuntimeError("Dataset " + dataset + " not defined")

import FWCore.ParameterSet.Config as cms

##########################
### Initialize process ###
##########################
process = cms.Process("RA3")
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(maxEvents))

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
    'newSecondaryVertoexTagInfos',
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
elif dataset == '53xPrompt':
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
    process.GlobalTag.globaltag = 'FT_53_V21_AN3::All'
elif isMC and is52x:
    process.GlobalTag.globaltag = 'START52_V16::All'
elif isMC and is53x:
    process.GlobalTag.globaltag = 'START53_V21::All'

#####################
### SusyNtuplizer ###
#####################
process.load("SusyAnalysis.SusyNtuplizer.susyNtuplizer_cfi")
process.susyNtuplizer.debugLevel = 0
process.susyNtuplizer.isFastSim = isFastSim

if dataset == '53x13July2012' or dataset == '53x24Aug2012':
    process.susyNtuplizer.storeLumiInfo = cms.bool(False)

if isRealData:
    process.susyNtuplizer.metCollectionTags.remove('genMetTrue')
    process.susyNtuplizer.jetFlavourMatchingTags = cms.PSet()
    process.susyNtuplizer.gridParams = cms.vstring()
elif isFastSim:
    process.susyNtuplizer.muonCollectionTags = cms.vstring("muons")
    process.susyNtuplizer.muonIDCollectionTags = cms.PSet()
    process.susyNtuplizer.metFilters.remove('CSCBeamHalo')
    process.susyNtuplizer.metFilters.remove('HcalNoise')

if is52x:
    process.susyNtuplizer.metFilters.remove('ManyStripClus53X')
    process.susyNtuplizer.metFilters.remove('TooManyStripClus53X')

if dataset != '53x22Jan2013':
    process.susyNtuplizer.metFilters.remove('HcalLaserRECOUserStep')

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
# Runs PU identificaton over particleFlow, then calculates isodeposits around
# the given particles using the NoPU collection.
# The end product is equivalent to running setupPFIso in
# CommonTools.ParticleFlow.Tools.pfIsolation, just with fewer modules to run.

process.photonPFIsoDepositCharged = process.phPFIsoDepositCharged.clone(src = 'photons')
process.photonPFIsoDepositNeutral = process.phPFIsoDepositNeutral.clone(src = 'photons')
process.photonPFIsoDepositGamma = process.phPFIsoDepositGamma.clone(src = 'photons')
process.gsfElectronPFIsoDepositCharged = process.elPFIsoDepositCharged.clone(src = 'gsfElectrons')
process.gsfElectronPFIsoDepositNeutral = process.elPFIsoDepositNeutral.clone(src = 'gsfElectrons')
process.gsfElectronPFIsoDepositGamma = process.elPFIsoDepositGamma.clone(src = 'gsfElectrons')
process.photonPFIsoValueCharged03 = process.phPFIsoValueCharged03PFId.clone()
process.photonPFIsoValueCharged03.deposits[0].src = cms.InputTag('photonPFIsoDepositCharged')
process.photonPFIsoValueNeutral03 = process.phPFIsoValueNeutral03PFId.clone()
process.photonPFIsoValueNeutral03.deposits[0].src = cms.InputTag('photonPFIsoDepositNeutral')
process.photonPFIsoValueGamma03 = process.phPFIsoValueGamma03PFId.clone()
process.photonPFIsoValueGamma03.deposits[0].src = cms.InputTag('photonPFIsoDepositGamma')
process.gsfElectronPFIsoValueCharged03 = process.elPFIsoValueCharged03PFId.clone()
process.gsfElectronPFIsoValueCharged03.deposits[0].src = cms.InputTag('gsfElectronPFIsoDepositCharged')
process.gsfElectronPFIsoValueNeutral03 = process.elPFIsoValueNeutral03PFId.clone()
process.gsfElectronPFIsoValueNeutral03.deposits[0].src = cms.InputTag('gsfElectronPFIsoDepositNeutral')
process.gsfElectronPFIsoValueGamma03 = process.elPFIsoValueGamma03PFId.clone()
process.gsfElectronPFIsoValueGamma03.deposits[0].src = cms.InputTag('gsfElectronPFIsoDepositGamma')

process.photonIsoDepositSequence = cms.Sequence(
    process.photonPFIsoDepositCharged +
    process.photonPFIsoDepositNeutral +
    process.photonPFIsoDepositGamma +
    process.photonPFIsoValueCharged03 +
    process.photonPFIsoValueNeutral03 +
    process.photonPFIsoValueGamma03
)
process.gsfElectronIsoDepositSequence = cms.Sequence(
    process.gsfElectronPFIsoDepositCharged +
    process.gsfElectronPFIsoDepositNeutral +
    process.gsfElectronPFIsoDepositGamma +
    process.gsfElectronPFIsoValueCharged03 +
    process.gsfElectronPFIsoValueNeutral03 +
    process.gsfElectronPFIsoValueGamma03
)

process.pfIsolationSequence = cms.Sequence(
    process.photonIsoDepositSequence +
    process.gsfElectronIsoDepositSequence
)

###############################################
### PFchs (charged hadron subtraction) jets ###
###############################################
# pfJets comes from PFBRECO.

process.ak5PFchsJets = process.pfJets.clone()

process.pfCHSJetSequence = cms.Sequence(
    process.ak5PFchsJets
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
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")

if isRealData:
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3Residual")
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")
else:
    process.caloJetMETcorr.jetCorrLabel = cms.string("ak5CaloL2L3")
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")

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

process.pfType0CorrectedMet = process.pfType1CorrectedMet.clone(
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfType0MetCorrection')
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
process.pfType01p2SysShiftCorrectedMet = process.pfType1p2CorrectedMet.clone(
    srcType1Corrections = cms.VInputTag(
        cms.InputTag('pfType0MetCorrection'),
        cms.InputTag('pfJetMETcorr', 'type1'),
        cms.InputTag('pfMEtSysShiftCorr')
    )
)

process.correctedMetSequence = cms.Sequence(
    process.produceCaloMETCorrections +
    process.producePFMETCorrections +
    process.pfMEtSysShiftCorrSequence +
    process.pfCandidateToVertexAssoc +
    process.pfType0MetCorrection +
    process.pfType0CorrectedMet +
    process.pfSysShiftCorrectedMet +
    process.pfType01SysShiftCorrectedMet +
    process.pfType01p2SysShiftCorrectedMet
)

#############################
### MVA-based electron ID ###
#############################
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(
    process.mvaTrigV0 *
    process.mvaNonTrigV0
)

#################
### PU jet ID ###
#################
from RecoJets.JetProducers.PileupJetIDParams_cfi import full_5x, full_5x_chs, cutbased
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetIdProducer, pileupJetIdProducerChs

process.recoPuJetId = pileupJetIdProducer.clone(
    jets = cms.InputTag("ak5PFJets"),
    algos = cms.VPSet(full_5x, cutbased),
    applyJec = cms.bool(True),
    inputIsCorrected = cms.bool(False),
    residualsTxt = cms.FileInPath("RecoJets/JetProducers/data/mva_JetID_v1.weights.xml") # not used, but has to point to an existing file
)

process.recoPuJetIdChs = pileupJetIdProducerChs.clone(
    jets = cms.InputTag("ak5PFchsJets"),
    algos = cms.VPSet(full_5x_chs, cutbased),
    applyJec = cms.bool(True),
    inputIsCorrected = cms.bool(False),
    residualsTxt = cms.FileInPath("RecoJets/JetProducers/data/mva_JetID_v1.weights.xml") # not used, but has to point to an existing file
)

process.recoPuJetIdSequence = cms.Sequence(
    process.recoPuJetId +
    process.recoPuJetIdChs
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
        jets = cms.InputTag("ak5PFchsJets"),
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
    jets = cms.InputTag("ak5PFchsJets"),
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
process.newSoftElectronTagInfos = process.softElectronTagInfos.clone(
    jets = "ak5PFJets"
    )
process.newSoftElectronBJetTags = process.softElectronBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("newSoftElectronTagInfos") )
    )
process.chsSoftElectronTagInfos = process.softElectronTagInfos.clone(
    jets = "ak5PFchsJets"
    )
process.chsSoftElectronBJetTags = process.softElectronBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("chsSoftElectronTagInfos") )
    )
# soft muon b-tag
process.newSoftMuonTagInfos = process.softMuonTagInfos.clone(
    jets = "ak5PFJets"
    )
process.newSoftMuonBJetTags = process.softMuonBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
    )
process.newSoftMuonByIP3dBJetTags = process.softMuonByIP3dBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
    )
process.newSoftMuonByPtBJetTags = process.softMuonByPtBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("newSoftMuonTagInfos") )
    )
process.chsSoftMuonTagInfos = process.softMuonTagInfos.clone(
    jets = "ak5PFchsJets"
    )
process.chsSoftMuonBJetTags = process.softMuonBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("chsSoftMuonTagInfos") )
    )
process.chsSoftMuonByIP3dBJetTags = process.softMuonByIP3dBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("chsSoftMuonTagInfos") )
    )
process.chsSoftMuonByPtBJetTags = process.softMuonByPtBJetTags.clone(
    tagInfos = cms.VInputTag( cms.InputTag("chsSoftMuonTagInfos") )
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
    process.softElectronCands *
    process.newSoftElectronTagInfos *
    process.newSoftElectronBJetTags
)
process.chsJetBtaggingEle = cms.Sequence(
    process.softElectronCands *
    process.chsSoftElectronTagInfos *
    process.chsSoftElectronBJetTags
)

process.newJetBtaggingMu = cms.Sequence(
    process.newSoftMuonTagInfos * (
        process.newSoftMuonBJetTags +
        process.newSoftMuonByIP3dBJetTags +
        process.newSoftMuonByPtBJetTags
    )
)
process.chsJetBtaggingMu = cms.Sequence(
    process.chsSoftMuonTagInfos * (
        process.chsSoftMuonBJetTags +
        process.chsSoftMuonByIP3dBJetTags +
        process.chsSoftMuonByPtBJetTags
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
from QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff import goodOfflinePrimaryVerticesQG, QGTagger

process.goodOfflinePrimaryVerticesQG = goodOfflinePrimaryVerticesQG

process.QGTaggerAK5 = QGTagger.clone(
    srcJets = cms.InputTag('ak5PFJets'),
    srcRho = cms.InputTag('kt6PFJets', 'rho'),
    srcRhoIso = cms.InputTag('kt6PFJetsRho25', 'rho')
)
process.QGTaggerAK5chs = QGTagger.clone(
    srcJets = cms.InputTag('ak5PFchsJets'),
    useCHS = cms.untracked.bool(True),
    srcRho = cms.InputTag('kt6PFJets', 'rho'),    
    srcRhoIso = cms.InputTag('kt6PFJetsRho25', 'rho')
)

if isRealData:
    process.QGTaggerAK5.jec = cms.untracked.string('ak5PFL1FastL2L3Residual')
    process.QGTaggerAK5chs.jec = cms.untracked.string('ak5PFchsL1FastL2L3Residual')
else:
    process.QGTaggerAK5.jec = cms.untracked.string('ak5PFL1FastL2L3')
    process.QGTaggerAK5chs.jec = cms.untracked.string('ak5PFchsL1FastL2L3')

process.QGTaggingSequence = cms.Sequence(
    process.goodOfflinePrimaryVerticesQG +
    process.QGTaggerAK5 +
    process.QGTaggerAK5chs
)

###################
### MET filters ###
###################
# HBHENoiseFilterResultProducer
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')

# HCAL laser events filter
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

if isFastSim:
    process.metFiltersSequence.remove(process.HBHENoiseFilterResultProducer)

if is52x:
    process.metFiltersSequence.remove(process.manystripclus53X)
    process.metFiltersSequence.remove(process.toomanystripclus53X)

######################
### NoPU & MVA MET ###
######################

####### SET THIS TO TRUE TO RUN THE SEQUENCE (CHECK THE README FILE FOR INSTRUCTIONS) ########
runNoPUMVAMetSequence = True
##############################################################################################

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
    process.susyNtuplizer.metCollectionTags.remove('pfNoPileUpMet')
    process.susyNtuplizer.metCollectionTags.remove('pfMVAMet')
    process.noPUMVAMetSequence = cms.Sequence()

#########################
### HLT result filter ###
#########################
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.HLTPaths = hltPaths

#####################
### Finalize path ###
#####################
process.p = cms.Path(
    process.hltHighLevel +
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
    process.metFiltersSequence +
    process.noPUMVAMetSequence +
    process.susyNtuplizer
)
