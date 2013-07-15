from runOverAOD import configure

sourceNames = [
    'root://xrootd.unl.edu//store/data/Run2012A/Photon/AOD/22Jan2013-v1/20000/FEF664CA-ED68-E211-A3FA-003048678B08.root'
]

hltPaths = [
    'HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v*',
    'HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v*',
    'HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v*',
    'HLT_Photon36_R9Id85_Photon22_R9Id85_v*',
    'HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v*'
]

process = configure( '53x22Jan2013', sourceNames, hltPaths, runNoPUMVAMetSequence=False)
