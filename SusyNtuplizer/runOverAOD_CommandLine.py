import FWCore.ParameterSet.Config as cms

from SUSYPhotonAnalysis.SusyNtuplizer.runOverAOD import configure
from FWCore.ParameterSet.VarParsing import VarParsing
# Call by e.g. "cmsRun runOverAOD_CommandLine.py dataset=53x22Jan2013 inputFiles=file1.root file2.root"


## get command line arguments
opt = VarParsing("analysis")
# 'analysis' contains some default options:
# maxEvents: Number of events to process (singleton integer)
# inputFiles: List of files to process as input (list string)
# secondaryFiles: List of secondary files (if needed; list string)
# outputFile: Name of output file (singleton string)
# secondaryOutput: Name of secondary output (if needed; singleton string)

# register an additional option
opt.register(
    "dataset",
    "", # default value
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string, # string, int, or float
    "Please choose a dataset token. For examples see config file."
)
opt.register(
    "hltPaths",
    [],
    VarParsing.multiplicity.list,
    VarParsing.varType.string,
    "Please choose a hltPath token. For examples see config file."
)
opt.register(
    "outputName",
    "susyEvents.root",
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Please choose an output file name."
)

opt.parseArguments()

process = configure(dataset = opt.dataset, sourceNames = opt.inputFiles, hltPaths = opt.hltPaths, maxEvents = opt.maxEvents, outputName = opt.outputName)
