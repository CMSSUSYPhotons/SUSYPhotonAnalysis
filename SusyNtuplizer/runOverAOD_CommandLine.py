from SUSYPhotonAnalysis.SusyNtuplizer.runOverAOD import configure
import FWCore.ParameterSet.VarParsing as VarParsing
# Call by e.g. "cmsRun runOverAOD_CommandLine.py dataset=53x22Jan2013 inputFiles=file1.root file2.root"


## get command line arguments
opt = VarParsing.VarParsing("analysis")
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
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string, # string, int, or float
    "Please choose a dataset token. For examples see config file."
)
opt.register(
    "hltPaths",
    [],
    VarParsing.VarParsing.multiplicity.list,
    VarParsing.VarParsing.varType.string,
    "Please choose a hltPath token. For examples see config file."
)

opt.parseArguments()

process = configure( opt.dataset, opt.inputFiles, opt.hltPaths, opt.maxEvents )
