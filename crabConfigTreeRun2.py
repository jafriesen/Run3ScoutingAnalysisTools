from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'scoutingTreeRun2_23May2023'

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'treeRun2.py'

config.Data.inputDataset = '/ScoutingCaloMuon/Run2018D-v1/RAW'
#config.Data.inputDataset = '/ScoutingPFMonitor/Run2022F-v1/RAW'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 100
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'scoutingTreeRun2_23May2023'

# These values only make sense for processing data
#    Select input data based on a lumi mask
config.Data.lumiMask = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'

# Where the output files will be transmitted to
config.Site.storageSite = 'T2_US_MIT'
