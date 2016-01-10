### this is an example for running on RECO
### the name must be changed crab.cfg for actual running

from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.General.requestName = 'PbPb_singleTrackTurnOns_v8'
config.General.workArea = 'PbPb_singleTrackTurnOns_v8'
config.General.transferOutputs = True
config.General.transferLogs = True
config.JobType.allowUndistributedCMSSW = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HLTandL1singletrack_crab_cfg.py'

config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/HIHardProbes/qwang-HIHardProbes_FullTrackSkim2015_v3-82d3c5ee469522058df894563cd74923/USER'
config.Data.splitting = 'FileBased'
config.Data.ignoreLocality = False
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Site.storageSite = 'T2_US_MIT'
