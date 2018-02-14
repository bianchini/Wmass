from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TEST'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_script.py', 'FrameworkJobReport.xml', 'tree_utils.py', 'tree_producer.py', 'crabhelper.py']
config.JobType.sendPythonFolder	 = True
config.JobType.outputFiles = ['tree.root']

config.section_("Data")
config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
#config.Data.totalUnits = 10
config.Data.publication = False
config.Data.outputDatasetTag = 'TEST'
config.Data.outLFNDirBase = '/store/user/bianchi/'

config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"


datasets = [
    ['/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 'WJets_FxFx_ext0'],
    #['/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2/MINIAODSIM', 'WJets_FxFx_ext1']
    ]

from CRABAPI.RawCommand import crabCommand

for data in datasets:
    config.Data.inputDataset = data[0]
    config.General.requestName = 'TreeProducer_'+data[1]
    crabCommand('submit', config = config)    
