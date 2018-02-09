import os

samples = [
    '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer15GS-MCRUN2_71_V1_ext2-v1/GEN-SIM',
    '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
    '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2/MINIAODSIM',
    '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer15GS-MCRUN2_71_V1-v1/GEN-SIM',
    '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM',
    '/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM',
    '/WplusToMuNu_NNPDF30_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISummer15wmLHEGS-MCRUN2_71_V1-v2/GEN-SIM',
    '/WplusToMuNu_NNPDF30_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM',
    '/WplusToMuNu_NNPDF30_TuneEE5C_13TeV-powheg-herwigpp/RunIISummer15wmLHEGS-MCRUN2_71_V1-v2/GEN-SIM',
    '/WplusToMuNu_NNPDF30_TuneEE5C_13TeV-powheg-herwigpp/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    ]

for sample in samples:
    print "**************************************"
    tags = sample.split('/')
    print tags[1], "==>", tags[-1]
    os.system('das_client.py --query="dataset='+sample+' | grep dataset.nevents, dataset.nfiles"')
    os.system('das_client.py --query="release dataset='+sample+'"')

