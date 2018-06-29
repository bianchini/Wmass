from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TEST'
config.General.transferLogs=True

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.inputFiles = ['crab_script.py', 'FrameworkJobReport.xml', 'tree_utils.py', 'fit_utils.py', 'template_fitter.py']
config.JobType.sendPythonFolder	 = True
#config.JobType.maxJobRuntimeMin = 1315
config.JobType.outputFiles = ['result'+'.root']

config.section_("Data")
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 500
#config.Data.totalUnits = 1
config.Data.publication = False
config.Data.outputDatasetTag = 'TEST'
config.Data.outLFNDirBase = '/store/user/bianchi/'

config.section_("Site")
config.Site.storageSite = "T2_IT_Pisa"

from sys import argv
import copy 

job_base = {'job_name'         : 'TEST', 
            'ntoys'            : 1,
            'dataset'          : 'random',
            'input_tag_fit'    : 'all_A0-4_forced_v4_finer_y_qt32_decorrelated', 
            'input_tag_templ'  : '_finer_y_qt32',
            'fixed_parameters' : ['pol', 'A', 'mass'], 
            'fit_mode'         : 'parametric',
            'reduce_y'         : -1,
            'prior_options'    : 'prior_options_noprior'
            }

bins_template_y = [ 0.,0.2,  0.4, 0.6, 0.8, 1.0, 1.2, 1.4,  1.6, 1.8,  2. ,  2.5,  3. , 3.5]

jobs = []
for dataset in ['random']:
    for fit_mode in ['parametric']:
        for reduce_y in [-6, 
                          #-4, -3
                          ]:
            job_new = copy.deepcopy(job_base)
            job_new['job_name'] = dataset+'_'+fit_mode+'_'+('y{:03.2f}'.format(bins_template_y[reduce_y])).replace('.', 'p')+'_'+'qt32'+'_release'
            job_new['fit_mode'] = fit_mode
            job_new['reduce_y'] = reduce_y
            job_new['dataset'] = dataset
            jobs.append(job_new)

    job_new = copy.deepcopy(job_base)
    job_new['job_name'] = dataset+'_'+'parametric'+'_'+('y{:03.2f}'.format(bins_template_y[-1])).replace('.', 'p')+'_'+'prior_options_y'+'_qt32'
    job_new['fit_mode'] = 'parametric'
    job_new['reduce_y'] = -1
    job_new['prior_options'] = 'prior_options_y'
    #jobs.append(job_new)


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    if argv[1]=='create':

        import os
        os.system('cp '+os.environ['CMSSW_BASE']+'/src/Wmass/python/template_fitter.py .')
        os.system('cp '+os.environ['CMSSW_BASE']+'/src/Wmass/python/tree_utils.py .')
        os.system('cp '+os.environ['CMSSW_BASE']+'/src/Wmass/python/fit_utils.py .')

        for job in jobs:
            fout = open('crab_script_'+job['job_name']+'.py', 'w') 
            fin = open('crab_script.py', 'r') 
            lines = fin.readlines()
            for line in lines:
                line = line.rstrip()

                if 'job_name' in line:
                    line += '"'+job['job_name']+'"' 
                elif 'ntoys =' in line:
                    line += str(job['ntoys'])
                elif 'dataset=' in line:
                    line += '"'+job['dataset']+'"'
                elif 'input_tag_fit' in line:
                    line += '"'+job['input_tag_fit']+'"'
                elif 'input_tag_templ' in line:
                    line += '"'+job['input_tag_templ']+'"'
                elif 'fit_mode' in line:
                    line += '"'+job['fit_mode']+'"'
                elif 'prior_options=' in line:
                    line += job['prior_options']
                elif 'fixed_parameters' in line:
                    line += '['
                    for fixes in job['fixed_parameters']:
                        line += ('"'+fixes+'"')
                        line += ', '
                    line += ']'
                elif 'reduce_y' in line:
                    line += str(job['reduce_y'])
                elif "#save_plots=[]" in line and job['dataset']=='random':
                    line = line.replace('#', '')
                elif "#save_plots=['norm', 'cov']" in line and job['dataset']=='asimov':
                    line = line.replace('#', '')

                line += '\n'
                fout.write(line)
            fout.close()
            fin.close()

            # scritpExe
            fout = open('crab_script_'+job['job_name']+'.sh', 'w') 
            fin = open('crab_script.sh', 'r') 
            lines = fin.readlines()
            for line in lines:
                line = line.rstrip()
                if 'python crab_script_' in line:
                    line += (job['job_name']+'.py $1')
                fout.write(line+'\n')
            fout.close()
            fin.close()

    if argv[1]=='submit':
        for job in jobs:
            config.General.requestName = 'fitter_'+job['job_name']
            config.JobType.scriptExe = 'crab_script_'+job['job_name']+'.sh'
            config.JobType.inputFiles = ['crab_script_'+job['job_name']+'.py', 'FrameworkJobReport.xml', 'tree_utils.py', 'fit_utils.py', 'template_fitter.py']
            config.JobType.outputFiles = ['result_CC_FxFx_Wplus_'+job['job_name']+'.root']
            if job['dataset']=='asimov':
                config.Data.totalUnits = 1     
                config.JobType.outputFiles.extend([ 'covariance_fit_'+job['job_name']+'.png', 'covariance_fit_'+job['job_name']+'.C',
                                                    'norm_resolution_'+job['job_name']+'.png', 'norm_resolution_'+job['job_name']+'.C'])
                
            crabCommand('submit', config = config)

    if argv[1]=='killall':
        import os
        for job in jobs:
            os.system('crab kill -d '+'crab_'+'fitter_'+job['job_name'])
            #os.system('rm -r '+'crab_'+'fitter_'+job['job_name'])
        #os.system('rm crab_script_*.py crab_script_*.sh')

    if argv[1]=='getoutput':
        import os
        for job in jobs:
            os.system('crab getoutput -d '+'crab_'+'fitter_'+job['job_name'])
