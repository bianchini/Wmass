#!/usr/bin/env python
import commands
import time
import re
import os
import string
from os import listdir
from os.path import isfile, join

import sys
sys.path.append('./')


################################################################################

def submit(path_to_script, script, num_events, ntoys, job_name, job_id):
    script_name = script+'_'+job_name+'_'+str(job_id)+'.sh'
    job_name    = script+'_'+job_name+'_'+str(job_id)
    out_name    = script+'_'+job_name+'_'+str(job_id)
    f = open(script_name,'w')
    f.write('#!/bin/bash\n\n')
    #f.write('mkdir /gpfs/ddn/cms/user/bianchi/\n')
    f.write('cd /home/users/bianchini/CMSSW_8_0_25/src/Wmass/test/\n')
    f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('\n')    
    f.write('python '+script+'.py '+str(num_events)+' '+str(ntoys)+' '+job_name+'\n')
    #f.write('mv /gpfs/ddn/cms/user/bianchi/'+out_name+'.root ./\n')    
    f.close()
    os.system('chmod +x '+script_name)

    submit_to_queue = 'bsub -q local -J '+job_name+' -o '+job_name+'.stdout -cwd `pwd` '+path_to_script+'/'+script_name
    print submit_to_queue
    os.system(submit_to_queue)
    time.sleep( 1.0 )

    print "@@@@@ END JOB @@@@@@@@@@@@@@@"


##########################################


path_to_script = '/home/users/bianchini/CMSSW_8_0_25/src/Wmass/test/'
for ijob, job in enumerate([
    ['run_unfolder', 1000000, 10, '1e6_pt_y_local' ], 
    ['run_unfolder', 1000000, 10, '1e6_pt_y_local' ], 
    ]):    
    submit(path_to_script, job[0], job[1], job[2], job[3], ijob)
