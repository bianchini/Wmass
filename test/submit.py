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

def submit(path_to_script, script, sample, first, last):
    script_name = script+'_'+sample+'_'+str(first)+'_'+str(last)+'.sh'
    job_name    = script+'_'+sample+'_'+str(first)+'_'+str(last)
    out_name    = script+'_'+sample+'_'+str(first)+'_'+str(last)
    f = open(script_name,'w')
    f.write('#!/bin/bash\n\n')
    f.write('mkdir /gpfs/ddn/cms/user/bianchi/\n')
    f.write('cd /home/users/bianchini/CMSSW_8_0_25/src/Wmass/test/\n')
    f.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('\n')    
    f.write('python '+script+'.py '+sample+' '+str(first)+' '+str(last)+'\n')
    f.write('mv /gpfs/ddn/cms/user/bianchi/'+out_name+'.root ./\n')    
    f.close()
    os.system('chmod +x '+script_name)

    submit_to_queue = 'bsub -q cms -J '+job_name+' -o '+job_name+'.stdout -cwd `pwd` '+path_to_script+'/'+script_name
    print submit_to_queue
    #os.system(submit_to_queue)
    time.sleep( 1.0 )

    print "@@@@@ END JOB @@@@@@@@@@@@@@@"


##########################################


path_to_script = '/home/users/bianchini/'
# [name, files_per_sample, njobs]
for sample in [
    ["hello", "TEST", 2,  2], 
    ]:
    for it in xrange(sample[3]):
        first = it*sample[2]+1
        last = (it+1)*sample[2]
        submit(path_to_script, sample[0], sample[1], first, last)
