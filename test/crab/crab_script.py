#!/usr/bin/env python
import os
import numpy as np

import sys
sys.path.append('./')

from unfolder import Unfolder
from template_parameters import params_test

#params_test['params_template']['pt'] = np.linspace(0.0, 20.0, 2) #01
#params_test['params_template']['y']  = np.linspace(0.0, 3.6, 2)  #10

job_name =
num_events =
ntoys = 5
fix =
prior_coeff = 
prior_xsec = 

unfolder = Unfolder(input_dir=(os.environ['CMSSW_BASE']+'/src/Wmass/data/'), 
                    params=params_test, 
                    mass=80.000, 
                    num_events=num_events, 
                    fix=fix, 
                    interp_deg=1, 
                    job_name=job_name, 
                    verbose=False, 
                    prior_coeff=prior_coeff, 
                    prior_xsec=prior_xsec)

for itoy in range(ntoys):
    print "Process toy n.", itoy
    unfolder.toy_data(itoy)
    unfolder.reset_parameters()
    unfolder.run()
unfolder.save_result()

print "DONE"
os.system("ls -lR")

