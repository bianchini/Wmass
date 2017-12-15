import numpy as np
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from unfolder import Unfolder
from template_parameters import params_test

#params_test['params_template']['pt'] = np.linspace(0.0, 20.0, 6)
#params_test['params_template']['pt'] = np.linspace(0.0, 20.0, 2)
#params_test['params_template']['pt'] = np.linspace(0.0, 4.0, 2)
#params_test['params_template']['y']  = np.linspace(0.0, 3.6, 10)
#params_test['params_template']['y']  = np.linspace(0.0, 3.6, 2)
#params_test['params_template']['y']  = np.linspace(0.0, 0.4, 2)

# predictable
np.random.seed(0)

num_events = 1e+6
ntoys = 1
job_name = 'TEST'
fix = [ 'A0', 'mass', 'A4' ]
if len(argv)==4:
    num_events = int(argv[1])
    ntoys = int(argv[2])
    job_name = argv[3]

unfolder = Unfolder(input_dir=(os.environ['CMSSW_BASE']+'/src/Wmass/data/'), 
                    params=params_test, 
                    mass=80.000, 
                    num_events=num_events, 
                    fix=fix, 
                    interp_deg=1, 
                    job_name=job_name, 
                    verbose=True, 
                    prior_coeff=0.3, prior_xsec=0.3)

for itoy in range(ntoys):
    print "Process toy n.", itoy
    unfolder.toy_data(itoy)
    unfolder.reset_parameters()
    unfolder.run()
unfolder.save_result()
