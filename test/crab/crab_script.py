#!/usr/bin/env python
import os
import numpy as np

import sys
sys.path.append('./')

from unfolder import Unfolder
from template_parameters import params_test

params_test['params_template']['pt'] = np.array([0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 26.0, 32.0])
params_test['params_template']['y']  = np.array([0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4])

#params_test['params_template']['pt'] = np.linspace(0.0, 20.0, 2) #01
#params_test['params_template']['y']  = np.linspace(0.0, 3.6, 2)  #10

job_name =
num_events =
ntoys = 1
fix =
prior_coeff = 
prior_xsec = 

decorrelate = 
do_taylor_expansion = 
n_taylor = 
rebin =

unfolder = Unfolder(input_dir=(os.environ['CMSSW_BASE']+'/src/Wmass/data/'), 
                    params=params_test, 
                    rebin=rebin,
                    mass=80.000, 
                    num_events=num_events, 
                    fix=fix, 
                    interp_deg=1, 
                    n_points=1000000,
                    job_name=job_name, 
                    verbose=True, 
                    prior_coeff=prior_coeff, prior_xsec=prior_xsec,
                    strategy=0,
                    decorrelate=decorrelate,
                    decorrelate_full=False,
                    do_semianalytic=True,
                    do_taylor_expansion=do_taylor_expansion,
                    n_taylor=n_taylor,
                    run_minos=False
                    )

for itoy in range(ntoys):
    print "Process toy n.", itoy
    unfolder.toy_data(itoy)
    unfolder.reset_parameters()
    unfolder.run()
unfolder.save_result()

print "DONE"
os.system("ls -lR")

