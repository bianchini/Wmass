#!/usr/bin/env python
import numpy as np
import copy
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_fitter import TemplateFitter

#np.random.seed(0)

ntoys =

# template
prior_options_base = {'prior'  : 'sum', 
                      'select' : [], 
                      'inflate': 1e+03, 
                      'decorrelate' : []
                      }

# prior for y>=2.00
prior_options_y = copy.deepcopy(prior_options_base)
prior_options_y['select'] = ['y2.50', 'y3.00', 'y3.50']

# no prior
prior_options_noprior = {}

templateFitter = TemplateFitter(DY='CC_FxFx'
                                ,charge='Wplus' 
                                ,var='WpreFSR'
                                ,job_name=
                                ,input_tag_fit=
                                ,input_tag_templ=
                                ,alternative_mc=''
                                ,mc_mass=80.419
                                ,num_events=1.5e+06
                                ,verbose=True
                                ,fixed_parameters=
                                ,prior_options=
                                ,reduce_qt=-1 
                                ,reduce_y=
                                ,reduce_pt=0
                                ,fit_mode=
                                ,interpolation='quadratic'
                                ,use_prefit=False
                                ,add_nonclosure=True
                                ,save_plots=[]
                                ,print_evals=True
                                ,run_on_crab=True
                                )

for i in range(ntoys):
    templateFitter.load_data( dataset= 
                              ,save_plots=[]
                              ,postfix='_'+str(i) )
    status = templateFitter.run(n_points=100000, run_minos=False, run_post_hesse=False)
    if status>0:
        continue
    templateFitter.update_results(print_results=False, 
                                  #save_plots=['norm', 'cov'], 
                                  #save_plots=[], 
                                  propagate_covariance=True)

templateFitter.close()

os.system("ls -lR")

