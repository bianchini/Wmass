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

# no prior
prior_options_noprior = {}

# prior for y>=2.00
prior_options_y = copy.deepcopy(prior_options_base)
prior_options_y['select'] = ['y2.50', 'y3.00', 'y3.50']

# prior on A0 only
prior_options_A0 = copy.deepcopy(prior_options_base)
prior_options_A0['select'] = ['A0']

# prior on A1 only
prior_options_A1 = copy.deepcopy(prior_options_base)
prior_options_A1['select'] = ['A1']

# prior on A2 only
prior_options_A2 = copy.deepcopy(prior_options_base)
prior_options_A2['select'] = ['A2']

# prior on A3 only
prior_options_A3 = copy.deepcopy(prior_options_base)
prior_options_A3['select'] = ['A3']

# prior on A4 only
prior_options_A4 = copy.deepcopy(prior_options_base)
prior_options_A4['select'] = ['A4']

# prior on A0,A1,A2 only
prior_options_A0A1A2 = copy.deepcopy(prior_options_base)
prior_options_A0A1A2['select'] = ['A0','A1','A2']

# prior on A0,A2 only
prior_options_A0A2 = copy.deepcopy(prior_options_base)
prior_options_A0A2['select'] = ['A0','A2']

# prior on A3,A4 only
prior_options_A3A4 = copy.deepcopy(prior_options_base)
prior_options_A3A4['select'] = ['A3','A4']

templateFitter = TemplateFitter(DY='CC_FxFx'
                                ,charge='Wplus' 
                                ,var='WpreFSR'
                                ,job_name=
                                ,input_tag_fit=
                                ,input_tag_templ=
                                ,alternative_mc=''
                                ,mc_mass=80.419
                                ,num_events=3.0e+07
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
                                ,use_gradient=False
                                ,random_start=False
                                ,save_plots=[]
                                ,print_evals=True
                                ,run_on_crab=True
                                )

for i in range(ntoys):
    templateFitter.load_data( dataset= 
                              ,save_plots=[]
                              ,postfix='_'+str(i) )
    status = templateFitter.run(n_points=1000000, strategy=2, tolerance=0.1, run_minos=False, run_post_hesse=False)
    if status>0:
        continue
    templateFitter.update_results(print_results=False, 
                                  #save_plots=['norm', 'cov', 'coeff'], 
                                  #save_plots=[], 
                                  propagate_covariance=True)

templateFitter.close()

os.system("ls -lR")

