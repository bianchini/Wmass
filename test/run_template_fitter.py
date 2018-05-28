import numpy as np
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_fitter import TemplateFitter

np.random.seed(0)

ntoys = 1

prior_options = {'prior'  : 'sum', 
                 'select' : ['A2', 'A3'], 
                 'inflate': 1e+03, 
                 'uncorrelate' : False}
prior_options = {}

templateFitter = TemplateFitter(DY='CC_FxFx', charge='Wplus', var='WpreFSR', job_name='TEST', mc_mass=80.419, 
                                num_events=1.5e+06,
                                verbose=False, 
                                fixed_parameters=['pol', 'A', 'mass'], 
                                prior_options=prior_options,
                                reduce_qt=-1, 
                                reduce_y=-8,
                                reduce_pt=0,
                                fit_mode='parametric',
                                interpolation='quadratic',
                                use_prefit=False,
                                add_nonclosure=True,
                                save_plots=[],
                                print_evals=True
                                )

for i in range(ntoys):
    templateFitter.load_data( dataset='asimov', save_plots=[], postfix='_'+str(i) )
    status = templateFitter.run(n_points=100000, run_minos=False, run_post_hesse=False)
    if status>0:
        continue
    templateFitter.update_results(print_results=True, save_plots=['coeff'], propagate_covariance=True)

templateFitter.close()
