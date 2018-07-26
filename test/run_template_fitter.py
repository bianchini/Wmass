import numpy as np
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_fitter import TemplateFitter

#np.random.seed(0)

ntoys = 50

prior_options = {'prior'  : 'sum', 
                 'select' : [], 
                 'inflate': 1e+03, 
                 'decorrelate' : []}
prior_options = {}

templateFitter = TemplateFitter(DY='CC_FxFx', charge='Wplus', var='WpreFSR', job_name='TEST_p3', 
                                input_tag_fit='all_A0-4_forced_v4_finer_y_qt20_decorrelated',
                                input_tag_templ='_finer_y_qt20',
                                alternative_mc='',
                                mc_mass=80.419, 
                                num_events=1.5e+06,
                                verbose=False, 
                                fixed_parameters=['pol', 'A', 'mass'],  
                                prior_options=prior_options,
                                reduce_qt=-1, 
                                reduce_y=-13,
                                reduce_pt=0,
                                fit_mode='binned',
                                do_analytic_fit=False,
                                interpolation='quadratic',
                                use_prefit=False,
                                add_nonclosure=True,
                                save_plots=[],
                                print_evals=True
                                )

for i in range(ntoys):
    templateFitter.load_data( dataset='random', save_plots=[], postfix='_'+str(i) )
    status = templateFitter.run(n_points=100000, strategy=2, tolerance=0.1, run_minos=False, run_post_hesse=False)
    if status>0:
        continue

    #templateFitter.scan(par_name='y0.00_y0.20_qt16.0_qt20.0_norm', step=0.2*9.50165e+03, par_low=1.38e+04-2*9.5e+03, par_high=1.38e+04+2*9.5e+03)

    templateFitter.update_results(print_results=True, 
                                  #save_plots=['cov', 'norm', 'coeff', 'polynom'], 
                                  #save_plots=['cov', 'norm'], 
                                  save_plots=[], 
                                  propagate_covariance=False)

templateFitter.close()
