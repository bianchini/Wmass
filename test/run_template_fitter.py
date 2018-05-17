import numpy as np
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_fitter import TemplateFitter

np.random.seed(0)

ntoys = 1

templateFitter = TemplateFitter(DY='CC_FxFx', charge='Wplus', var='WpreFSR', job_name='TEST', mc_mass=80.419, 
                                num_events=1.5e+06,
                                verbose=False, 
                                fixed_parameters=['pol', 'mass', 'A'], 
                                use_prior=False, 
                                reduce_qt=-1, 
                                reduce_y=-1, 
                                do_parametric=True,
                                use_prefit=False,
                                add_nonclosure=False
                                )

for i in range(ntoys):
    templateFitter.load_data( dataset='asymov', postfix='_'+str(i) )
    templateFitter.run(n_points=100000, run_minos=False, run_post_hesse=False)
    templateFitter.update_results()

templateFitter.close()
