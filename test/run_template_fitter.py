import numpy as np
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_fitter import TemplateFitter

#np.random.seed(0)

templateFitter = TemplateFitter(DY='CC_FxFx', charge='Wplus', var='WpreFSR', job_name='TEST', mc_mass=80.419, 
                                verbose=True, 
                                fixed_parameters=['pol', 'mass', 'A'], 
                                use_prior=False, 
                                reduce_qt=-9, 
                                reduce_y=-8, 
                                debug=False, 
                                do_parametric=False,
                                dataset='rnd'
                                )

templateFitter.run(n_points=50000, run_minos=False, run_post_hesse=True)
