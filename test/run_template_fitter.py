import numpy as np
import os
import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_fitter import TemplateFitter

templateFitter = TemplateFitter(DY='CC', charge='Wplus', var='WpreFSR', job_name='TEST', mc_mass=80.419, verbose=True, fixed_parameters=[])
templateFitter.run(n_points=50000, run_minos=False)
