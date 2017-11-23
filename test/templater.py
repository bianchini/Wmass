import math

import sys
sys.path.append('./')
sys.path.append('../python/')
 
from make_templates import TemplateMaker
from mix_data import MixData

from template_parameters import *

template = TemplateMaker(out_dir='../data/')
#template.set_CS_grid_mesh(grid_px_cos=200, grid_px_phi=200)
#template.make_grid_lab( params = params_test, coeff=params_test['params_CS']['coeff'], mass=params_test['params_CS']['mass'], ntoys=100000 )

mixdata = MixData(input_dir='../data/', output_dir='../data/')
mixdata.mix(mass=80.000, pdf=pdf_test, params=params_test)



