import math

import sys
sys.path.append('./')
sys.path.append('../python/')
 
from make_templates import TemplateMaker
from mix_data import MixData

from template_parameters import *
from template_parameters import params_test as cfg

template = TemplateMaker(out_dir='../data/')
#template.set_CS_grid_mesh(grid_px_cos=200, grid_px_phi=200)
#template.make_grid_lab( params = cfg, coeff=cfg['params_CS']['coeff'], mass=cfg['params_CS']['mass'], ntoys=200000 )


mixdata = MixData(input_dir='../data/', output_dir='../data/')
mixdata.add_shapes(cfg, symmetrise=True, make_templates=True)
mixdata.mix_all_bins()



