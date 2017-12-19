import math
from sys import argv

import sys
sys.path.append('./')
sys.path.append('../python/')
 
from make_templates import TemplateMaker
from mix_data import MixData

from template_parameters import *
from template_parameters import params_test as cfg

if argv[1]=='grid':
    template = TemplateMaker(out_dir='../root/TEST/')
    template.set_CS_grid_mesh(grid_px_cos=200, grid_px_phi=200)
    template.make_grid_lab( params=cfg, coeff=cfg['params_CS']['coeff'], mass=cfg['params_CS']['mass'], ntoys=100000 )
elif argv[1]=='mix':
    mixdata = MixData(input_dir='../root/TEST/', output_dir='../data/TEST/')
    mixdata.add_shapes(cfg, symmetrise=True, make_templates=True)
    mixdata.mix_all_bins()
else:
    print "Unsupported option"



