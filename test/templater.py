import math
from sys import argv

import sys
sys.path.append('./')
sys.path.append('../python/')
 
from make_templates import TemplateMaker
from mix_data import MixData

from template_parameters import *
from template_parameters import params_test as cfg

#cfg['params_W']['pt'] = np.linspace(0.0, 5.0, 11)
#cfg['params_W']['pt'] = np.linspace(25.5, 32.0, 14)
#cfg['params_W']['pt'] = np.linspace(25.5, 32.0, 14)

params_test['params_template']['pt'] = np.array([4.0, 8.0, 12.0, 16.0, 20.0, 26.0, 32.0])

if argv[1]=='grid':
    template = TemplateMaker(out_dir='../root/TEST/')
    template.set_CS_grid_mesh(grid_px_cos=200, grid_px_phi=200)
    template.make_grid_lab( params=cfg, coeff=cfg['params_CS']['coeff'], mass=cfg['params_CS']['mass'], ntoys=1000000 )
elif argv[1]=='mix':
    mixdata = MixData(input_dir='../root/TEST/', output_dir='../data/TEST/')
    mixdata.add_shapes(cfg, symmetrise=True, make_templates=True)
    mixdata.mix_all_bins()
else:
    print "Unsupported option"



