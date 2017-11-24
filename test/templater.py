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
#template.make_grid_lab( params = cfg, coeff=cfg['params_CS']['coeff'], mass=cfg['params_CS']['mass'], ntoys=10000 )


mixdata = MixData(input_dir='../data/', output_dir='../data/')
if True:
    for im,m in enumerate(cfg['params_template']['mass']):
        for ic0,c0 in enumerate(cfg['params_template']['A0']):
            for ic1,c1 in enumerate(cfg['params_template']['A1']):
                for ic2,c2 in enumerate(cfg['params_template']['A2']):                    
                    for ic3,c3 in enumerate(cfg['params_template']['A3']):
                        for ic4,c4 in enumerate(cfg['params_template']['A4']):    
                            for ipt in range(len(cfg['params_template']['pt'])-1):
                                pt_bin=[ cfg['params_template']['pt'][ipt], cfg['params_template']['pt'][ipt+1] ]
                                for iy in range(len(cfg['params_template']['y'])-1):
                                    y_bin=[ cfg['params_template']['y'][iy], cfg['params_template']['y'][iy+1] ]
                                    print "Call MixData for", pt_bin, y_bin
                                    mixdata.mix(mass=m, pdf=pdf_test, params=cfg, pt_bin=pt_bin, y_bin=y_bin, coeff=[c0,c1,c2,c3,c4])



