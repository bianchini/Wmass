import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv

import os.path
from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

from fit_utils import draw_y_slice, draw_qt_slice, get_covariance

if argv[1]=='fit':
    for coeff in ['A0','A1','A2','A3','A4','A5', 'A6', 'A7']:
    #for coeff in ['A0']:
        draw_y_slice(fname='./tree_histos.root', var='Wdress', coeff=coeff, weight_name=1, do_fit=True)
        draw_qt_slice(fname='./tree.root', var='Wdress', coeff=coeff, weight_name=0, do_fit=False)

elif argv[1]=='cov':
    weights = {}
    weights['scale'] = ([0] + [1,2,3,4,6,8])
    weights['pdf'] = range(9,109)
    for q in ['Wplus', 'Wminus']:
        for coeff in ['A0','A1','A2','A3','A4']:
            get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=[coeff], add_stat_uncert=True, postfix='_'+coeff)
