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
    weights = [0,1,2,3,4,6,8] + range(9,109)
    for q in ['Wplus', 'Wminus']:
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A0','A1', 'A2', 'A3', 'A4'], add_stat_uncert=True, postfix='_all')
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A3', 'A4'], add_stat_uncert=True, postfix='_A3A4')
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A0'], add_stat_uncert=True, postfix='_A0')
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A1'], add_stat_uncert=True, postfix='_A1')
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A2'], add_stat_uncert=True, postfix='_A2')
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A3'], add_stat_uncert=True, postfix='_A3')
        get_covariance(fname='./tree_histos.root', var='Wdress', q=q, weights=weights, coefficients=['A4'], add_stat_uncert=True, postfix='_A4')
