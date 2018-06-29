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

from plot_utils import *

if argv[1]=='profile_toys':
    files = [
        ['CC_FxFx_Wplus_random_parametric_y1p40_qt32', 'W^{+}, toys, param-1D, |y| #leq 1.4, q_{T} #leq 32 GeV'],
        ['CC_FxFx_Wplus_random_parametric_y1p80_qt32', 'W^{+}, toys, param-1D, |y| #leq 1.8, q_{T} #leq 32 GeV'],
        ['CC_FxFx_Wplus_random_parametric_y2p50_qt32', 'W^{+}, toys, param-1D, |y| #leq 2.5, q_{T} #leq 32 GeV'],
        ['CC_FxFx_Wplus_random_parametric_y3p50_prior_options_y_qt32', 'W^{+}, toys, param-1D, |y| #leq 3.5, q_{T} #leq 32 GeV (prior 2.0#leq|y|#leq3.5)'],                                                               
        ]
    for v in ['bias', 
              'rms',
              'biasANDrms'
              ]:
        profile_toys( files=files, alphas=['norm'], var=v, postfix='_norm')
        profile_toys( files=files, alphas=['A0'],   var=v, postfix='_A0')
        profile_toys( files=files, alphas=['A1'],   var=v, postfix='_A1')
        profile_toys( files=files, alphas=['A2'],   var=v, postfix='_A2')
        profile_toys( files=files, alphas=['A3'],   var=v, postfix='_A3')
        profile_toys( files=files, alphas=['A4'],   var=v, postfix='_A4')
