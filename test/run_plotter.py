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
        ['CC_FxFx_Wplus_random_parametric_y2p00_qt32_mass', 'W^{+}, toys, param-1D, |y| #leq 2.0, 30M, mass'],
        ['CC_FxFx_Wplus_asimov_parametric_y2p00_qt32_mass', 'W^{+}, asimov, param-1D, |y| #leq 2.0, 30M, mass'],
        ]

    for v in [#'bias', 
              #'rms',
              'biasANDrms'
              ]:
        profile_toys( files=files, alphas=['mass'],   var=v, postfix='_mass', save_pulls=True,  truth='val')
        profile_toys( files=files, alphas=['minuit'], var=v, postfix='',      save_pulls=False, truth='')
        #profile_toys( files=files, alphas=['norm'],   var=v, postfix='_norm', save_pulls=True,  truth='val')
        #profile_toys( files=files, alphas=['A0'],     var=v, postfix='_A0',   save_pulls=True,  truth='fit')
        #profile_toys( files=files, alphas=['A1'],     var=v, postfix='_A1',   save_pulls=True,  truth='fit')
        #profile_toys( files=files, alphas=['A2'],     var=v, postfix='_A2',   save_pulls=True,  truth='fit')
        #profile_toys( files=files, alphas=['A3'],     var=v, postfix='_A3',   save_pulls=True,  truth='fit')
        #profile_toys( files=files, alphas=['A4'],     var=v, postfix='_A4',   save_pulls=True,  truth='fit')
else:
    print "Option not available"
