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
        #['CC_FxFx_Wplus_random_parametric_y2p00_qt32_mass_A1A3pol3', 'W^{+}, toys, param-1D, |y| #leq 2.0, 30M, mass A_{1,3} pol_{3}',   'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3pol3_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p00_qt32_mass_A1A3pol3', 'W^{+}, asimov, param-1D, |y| #leq 2.0, 30M, mass A_{1,3} pol_{3}', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3pol3_decorrelated'],
        ['CC_FxFx_Wplus_random_parametric_y2p00_qt32_mass_A1A3pol3_prior', 'W^{+}, toys, param-1D, |y| #leq 2.0, 30M, mass A_{1,3} pol_{3} prior',   'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3pol3_decorrelated'],
        ['CC_FxFx_Wplus_asimov_parametric_y2p00_qt32_mass_A1A3pol3_prior', 'W^{+}, asimov, param-1D, |y| #leq 2.0, 30M, mass A_{1,3} pol_{3} prior', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3pol3_decorrelated'],
        #['CC_FxFx_Wplus_random_parametric_y2p00_qt32_mass_A1A3A4pol3', 'W^{+}, toys, param-1D, |y| #leq 2.0, 30M, mass A_{1,3,4} pol_{3}',   'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3A4pol3_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p00_qt32_mass_A1A3A4pol3', 'W^{+}, asimov, param-1D, |y| #leq 2.0, 30M, mass A_{1,3,4} pol_{3}', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3A4pol3_decorrelated'],
        #['CC_FxFx_Wplus_random_parametric_y2p00_qt32_mass_A1A3A4pol3_prior', 'W^{+}, toys, param-1D, |y| #leq 2.0, 30M, mass A_{1,3,4} pol_{3} prior', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3A4pol3_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p00_qt32_mass_A1A3A4pol3_prior', 'W^{+}, asimov, param-1D, |y| #leq 2.0, 30M, mass A_{1,3,4} pol_{3} prior','CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_A1A3A4pol3_decorrelated']
        ]

    for v in [#'bias', 
              #'rms',
              'biasANDrms'
              ]:
        profile_toys( files=files, alphas=['mass'],   var=v, postfix='_mass', save_pulls=True,  truth='val')
        profile_toys( files=files, alphas=['minuit'], var=v, postfix='',      save_pulls=False, truth='')
        profile_toys( files=files, alphas=['norm'],   var=v, postfix='_norm', save_pulls=True,  truth='val')
        for coeff in ['A0', 'A1', 'A2', 'A3', 'A4']:
            profile_toys( files=files, alphas=[coeff], var=v, postfix='_'+coeff,   save_pulls=True,  truth='fit')

elif argv[1]=='derivatives':
    derivative_templates(charge='Wplus', var='WpreFSR', coeff_eval='val', tag='finer_y_qt32', bin=(0,0,1) )

elif argv[1]=='weighted_templates':
    weighted_templates(charge='Wplus', weights=range(109))

else:
    print "Option not available"
