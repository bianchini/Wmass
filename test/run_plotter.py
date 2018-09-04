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
        #['CC_FxFx_Wplus_random_parametric_y2p50_qt32_default', 'W^{+}, toys, param-1D, |y| #leq 2.5, 30M, mass', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p50_qt32_default', 'W^{+}, asimov, param-1D, |y| #leq 2.5, 30M, mass', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_random_parametric_y3p00_qt32_default', 'W^{+}, toys, param-1D, |y| #leq 3.0, 30M, mass', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y3p00_qt32_default', 'W^{+}, asimov, param-1D, |y| #leq 3.0, 30M, mass', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_random_parametric_y3p50_qt32_default', 'W^{+}, toys, param-1D, |y| #leq 3.5, 30M, mass', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y3p50_qt32_default', 'W^{+}, asimov, param-1D, |y| #leq 3.5, 30M, mass', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p00_qt32_default_pt55', 'W^{+}, asimov, param-1D, |y| #leq 2.0, 30M, mass, p_{T}<55 GeV', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        ['CC_FxFx_Wplus_random_parametric_y2p00_qt32_default_pt55', 'W^{+}, toys, param-1D, |y| #leq 2.0, 30M, mass, p_{T}<55 GeV', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        ]

    for v in [#'bias', 
              #'rms',
              'biasANDrms'
              ]:
        #profile_toys( files=files, alphas=['mass'],   var=v, postfix='_mass', save_pulls=True,  truth='val')
        profile_toys( files=files, alphas=['minuit'], var=v, postfix='',      save_pulls=False, truth='')
        #profile_toys( files=files, alphas=['norm'],   var=v, postfix='_norm', save_pulls=True,  truth='val')
        #for coeff in ['A0', 'A1', 'A2', 'A3', 'A4']:
        #    profile_toys( files=files, alphas=[coeff], var=v, postfix='_'+coeff,   save_pulls=True,  truth='fit')

elif argv[1]=='derivatives':
    bins_template_y = [ 0.,  0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,  1.6, 1.8,  2.]
    #bins_template_y = [ 0.,  0.2 ]
    bins_template_qt= [ 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 32.0 ]
    #bins_template_qt= [ 0.0, 4.0 ]
    for iqt,qt in enumerate(bins_template_qt):
        for iy,y in enumerate(bins_template_y):
            derivative_templates(charge='Wplus', var='WpreFSR', coeff_eval='val', tag='finer_y_qt32', bin=(iqt,iy,1), bins_qt=[], bins_y=[])

elif argv[1]=='weighted_templates':
    weighted_templates(charge='Wplus', weights=[4,8]#weights=range(109)
                       )

else:
    print "Option not available"
