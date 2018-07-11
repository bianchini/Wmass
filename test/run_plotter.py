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
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_release', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_random_parametric_y1p60', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 60 GeV'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_4sigmas', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, 4#sigma'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_3sigmas', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, 3#sigma']
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_noclosure', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, non-closure']
        #['CC_FxFx_Wplus_asimov_parametric_y1p60_qt32_noclosure', 'W^{+}, asimov, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, non-closure']
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A0', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{0}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A1', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{1}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A2', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{2}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A3', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{3}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A4', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{4}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A0A1A2', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{0,1,2}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A3A4', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{3,4}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_A0A2', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, prior on A_{0,2}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_floatA0A2', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, linear on A_{0,2}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_floatA0A2_pol3', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, pol_{3} on A_{i}'],
        #['CC_FxFx_Wplus_random_parametric_y1p60_qt32_floatA0A2_pol2', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, float A_{0,2}, pol_{2} on A_{i}'],
        ['CC_FxFx_Wplus_random_parametric_y1p60_qt32_strategy2', 'W^{+}, toys, param-1D, |y| #leq 1.6, q_{T} #leq 32 GeV, strategy 2'],
        #['CC_FxFx_Wplus_random_parametric_y1p40_qt32', 'W^{+}, toys, param-1D, |y| #leq 1.4, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_asimov_parametric_y1p40_qt32', 'W^{+}, asimov, param-1D, |y| #leq 1.4, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_random_parametric_y1p80_qt32', 'W^{+}, toys, param-1D, |y| #leq 1.8, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_asimov_parametric_y1p80_qt32', 'W^{+}, asimov, param-1D, |y| #leq 1.8, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_random_parametric_y2p50_qt32', 'W^{+}, toys, param-1D, |y| #leq 2.5, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p50_qt32', 'W^{+}, asimov, param-1D, |y| #leq 2.5, q_{T} #leq 32 GeV'],
        #['CC_FxFx_Wplus_random_parametric_y3p50_prior_options_y_qt32', 'W^{+}, toys, param-1D, |y| #leq 3.5, q_{T} #leq 32 GeV (prior 2.0#leq|y|#leq3.5)'],
        #['CC_FxFx_Wplus_asimov_parametric_y3p50_prior_options_y_qt32', 'W^{+}, asimov, param-1D, |y| #leq 3.5, q_{T} #leq 32 GeV (prior 2.0#leq|y|#leq3.5)'],
        ]
    for v in ['bias', 
              'rms',
              #'biasANDrms'
              ]:
        profile_toys( files=files, alphas=['norm'], var=v, postfix='_norm')
        profile_toys( files=files, alphas=['A0'],   var=v, postfix='_A0')
        profile_toys( files=files, alphas=['A1'],   var=v, postfix='_A1')
        profile_toys( files=files, alphas=['A2'],   var=v, postfix='_A2')
        profile_toys( files=files, alphas=['A3'],   var=v, postfix='_A3')
        profile_toys( files=files, alphas=['A4'],   var=v, postfix='_A4')

else:
    print "Option not available"
