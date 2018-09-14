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
        #['CC_FxFx_Wplus_asimov-scaled-full_parametric_y2p00_qt32_noprior_weightsall',  'W^{+}, |y| #leq 2.0, q_{T}<32 GeV, 30M events', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov-scaled-full_parametric_y2p50_qt32_noprior_weightsall',  'W^{+}, |y| #leq 2.0, q_{T}<32 GeV, 30M events', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_random_parametric_y2p50_qt32_noprior_A1A3pol3_pt55',  'W^{+}, |y| #leq 2.5, q_{T}<32 GeV, 30M events, A_{1,3} pol_{3}, p_{T}<55', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov_parametric_y2p50_qt32_noprior_A1A3pol3_pt55',  'W^{+}, |y| #leq 2.5, q_{T}<32 GeV, 30M events, A_{1,3} pol_{3}, p_{T}<55', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated'],
        #['CC_FxFx_Wplus_asimov-scaled-full_parametric_y2p50_qt32_A1A3pol3_noprior_weightsall', '', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32_decorrelated']
        ['CC_FxFx_Wplus_random_parametric_y2p00_qt32v2_noprior_pt55', '', 'CC_FxFx_Wplus_WpreFSR_all_A0-4_forced_v4_finer_y_qt32v2_decorrelated']
        ]
    for fname in files:
        for v in ['biasANDrms']:
            profile_toys( fname=fname, alphas=['mass'],   var=v, postfix='_mass', save_pulls=True,  truth='val', do_fit=True)
            profile_toys( fname=fname, alphas=['minuit'], var=v, postfix='',      save_pulls=False, truth='', do_fit=True)
            profile_toys( fname=fname, alphas=['norm'],   var=v, postfix='_norm', save_pulls=True,  truth='val', do_fit=True)
            for coeff in ['A0', 'A1', 'A2', 'A3', 'A4']:
                profile_toys( fname=fname, alphas=[coeff], var=v, postfix='_'+coeff, save_pulls=True,  truth='fit', do_fit=True)

elif argv[1]=='derivatives':
    bins_template_y = [ 0.,  0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,  1.6, 1.8,  2.]
    #bins_template_y = [ 0.,  0.2 ]
    bins_template_qt= [ 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 32.0 ]
    #bins_template_qt= [ 0.0, 4.0 ]
    for iqt,qt in enumerate(bins_template_qt):
        for iy,y in enumerate(bins_template_y):
            derivative_templates(charge='Wplus', var='WpreFSR', coeff_eval='val', tag='finer_y_qt32', bin=(iqt,iy,1), bins_qt=[], bins_y=[])

elif argv[1]=='weighted_templates':
    weighted_templates(charge='Wplus', weights=range(109) )
                       
elif argv[1]=='profile_systs':
    fnames = [ 
        ['CC_FxFx_Wplus_asimov-scaled-in_norm_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #leq 2.0, q_{T}<32 GeV, 30M events: fit M_{W}'],
        ['CC_FxFx_Wplus_asimov-scaled-in_release_parametric_y2p00_qt32_noprior_weightsall',  'W^{+}, |y| #leq 2.0, q_{T}<32 GeV, 30M events: fit (#nu_{ij}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-in_norm_parametric_y2p00_qt32_noprior_weightsall',  'W^{+}, |y| #leq 2.0, q_{T}<32 GeV, 30M events: fit (A_{i}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-in_parametric_y2p00_qt32_noprior_weightsall',  'W^{+}, |y| #leq 2.0, q_{T}<32 GeV, 30M events: fit (#nu_{ij}, A_{i}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-out-qt_norm_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #leq 2.0, q_{T}>32 GeV, 30M events: fit M_{W}'],
        ['CC_FxFx_Wplus_asimov-scaled-out-qt_norm_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #leq 2.0, q_{T}>32 GeV, 30M events: fit (A_{i}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-out-qt_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #leq 2.0, q_{T}>32 GeV, 30M events: fit (#nu_{ij}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-out-qt_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #leq 2.0, q_{T}>32 GeV, 30M events: fit (#nu_{ij}, A_{k}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-out-y_norm_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 2.0, q_{T}#geq0 GeV, 30M events: fit M_{W}'], 
        ['CC_FxFx_Wplus_asimov-scaled-out-y_norm_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 2.0, q_{T}#geq0 GeV, 30M events: fit (A_{i}, M_{W})'],
        ['CC_FxFx_Wplus_asimov-scaled-out-y_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 2.0, q_{T}#geq0 GeV, 30M events: fit (#nu_{ij}, M_{W})'],  
        ['CC_FxFx_Wplus_asimov-scaled-out-y_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 2.0, q_{T}#geq0 GeV, 30M events: fit (#nu_{ij}, A_{k}, M_{W})'],  
        ['CC_FxFx_Wplus_asimov-scaled-full_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events: fit (#nu_{ij}, A_{k}, M_{W})'],   
        ['CC_FxFx_Wplus_asimov-scaled-full_norm_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events: fit (M_{W})'],   
        ['CC_FxFx_Wplus_asimov-scaled-full_norm_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events: fit (A_{k}, M_{W})'],   
        ['CC_FxFx_Wplus_asimov-scaled-full_release_parametric_y2p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events: fit (#nu_{ij}, M_{W})'],   
        ['CC_FxFx_Wplus_asimov-scaled-full_parametric_y2p50_qt32_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events: fit (#nu_{ij}, A_{k}, M_{W})'],   
        ['CC_FxFx_Wplus_asimov-scaled-full_parametric_y3p00_qt32_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events: fit (#nu_{ij}, A_{k}, M_{W})'],   
        ['CC_FxFx_Wplus_asimov-scaled-full_parametric_y2p50_qt32_A1A3pol3_noprior_weightsall', 'W^{+}, |y| #geq 0.0, q_{T} #geq 0 GeV, 30M events, A_{1,3} pol_{3}: fit (#nu_{ij}, A_{k}, M_{W})'],   
              ]
    for fname in fnames:
        profile_systs( fname=fname,  mass_true=80.419 )

else:
    print "Option not available"
