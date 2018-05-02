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

from fit_utils import *
from plot_utils import *

if argv[1]=='pole':
    for var in ['lhe', 'Wdress', 'WpreFSR' ]:
        get_pole_mass(var=var, running=0, tag='qt')
        get_pole_mass(var=var, running=0, tag='y')

elif argv[1]=='plot_qt':
    plot_qt(var='Wdress')

elif argv[1]=='pt_bias':
    qt_max = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
    #pt_bias_vs_qt( q='Wplus', var='Wdress', qt_max=qt_max, pt_max=65.0, pt_min=25.0, eta_max=2.5, verbose=False, debug=False )
    print_pt_bias_vs_qt(pt_min=25.0, pt_max=[55.0, 60.0, 65.0, 70.0], qt_max=qt_max, q='Wplus', var='Wdress')

elif argv[1]=='fit':
    for coeff in ['A0','A1','A2','A3','A4','A5', 'A6', 'A7']:
        draw_y_slice(fname='../root/tree_histos_CC.root', var='Wdress', coeff=coeff, weight_name=0, do_fit=True)
        #draw_qt_slice(fname='../root/tree_histos_NC.root', var='Wdress', coeff=coeff, weight_name=0, do_fit=False)

elif argv[1]=='cov-all':
    weights = {}
    weights['scale'] = ([0] + [1,2,3,4,6,8])
    weights['pdf'] = range(9,109)
    for DY in ['CC_FxFx', 'CC_MG5']:
        for q in ['Wplus', 'Wminus']:
            for var in ['Wdress', 'Wbare', 'WpreFSR']:
                get_covariance(fname='../root/tree_histos1_'+DY+'.root', DY=DY, var=var, q=q, weights=weights, 
                               coefficients=['A0','A1','A2','A3','A4','A5','A6','A7'], 
                               add_stat_uncert=True, postfix='all_A0-7',
                               save_corr=True, save_coeff=True, save_tree=True, save_pkl=True)
            
elif argv[1]=='cov':
    weights = {}
    weights['scale'] = ([0] + [1,2,3,4,6,8])
    weights['pdf'] = range(9,109)
    for DY in ['CC_FxFx', 'CC_MG5']:
        for q in ['Wplus']:
            for var in ['Wdress', 'Wbare', 'WpreFSR']:
                for coeff in ['A0']: #,'A1','A2','A3','A4','A5','A6','A7']: 
                    get_covariance(fname='../root/tree_histos1_'+DY+'.root', DY=DY, var=var, q=q, weights=weights, 
                                   coefficients=[coeff], 
                                   add_stat_uncert=True, postfix=coeff,
                                   save_corr=True, save_coeff=True, save_tree=True, save_pkl=True)

elif argv[1]=='closure':
    for DY in ['CC_FxFx', 'CC_MG5']:
        for q in ['Wplus', 'Wminus']:
            for var in ['Wdress', 'Wbare', 'WpreFSR']:
                for ceval in ['val']:
                    # MC
                    plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                                      byInvPdf=False,
                                      save_2D=True, save_pdf=False, save_summary=False, 
                                      do_toy=False, extra_variance_toy=0.0)
                    # toy
                    plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                                      byInvPdf=False,
                                      save_2D=False, save_pdf=False, save_summary=True, 
                                      do_toy=True, extra_variance_toy=0.001)
                    # MC byInvPdf
                    plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                                      byInvPdf=True,
                                      save_2D=True, save_pdf=False, save_summary=True, 
                                      do_toy=False, extra_variance_toy=0.0)            

elif argv[1]=='compare':
    compare_fit_results(DYs=['CC_FxFx'], 
                        charges=['Wplus'], 
                        variables=['Wbare', 'Wdress', 'WpreFSR'], 
                        coefficients=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                        fit_range=[0.0, 50.0], postfix='all_A0-7', 
                        compare='CC_FxFx_Wplus_dressing')
    compare_fit_results(DYs=['CC_FxFx'], 
                        charges=['Wminus'], 
                        variables=['Wbare', 'Wdress', 'WpreFSR'], 
                        coefficients=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                        fit_range=[0.0, 50.0], postfix='all_A0-7', 
                        compare='CC_FxFx_Wminus_dressing')
    compare_fit_results(DYs=['CC_FxFx', 'CC_MG5'], 
                        charges=['Wplus'], 
                        variables=['Wdress'], 
                        coefficients=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                        fit_range=[0.0, 50.0], postfix='all_A0-7', 
                        compare='CC_Wplus_Wdress_generator')
    compare_fit_results(DYs=['CC_FxFx', 'CC_MG5'], 
                        charges=['Wminus'], 
                        variables=['Wdress'], 
                        coefficients=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                        fit_range=[0.0, 50.0], postfix='all_A0-7', 
                        compare='CC_Wminus_Wdress_generator')
