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
    # take scales as uncorrelated
    scale_ids = ([0,1,2,3,4,6,8]) 
    weights['scale'] = [[0,0]]*len(scale_ids)*len(scale_ids)
    count = 0
    for isid1,sid1 in enumerate(scale_ids):
        for isid2,sid2 in enumerate(scale_ids):
            if not veto_scale([sid1,sid2]):
                weights['scale'][count] = [sid1,sid2]
                count += 1
    weights['scale'] =  weights['scale'][0:count]

    # take pdfs as correlated
    pdf_ids = range(9,109)
    weights['pdf'] = [[0,0]]*len(pdf_ids)
    for isid1,sid1 in enumerate(pdf_ids):
        weights['pdf'][isid1] = [sid1,sid1]

    #weights['scale'] = ([0] + [1,2,3,4,6,8])
    #weights['pdf'] = range(9,109)
    for DY in ['CC_FxFx']:
        for q in ['Wplus']:
            for var in ['WpreFSR']:
                get_covariance(fname='../root/tree_histos1_'+DY+'.root', DY=DY, var=var, q=q, weights=weights, 
                               coefficients=['A0','A1','A2','A3','A4'], 
                               fix_to_zero={'A0': [0,1], 'A1': [0], 'A2': [0,1], 'A3': [0], 'A4': [], 'A5':[0], 'A6': [0], 'A7': [0]},
                               forced_orders={'A0': 3, 'A1': 4, 'A2': 3, 'A3': 4, 'A4': 4, 'A5': 0, 'A6':0, 'A7':0},
                               add_stat_uncert=True, postfix='all_A0-4_forced_v4_finer_y_qt_decorrelated',
                               save_corr=True, save_coeff=True, save_tree=True, save_pkl=True,
                               #np_bins_template_qt = np.array([   0.,  4., 8.,  12.,  16.,  20.,  24.,  32.,   40.,   60. ]),
                               np_bins_template_qt = np.array([  0.,  2.,   4.,  6.,   8.,  10.,   12.,  14.,  16., 20.,   24.,  32.,   40.,   60. ]),
                               np_bins_template_y  = np.array([-3.5, -3. , -2.5, -2., -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2,  0. , 
                                                                0.2,  0.4, 0.6, 0.8, 1.0, 1.2, 1.4,  1.6, 1.8,  2. ,  2.5,  3. ,  3.5]),
                               #np_bins_template_y  = np.array([-3.5, -3. , -2.5, -2., -1.6, -1.2, -0.8, -0.4,  0. , 0.4, 0.8, 1.2,  1.6,  2. ,  2.5,  3. ,  3.5]),
                               plot_updown=False,
                               decorrelate_scale=True
                               )
            
elif argv[1]=='cov':
    weights = {}

    scale_ids = ([0,1,2,3,4,6,8]) 
    weights['scale'] = [[0,0]]*len(scale_ids)*len(scale_ids)
    count = 0
    for isid1,sid1 in enumerate(scale_ids):
        for isid2,sid2 in enumerate(scale_ids):
            if not veto_scale([sid1,sid2]):
                weights['scale'][count] = [sid1,sid2]
                count += 1
    weights['scale'] =  weights['scale'][0:count]
    #weights['scale'] = ([0] + [1,2,3,4,6,8])
    #weights['pdf'] = range(9,109)
    #weights['scale'] = ([0,0], [1,1], [2,2], [3,3], [4,4], [5,5])
    #weights['pdf'] = [[9,9], [10,10], [11,11], [12,12], [13,13]]

    pdf_ids = range(9,109)
    weights['pdf'] = [[0,0]]*len(pdf_ids)
    for isid1,sid1 in enumerate(pdf_ids):
        weights['pdf'][isid1] = [sid1,sid1]

    for DY in ['CC_FxFx']:
        for q in ['Wplus']:
            for var in ['WpreFSR']:
                for coeff in [ 'A0'
                              #'A3','A4'
                              ]: 
                    get_covariance(fname='../root/tree_histos1_'+DY+'.root', DY=DY, var=var, q=q, weights=weights, 
                                   coefficients=[coeff], 
                                   fix_to_zero={'A0': [0,1], 'A1': [0], 'A2': [0,1], 'A3': [0], 'A4': [], 'A5':[0], 'A6': [0], 'A7': [0]},
                                   forced_orders={'A0': 3, 'A1': 4, 'A2': 3, 'A3': 4, 'A4': 4, 'A5': 0, 'A6':0, 'A7':0},
                                   #forced_orders={},
                                   add_stat_uncert=True, postfix=coeff+'_TEST_correlate',
                                   save_corr=True, save_coeff=True, save_tree=True, save_pkl=True,
                                   np_bins_template_qt = np.array([   0.,  4., 8.,  12.,  16.,  20.,  24.,  32.,   40.,   60. ]),
                                   np_bins_template_y  = np.array([-3.5, -3. , -2.5, -2., -1.6, -1.2, -0.8, -0.4,  0. ,  
                                                                    0.4, 0.8, 1.2,  1.6,  2. ,  2.5,  3. ,  3.5]),
                                   #np_bins_template_y  = np.array([ -0.4,  0. , 0.4])
                                   plot_updown=False,
                                   decorrelate_scale=False
                                   )

elif argv[1]=='closure':
    for DY in ['CC_FxFx']:
        for q in ['Wplus']:
            for var in ['WpreFSR']:
                for ceval in ['val']:
                    plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, 
                                      coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                                      byInvPdf=False,
                                      save_2D=True, save_pdf=False, save_summary=True, 
                                      do_toy=False, extra_variance_toy=0.0,
                                      do_pt_eta=True)

                    
                    # MC
                    #plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                    #                  byInvPdf=False,
                    #                  save_2D=False, save_pdf=False, save_summary=True, 
                    #                  do_toy=False, extra_variance_toy=0.0)
                    # toy
                    #plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                    #                  byInvPdf=False,
                    #                  save_2D=False, save_pdf=False, save_summary=True, 
                    #                  do_toy=True, extra_variance_toy=0.001)
                    # MC byInvPdf
                    #plot_closure_test(charge=q,  DY=DY, var=var, coeff_eval=ceval, coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'], 
                    #                  byInvPdf=True,
                    #                  save_2D=True, save_pdf=False, save_summary=True, 
                    #                  do_toy=False, extra_variance_toy=0.0)            

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

elif argv[1]=='templates':
    #np_bins_template_qt = np.array([   0.,     4.,   8.,   12.,   16., 20.,   24.,  32.,   40.,   60. ])
    np_bins_template_qt = np.array([   0.,  2.,   4.,  6.,   8.,  10.,   12.,  14.,  16., 20.,   24.,  32.,   40.,   60. ])
    #np_bins_template_y  = np.array([-3.5, -3. , -2.5, -2., -1.6, -1.2, -0.8, -0.4,  0. ,  0.4, 0.8, 1.2,  1.6,  2. ,  2.5,  3. ,  3.5])
    np_bins_template_y  = np.array([-3.5, -3. , -2.5, -2., -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2,  0. , 
                                     0.2,  0.4, 0.6, 0.8, 1.0, 1.2, 1.4,  1.6, 1.8,  2. ,  2.5,  3. ,  3.5]),

    merge_templates(charges=['Wplus'], var=['WpreFSR'], coeff_eval=['val'], 
                    masses=[79.919, 80.419, 80.919], coeff=['A0', 'A1', 'A2', 'A3', 'A4'], 
                    np_bins_template_qt=np_bins_template_qt, 
                    np_bins_template_y=np_bins_template_y,
                    rebin=(2,4),
                    input_tag='CC_FxFx_extra_masses',
                    postfix='_finer_y_qt'
                    )

    #merge_templates(charges=['Wminus'], var=['WpreFSR'], coeff_eval=['val'], 
    #                masses=[80.419], coeff=['A0', 'A1', 'A2', 'A3', 'A4'], 
    #                np_bins_template_qt=np_bins_template_qt, 
    #                np_bins_template_y=np_bins_template_y,
    #                rebin=(2,4)
    #                )

