import os.path
from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import math

import sys
sys.path.append('./')
sys.path.append('../python/')

import pickle
from pprint import pprint
from template_parameters import params_test as params
from array import array
import numpy as np
ROOT.TGaxis.SetMaxDigits(2)

job =  '1e7_corr_cube_rebin4-1'
j = 1
toy = 0
dir_name = 'TEST'
params['params_template']['pt'] = np.array([0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 26.0, 32.0])
params['params_template']['y']  = np.array([0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4])
input_shapes_y  = params['params_template']['y']
input_shapes_pt  = params['params_template']['pt']
os.system('mkdir plots/'+dir_name)

if not os.path.isfile('crab/crab_unfolder_'+job+'/results/result_'+job+'_'+str(j)+'.pkl'):
    exit(1)
f = open('crab/crab_unfolder_'+job+'/results/result_'+job+'_'+str(j)+'.pkl', 'r')
res = pickle.load(f)

bins = array( 'f',  input_shapes_pt )

def plot_dsigma_dpt(res, pt_bin_name, y_bin_name, toy):
    fit = res[pt_bin_name+'_'+y_bin_name]['fit'][toy] 
    err = res[pt_bin_name+'_'+y_bin_name]['err'][toy] 
    true = res[pt_bin_name+'_'+y_bin_name]['true'][toy]    
    return (fit, err, true)

def plot_A(res, pt_bin_name, pt, y_bin_name, toy,  coeff, deg):
    val = 0.0
    err2 = 0.0
    if coeff<4:        
        for d in range(deg):
            val += math.pow(pt, d+2)*res['coeff'+str(d)+'_'+y_bin_name+'_A'+str(coeff)]['fit'][toy] 
            err2 += math.pow(math.pow(pt, d+2)*res['coeff'+str(d)+'_'+y_bin_name+'_A'+str(coeff)]['err'][toy],2) 
    else:
        for d in range(deg+1):
            if d==0:
                val += res['coeff'+str(d)+'_'+y_bin_name+'_A'+str(coeff)]['fit'][toy]
                err2 += math.pow(res['coeff'+str(d)+'_'+y_bin_name+'_A'+str(coeff)]['err'][toy],2)
                continue
            val += math.pow(pt, d+1)*res['coeff'+str(d)+'_'+y_bin_name+'_A'+str(coeff)]['fit'][toy] 
            err2 += math.pow(math.pow(pt, d+1)*res['coeff'+str(d)+'_'+y_bin_name+'_A'+str(coeff)]['err'][toy],2) 
    return (val, math.sqrt(err2), 0.0)

    err = res[pt_bin_name+'_'+y_bin_name]['err'][toy] 
    true = res[pt_bin_name+'_'+y_bin_name]['true'][toy]    
    return (fit, err, true)

for iy in range(len(input_shapes_y)-1):
    y_bin=[ input_shapes_y[iy], input_shapes_y[iy+1] ]
    y_bin_name = 'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
    print y_bin_name

    ROOT.TGaxis.SetMaxDigits(2)
    p_fit = ROOT.TH1F('pt_'+y_bin_name+'_fit', '', len(bins)-1, bins)
    p_gen = ROOT.TH1F('pt_'+y_bin_name+'_gen', '', len(bins)-1, bins)

    for ipt in range(len(input_shapes_pt)-1):
        pt_bin=[ input_shapes_pt[ipt], input_shapes_pt[ipt+1] ]
        pt_bin_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])

        fit, err, true = 0.0, 0.0, 0.0
        if argv[1]=='pt':
            (fit, err, true) = plot_dsigma_dpt(res, pt_bin_name, y_bin_name, toy)
        elif 'A' in argv[1]:
            (fit, err, true) = plot_A(res, pt_bin_name, 0.5*(pt_bin[0]+pt_bin[1]), y_bin_name, toy, int(argv[1][1]), (1 if 'quad' in job else 2) )
        
        p_fit.SetBinContent(ipt+1, fit )    
        p_fit.SetBinError(ipt+1, err)
        if argv[1]=='pt':
            print "\t --- ", pt_bin_name, " ==> accuracy = ", '{:03.1f}'.format((err/abs(true)*1e+02) if abs(true) > 0.0 else 0.0), "%"
        else:
            print "\t ", argv[1] , "--- ", pt_bin_name, " ==> accuracy = ", '{:03.1f}'.format(err*1e+02), "%"            
        p_gen.SetBinContent(ipt+1, true )
    
    c = ROOT.TCanvas("c_"+job+"_"+'pt_'+y_bin_name, "canvas", 600, 600) 
    leg = ROOT.TLegend(0.55,0.70,0.88,0.88, "","brNDC");
    leg.SetHeader('y_{W} #in ['+'{:03.2f}'.format(y_bin[0])+','+'{:03.2f}'.format(y_bin[1])+']')
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.05);
    leg.SetFillColor(10);

    title = 'W_pt_'+y_bin_name
    var = 'p_{T}^{W} [GeV]'
    p_fit.SetStats(0) 
    if not 'decorr' in job:
        p_fit.SetMinimum(0.0)
    if 'A' in argv[1]:
        p_fit.SetMinimum(-10.0)
        p_fit.SetMaximum(+10.0)

    p_fit.SetTitle('')
    p_fit.GetXaxis().SetTitle(var.title())
    p_fit.GetXaxis().SetTitleSize(25)
    p_fit.GetXaxis().SetTitleFont(43)
    p_fit.GetXaxis().SetTitleOffset(1.0)
    p_fit.GetXaxis().SetLabelFont(43) 
    p_fit.GetXaxis().SetLabelSize(20)
    if argv[1]=='pt':
        p_fit.GetYaxis().SetTitle('events')
    else:
        p_fit.GetYaxis().SetTitle(argv[1])
    p_fit.GetYaxis().SetTitleSize(25)
    p_fit.GetYaxis().SetTitleFont(43)
    p_fit.GetYaxis().SetTitleOffset(1.0)
    p_fit.GetYaxis().SetLabelFont(43) 
    p_fit.GetYaxis().SetLabelSize(20)
    p_fit.SetLineWidth(4)
    p_fit.SetLineStyle(ROOT.kSolid)
    p_fit.SetLineColor(ROOT.kRed)
    p_gen.SetLineWidth(4)
    p_gen.SetLineStyle(ROOT.kDashed)
    p_gen.SetLineColor(ROOT.kBlack)
    leg.AddEntry(p_gen, 'True', 'l')
    leg.AddEntry(p_fit, 'Fit', 'lp')
    p_fit.Draw("HISTE")
    p_gen.Draw("HISTSAME")
    leg.Draw()
    #raw_input()
        
    c.SaveAs('plots/'+dir_name+'/'+job+'_'+argv[1]+'_'+y_bin_name+'.png')
    c.IsA().Destructor( c )
