from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 8,7

import math

def fit_BW(tree=ROOT.TTree(), rel=1, var='', cut='', tag=''):

    h = ROOT.TH1F('h', 'h', 70, 50., 110.)
    h.Sumw2()
    h.SetStats(0)
    tree.Draw(var+'_mass>>h', 'weight*('+cut+')')

    fit = None
    if rel:
        fit = ROOT.TF1('fit', '[0]*([2]*[2]/(TMath::Power(x*x-[2]*[2],2) + x*x*x*x*([1]*[1])/([2]*[2])) )',  50., 110.)
        fit.SetParameter(0, h.Integral())
        fit.SetParameter(1,  2.000 )
        fit.SetParameter(2, 80.000 )
    else:
        fit = ROOT.TF1('fit', '[0]*([1]*[1]/4./(TMath::Power(x-[2],2) + [1]*[1]/4) )',  50., 110.)
        fit.SetParameter(0, h.Integral())
        fit.SetParameter(1,  2.000 )
        fit.SetParameter(2, 80.000 )

    fit.SetNpx(10000)
    r = h.Fit('fit', 'S')
    print 'Chi2 = ', r.Chi2()
    r.Print('V')

    c = ROOT.TCanvas("c", "canvas", 600, 600) 
    ROOT.gPad.SetLogy()
    leg = ROOT.TLegend(0.10,0.75,0.55,0.88, "","brNDC")
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)
    leg.SetFillColor(10)
    h.SetMaximum(h.GetMaximum()*10)
    h.SetTitle(cut)
    h.SetXTitle('M_{W} [GeV]')
    h.Draw()
    leg.AddEntry(h, var+' M_{W}', 'P')
    leg.AddEntry(fit, '#splitline{'+('Rel. ' if rel else 'Non-rel. ')+'BW M_{W}='+'{:05.3f}'.format(r.Parameter(2))+' #pm '+'{:05.3f}'.format(r.ParError(2))+', #Gamma_{W}='+'{:04.3f}'.format(r.Parameter(1))+' #pm '+'{:04.3f}'.format(r.ParError(1))+'}{#chi^{2}/ndof='+'{:02.0f}'.format(r.Chi2())+'/'+str(r.Ndf())+'}', 'L')
    leg.Draw()
    #raw_input()
    if tag!='':
        c.SaveAs('fit_BW_'+('rel_' if rel else 'non-rel_')+var+'_'+tag+'.png')
    c.IsA().Destructor( c )
    leg.IsA().Destructor( leg )
    return (r.Parameter(2), r.ParError(2), r.Parameter(1), r.ParError(1))

def run_fromPrompt():
    infile = ROOT.TFile('./tree_'+argv[3]+'.root')
    tree = infile.Get('tree')
    res = fit_BW(tree=tree, rel= int(argv[1]), var=argv[2], cut=argv[4], tag=argv[5])

def run_all():
    infile = ROOT.TFile('./tree_'+argv[3]+'.root')
    tree = infile.Get('tree')

    cuts = [
        'abs(lhe_y)>=0.0 & abs(lhe_y)<1.0 & charge<0',
        'abs(lhe_y)>=1.0 & abs(lhe_y)<2.0 & charge<0',
        'abs(lhe_y)>=2.0 & abs(lhe_y)<3.0 & charge<0',
        'abs(lhe_y)>=3.0 & abs(lhe_y)<4.0 & charge<0',
        'abs(lhe_y)>=4.0 & abs(lhe_y)<5.0 & charge<0',
        'abs(lhe_y)>=0.0 & abs(lhe_y)<1.0 & charge>0',
        'abs(lhe_y)>=1.0 & abs(lhe_y)<2.0 & charge>0',
        'abs(lhe_y)>=2.0 & abs(lhe_y)<3.0 & charge>0',
        'abs(lhe_y)>=3.0 & abs(lhe_y)<4.0 & charge>0',
        'abs(lhe_y)>=4.0 & abs(lhe_y)<5.0 & charge>0',
        #'dressFSR_qt>=0.0 & dressFSR_qt<10 & charge>0', 
        #'dressFSR_qt>=10.0 & dressFSR_qt<20 & charge>0', 
        #'dressFSR_qt>=20.0 & dressFSR_qt<200 & charge>0', 
        #'dressFSR_qt>=0.0 & dressFSR_qt<10 & charge<0', 
        #'dressFSR_qt>=10.0 & dressFSR_qt<20 & charge<0', 
        #'dressFSR_qt>=20.0 & dressFSR_qt<200 & charge<0', 
        ]

    width = 0.5
    x = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
    #x = np.array([0.5,1.5,2.5])
    y = np.zeros(len(cuts))
    ye = np.zeros(len(cuts))
    for ic,c in enumerate(cuts):
        tag = 'y'+str(x[ic%5])+str(2*(ic/5)-1)
        res = fit_BW(tree=tree, rel=int(argv[1]), var=argv[2], cut=c, tag=tag)
        y[ic] = res[0]
        ye[ic] = res[1]
    plt.figure()
    fig, ax = plt.subplots()
    ax.errorbar(x, y[:(len(cuts)/2)], xerr=width, yerr=ye[:(len(cuts)/2)], fmt='o', color='b', label='$W^-$')
    ax.errorbar(x, y[(len(cuts)/2):], xerr=width, yerr=ye[(len(cuts)/2):], fmt='o', color='g', label='$W^+$')
    plt.axis([x[0]-width, x[-1]+width, 80.250, 80.500])
    plt.grid(True)
    legend = ax.legend(loc='upper center', shadow=False, fontsize='x-large')
    plt.xlabel('$|y|$', fontsize=20)
    plt.ylabel('$M_{W}$ pole', fontsize=20)
    plt.title(argv[2], fontsize=20)
    plt.show()
    plt.savefig('pole_'+argv[2]+'_mass_vs_y.png')
    plt.close()
    infile.Close()

###################
run_fromPrompt()
#run_all()
###################
