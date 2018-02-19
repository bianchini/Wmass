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

def fit_BW(tree=ROOT.TTree(), running=0, var='', cut='', tag=''):

    h = ROOT.TH1F('h', 'h', 140, 50., 110.)
    h.Sumw2()
    h.SetStats(0)
    tree.Draw(var+'_mass>>h', 'weights[0]*('+cut+')')

    fit = None
    if running:
        fit = ROOT.TF1('fit', '[0]*(x*x*[1]*[1]/(TMath::Power(x*x-[2]*[2],2) + x*x*x*x*([1]*[1])/([2]*[2])) )',  h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        fit.SetParameter(0, h.Integral()*1.0)
        fit.SetParameter(1,  2.000 )
        fit.SetParameter(2, 80.000 )
    else:
        fit = ROOT.TF1('fit', '[0]*( 1.0/(TMath::Power(x*x-[2]*[2],2) + [2]*[2]*[1]*[1]) )',  h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
        fit.SetParameter(0, h.Integral()*1000.)
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
    leg.AddEntry(fit, '#splitline{'+('Run.-width ' if running else '')+'BW M_{W}='+'{:05.3f}'.format(r.Parameter(2))+' #pm '+'{:05.3f}'.format(r.ParError(2))+', #Gamma_{W}='+'{:04.3f}'.format(r.Parameter(1))+' #pm '+'{:04.3f}'.format(r.ParError(1))+'}{#chi^{2}/ndof='+'{:02.0f}'.format(r.Chi2())+'/'+str(r.Ndf())+'}', 'L')
    leg.Draw()
    #raw_input()
    if tag!='':
        c.SaveAs('fit_BW_'+('run-width_' if running else 'nonrun-width_')+var+'_'+tag+'.png')
    c.IsA().Destructor( c )
    leg.IsA().Destructor( leg )
    return (r.Parameter(2), r.ParError(2), r.Parameter(1), r.ParError(1))

def run_fromPrompt( var='Wdress', running=0, cut='(mu_charge==+13 & isFromW==2)', tag='' ):
    infile = ROOT.TChain('tree')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_1.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_2.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_3.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_4.root')
    res = fit_BW(tree=infile, running=running, var=var, cut=cut, tag=tag)

def run_all(var='Wdress', running=0, tag='' ):

    infile = ROOT.TChain('tree')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_1.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_2.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_3.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_4.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_5.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_6.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_7.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_8.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_9.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_10.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_11.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_12.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_13.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_14.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_15.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_16.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_17.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_18.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_19.root')
    infile.Add('/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180214_160748/0000/tree_20.root')

    cuts = [
        'abs('+var+'_y)>=0.0 & abs('+var+'_y)<1.0 & mu_charge==+13',
        'abs('+var+'_y)>=1.0 & abs('+var+'_y)<2.0 & mu_charge==+13',
        'abs('+var+'_y)>=2.0 & abs('+var+'_y)<3.0 & mu_charge==+13',
        'abs('+var+'_y)>=3.0 & abs('+var+'_y)<4.0 & mu_charge==+13',
        'abs('+var+'_y)>=4.0 & abs('+var+'_y)<5.0 & mu_charge==+13',
        'abs('+var+'_y)>=0.0 & abs('+var+'_y)<1.0 & mu_charge==-13',
        'abs('+var+'_y)>=1.0 & abs('+var+'_y)<2.0 & mu_charge==-13',
        'abs('+var+'_y)>=2.0 & abs('+var+'_y)<3.0 & mu_charge==-13',
        'abs('+var+'_y)>=3.0 & abs('+var+'_y)<4.0 & mu_charge==-13',
        'abs('+var+'_y)>=4.0 & abs('+var+'_y)<5.0 & mu_charge==-13',
        ]

    width = 0.50
    x = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
    y = np.zeros(len(cuts))
    ye = np.zeros(len(cuts))
    for ic,c in enumerate(cuts):
        posttag = 'y'+str(x[ic%(len(cuts)/2)])+str(2*(ic/(len(cuts)/2))-1)
        res = fit_BW(tree=infile, running=running, var=var, cut=c, tag=(tag+'_'+posttag))
        y[ic] = res[0]
        ye[ic] = res[1]
    plt.figure()
    fig, ax = plt.subplots()
    ax.errorbar(x, y[:(len(cuts)/2)], xerr=width, yerr=ye[:(len(cuts)/2)], fmt='o', color='b', label='$W^-$')
    ax.errorbar(x, y[(len(cuts)/2):], xerr=width, yerr=ye[(len(cuts)/2):], fmt='o', color='g', label='$W^+$')
    ax.plot([0,5], [80.419,80.419], 'r--', label='$M_{W}$ in MC')
    plt.axis([x[0]-width, x[-1]+width, 80.300, 80.480])
    plt.grid(True)
    legend = ax.legend(loc='upper center', shadow=False, fontsize='x-large')
    plt.xlabel('$|y|$', fontsize=20)
    plt.ylabel('$M_{W}$ pole', fontsize=20)
    plt.title(var+' mass', fontsize=20)
    plt.show()
    plt.savefig('pole_'+('run-width_' if running else 'nonrun-width_')+var+'_mass_vs_y.png')
    plt.close()
    #infile.Close()

###################
#run_fromPrompt()
#run_all(var='lhe',    running=0, tag='')
#run_all(var='Wdress', running=0, tag='')
run_all(var='WpreFSR', running=0, tag='')
###################
