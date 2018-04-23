import ROOT
import math
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 8,7
from array import array
import pickle
from scipy import stats

import os.path

def readFiles(num_files=2):    
    infile = ROOT.TChain('tree')
    path = '/gpfs/ddn/srm/cms/store/user/bianchi/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/TEST/180220_102334/0000/'
    for i in range(0, num_files):
        if not os.path.isfile(path+'/tree_'+str(i)+'.root'):
            continue
        print 'Adding file...tree_'+str(i)+'.root'
        infile.Add(path+'/tree_'+str(i)+'.root')
    print 'Total entries:', infile.GetEntries()
    return infile


def get_pole_mass(var='Wdress', running=0, tag='y' ):

    import os.path
    from fit_utils import fit_BreitWigner

    infile = readFiles(num_files=20)

    cuts = []
    if tag=='y':
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
    else:
        cuts = [
            'abs('+var+'_qt)>= 0.0 & abs('+var+'_qt)<  5.0 & mu_charge==+13',
            'abs('+var+'_qt)>= 5.0 & abs('+var+'_qt)< 10.0 & mu_charge==+13',
            'abs('+var+'_qt)>=10.0 & abs('+var+'_qt)< 15.0 & mu_charge==+13',
            'abs('+var+'_qt)>=15.0 & abs('+var+'_qt)< 20.0 & mu_charge==+13',
            'abs('+var+'_qt)>=20.0 & abs('+var+'_qt)<200.0 & mu_charge==+13',
            'abs('+var+'_qt)>= 0.0 & abs('+var+'_qt)<  5.0 & mu_charge==-13',
            'abs('+var+'_qt)>= 5.0 & abs('+var+'_qt)< 10.0 & mu_charge==-13',
            'abs('+var+'_qt)>=10.0 & abs('+var+'_qt)< 15.0 & mu_charge==-13',
            'abs('+var+'_qt)>=15.0 & abs('+var+'_qt)< 20.0 & mu_charge==-13',
            'abs('+var+'_qt)>=20.0 & abs('+var+'_qt)<200.0 & mu_charge==-13',
            ]

    width = 0.50
    x = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
    y = np.zeros(len(cuts))
    ye = np.zeros(len(cuts))
    for ic,c in enumerate(cuts):
        posttag = 'y'+str(x[ic%(len(cuts)/2)])+str(2*(ic/(len(cuts)/2))-1)
        res = fit_BreitWigner(tree=infile, running=running, var=var, cut=c, tag=(tag+'_'+posttag))
        y[ic] = res[0]
        ye[ic] = res[1]
    plt.figure()
    fig, ax = plt.subplots()
    ax.errorbar(x, y[:(len(cuts)/2)], xerr=width, yerr=ye[:(len(cuts)/2)], fmt='o', color='b', label='$W^+$')
    ax.errorbar(x, y[(len(cuts)/2):], xerr=width, yerr=ye[(len(cuts)/2):], fmt='o', color='g', label='$W^-$')
    ax.plot([0,5], [80.419,80.419], 'r--', label='$M_{W}$ in MC')
    plt.axis([x[0]-width, x[-1]+width, 80.300, 80.480])
    plt.grid(True)
    legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
    if tag=='y':
        plt.xlabel('$|y|$', fontsize=20)
    else:
        plt.xlabel('$q_{T}$ bin', fontsize=20)
    plt.ylabel('$M_{W}$ pole', fontsize=20)
    plt.title(var+' mass', fontsize=20)
    plt.show()
    plt.savefig('plots/pole_'+('run-width_' if running else 'nonrun-width_')+var+'_mass_vs_'+tag+'.png')
    plt.close()
    #infile.Close()

def plot_qt_slice(tree=ROOT.TTree(), h=None, it=0, var='Wdress', plot='qt', cut='mu_charge==+13', label='', leg=None, add_overflow=True, normal=1.0):
    h.Sumw2()
    h.SetStats(0)
    h.SetLineWidth(3 if plot =='qt' else 2)
    h.SetLineStyle(ROOT.kSolid if plot == 'qt' else ROOT.kDashed)
    h.SetLineColor(it+1)
    tree.Draw(var+'_'+plot+'>>h_'+var+'_'+plot+'_'+str(it), 'weights[0]*('+cut+')')
    norm = h.Integral(0, h.GetNbinsX()+1, "")
    if add_overflow:
        h.SetBinContent(h.GetNbinsX(), h.GetBinContent(h.GetNbinsX()) + h.GetBinContent(h.GetNbinsX()+1) )
    #for ib,b in enumerate(bin_widths):
    #    binc = h.GetBinContent(ib+1)
    #     h.SetBinContent(ib+1, binc/b )
    print cut, h.Integral()
    h.SetTitle('')
    h.SetXTitle('q_{T} or p_{T} [GeV]')    
    if label==None:
        leg.AddEntry(h, '{:.2E}'.format(h.Integral()/normal),  'PL')
    else:
        leg.AddEntry(h, label, 'PL')
    return norm


def plot_qt(q='Wplus', var='Wdress', qt_max=50, bin_size=2.0, pt_max=60, pt_min=24):

    bins_qt = array( 'f',  np.linspace(0.0, qt_max+bin_size, int(qt_max/bin_size)+2) )
    bins_pt = array( 'f',  np.linspace(0.0, pt_max, int(pt_max/bin_size)+1) )

    infile = readFiles(num_files=20)

    charge_cut = 'mu_charge==-13' if q=='Wplus' else 'mu_charge==+13'
    cuts = [
        'abs('+var+'_y)>=0.0 & abs('+var+'_y)<1.0 &'+charge_cut,
        'abs('+var+'_y)>=1.0 & abs('+var+'_y)<2.0 &'+charge_cut,
        'abs('+var+'_y)>=2.0 & abs('+var+'_y)<3.5 &'+charge_cut,
        'abs('+var+'_y)>=3.0 & abs('+var+'_y)<10 &'+charge_cut,
        ]

    labels = [
        '|y|<1.0',
        '1.0#leq|y|<2.0',
        '2.0#leq|y|<3.5',
        '3.5#leq|y|',
        ]

    leg1 = ROOT.TLegend(0.55,0.60,0.85,0.88, "","brNDC")
    leg1.SetFillStyle(0)
    leg1.SetBorderSize(0)
    leg1.SetTextSize(0.03)
    leg1.SetFillColor(10)
    leg2 = ROOT.TLegend(0.55,0.60,0.85,0.88, "","brNDC")
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.SetTextSize(0.03)
    leg2.SetFillColor(10)

    histos = {}
    for icut,cut in enumerate(cuts):
        histos['qt_'+str(icut)]    = ROOT.TH1F('h_'+var+'_qt_'+str(icut),    '', len(bins_qt)-1, bins_qt)        
        histos['mu_pt_'+str(icut)] = ROOT.TH1F('h_'+var+'_mu_pt_'+str(icut), '', len(bins_pt)-1, bins_pt)                
        norm_qt    = plot_qt_slice(tree=infile, h=histos['qt_'+str(icut)],    it=icut, var=var, plot='qt',    cut=cut, label='q_{T}: '+labels[icut], leg=leg1, add_overflow=True)
        norm_mu_pt = plot_qt_slice(tree=infile, h=histos['mu_pt_'+str(icut)], it=icut, var=var, plot='mu_pt', cut=cut+' & '+var+'_qt>'+str(qt_max), leg=leg2, label=None, normal=norm_qt, add_overflow=False)
        histos['qt_'+str(icut)].Scale(1./norm_qt)
        histos['mu_pt_'+str(icut)].Scale(1./norm_qt)

    c = ROOT.TCanvas("c", "canvas", 1200, 600)
    c.Divide(2,1)

    for icut,cut in enumerate(cuts):
        if icut==0:            
            c.cd(1)
            histos['qt_'+str(icut)].SetTitle(q)
            histos['qt_'+str(icut)].SetMaximum(10.0)
            histos['qt_'+str(icut)].SetMinimum(1e-03)
            histos['qt_'+str(icut)].GetXaxis().SetRangeUser(0.0, qt_max+bin_size)
            ROOT.gPad.SetLogy()
            histos['qt_'+str(icut)].Draw("HISTE")
            c.cd(2)
            histos['mu_pt_'+str(icut)].SetTitle(q+': q_{T}>'+str(qt_max))
            histos['mu_pt_'+str(icut)].SetMaximum(1e-01)
            histos['mu_pt_'+str(icut)].SetMinimum(1e-04)
            histos['mu_pt_'+str(icut)].GetXaxis().SetRangeUser(pt_min, pt_max)
            ROOT.gPad.SetLogy()
            histos['mu_pt_'+str(icut)].Draw("HISTE")
        else:
            c.cd(1)
            histos['qt_'+str(icut)].Draw("HISTESAME")            
            c.cd(2)
            histos['mu_pt_'+str(icut)].Draw("HISTESAME")
    c.cd(1)
    leg1.Draw()
    c.cd(2)
    leg2.Draw()

    c.cd()
    c.Update()
    c.Draw()
    c.SaveAs('plots/'+q+'_qt_spectrum_vs_y_'+'{:02.0f}'.format(qt_max)+'.png')
    c.IsA().Destructor( c )
    leg1.IsA().Destructor( leg1 )
    leg2.IsA().Destructor( leg2 )


def pt_bias_vs_qt( q='Wplus', var='Wdress', qt_max=[], pt_max=60, pt_min=24, eta_max=2.5, verbose=False, debug=False ):

    infile = readFiles(num_files=5)

    weights = {}
    if debug:        
        weights['scale'] = ([0, 1])
        weights['pdf']   = ([0] + range(9,10))
    else:
        weights['scale'] = [0]+[1,2,3,4,6,8]
        weights['pdf']   = [0]+range(9,109)

    bins_pt = array( 'f',  np.linspace(pt_min, pt_max, 50) )
    h = ROOT.TH1F('h', '', len(bins_pt)-1, bins_pt) 
    h.Sumw2()

    cut  = 'mu_charge==-13' if q=='Wplus' else 'mu_charge==+13'
    cut += ' & '+var+'_mu_pt>='+'{:03.1f}'.format(pt_min)+' & '+var+'_mu_pt<'+'{:03.1f}'.format(pt_max)
    cut += ' & TMath::Abs('+var+'_mu_eta)<='+'{:03.1f}'.format(eta_max)
    banner = 'pt in ['+'{:05.3f}'.format(pt_min)+', '+'{:05.3f}'.format(pt_max)+']: '

    shifts = {}
    for qt in qt_max:
        cut_out = (cut + ' & '+var+'_qt>='+'{:03.1f}'.format(qt) )
        cut_in  = (cut + ' & '+var+'_qt <'+'{:03.1f}'.format(qt) )
        h.Reset()
        infile.Draw(var+'_mu_pt'+'>>h', 'weights['+str(0)+']*('+cut_in+')')
        norm_in = h.Integral()
        h.Reset()
        infile.Draw(var+'_mu_pt'+'>>h', 'weights['+str(0)+']*('+cut_out+')')
        norm_out = h.Integral()
        mean_out = h.GetMean()
        mean_out_err = h.GetMeanError()
        shifts['qt{:02.0f}'.format(qt)+'_'+'mean']     = mean_out
        shifts['qt{:02.0f}'.format(qt)+'_'+'mean_err'] = mean_out_err
        frac = norm_out/(norm_out+norm_in)
        shifts['qt{:02.0f}'.format(qt)+'_'+'out_frac'] = frac
        for syst in ['pdf', 'scale']:
            if verbose:
                print 'Doing syst:', syst
            data = np.zeros( (len(weights[ syst ])) )
            for iw,w in enumerate(weights[ syst ]):
                h.Reset()
                infile.Draw(var+'_mu_pt'+'>>h', 'weights['+str(w)+']*('+cut_out+')')
                mean = h.GetMean()  
                mean_err = h.GetMeanError()  
                if verbose:
                    print '\t mean('+str(w)+') = '+'{:05.3f}'.format(mean), ' +/- ', '{:05.3f}'.format(mean_err)
                data[iw] = mean
            shifts['qt{:02.0f}'.format(qt)+'_'+syst+'_std'] = np.std(data) 
            shifts['qt{:02.0f}'.format(qt)+'_'+syst+'_max'] = max( abs(np.amax(data-mean_out)), abs(np.amin(data-mean_out)) )
            if verbose:
                print syst, 'std: '+'{:05.3f}'.format( shifts['qt{:02.0f}'.format(qt)+'_'+syst+'_std'] ), ', max: ', '{:05.3f}'.format( shifts['qt{:02.0f}'.format(qt)+'_'+syst+'_max'])

        err2 = 0.0
        err2 += math.pow(shifts['qt{:02.0f}'.format(qt)+'_'+'pdf_std'],  2.0)
        err2 += math.pow(shifts['qt{:02.0f}'.format(qt)+'_'+'scale_max'],2.0)
        err = math.sqrt(err2)
        print 'qt >= '+'{:05.3f}'.format(qt)+', '+banner+'{:05.3f}'.format(err) + '*' + '{:04.3f}'.format(frac) + ' = ' + '{:05.3f}'.format(err*frac)
        shifts['qt{:02.0f}'.format(qt)+'_'+'tot'] = err*frac
        
    pickle.dump(shifts, open('plots/'+q+'_'+var+'_pt'+'{:02.0f}'.format(pt_min)+'-'+'{:02.0f}'.format(pt_max)+'_bias_vs_qt.pkl','wb') )
    
    h.IsA().Destructor( h )

def print_pt_bias_vs_qt(pt_min=25.0, pt_max=[], qt_max=[], q='Wplus', var='Wdress'):
    plt.figure()
    fig, ax = plt.subplots()
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    fmts = ['o', 'v', '^', '>', '<', '.']
    for ipt, pt in enumerate(pt_max):
        f = open('plots/'+q+'_'+var+'_pt'+'{:02.0f}'.format(pt_min)+'-'+'{:02.0f}'.format(pt)+'_bias_vs_qt.pkl')
        res = pickle.load(f)
        x = np.zeros(len(qt_max))
        y = np.zeros(len(qt_max))
        ye = np.zeros(len(qt_max))
        for iqt,qt in enumerate(qt_max):
            y[iqt]  = res['qt{:02.0f}'.format(qt)+'_tot']/res['qt{:02.0f}'.format(qt)+'_mean']
            ye[iqt] = res['qt{:02.0f}'.format(qt)+'_mean_err']/res['qt{:02.0f}'.format(qt)+'_mean']
            x[iqt] = qt_max[iqt]
        ax.errorbar(x, y, xerr=0.0, yerr=0.0, fmt=fmts[ipt], color=colors[ipt], label='$p_{T}<'+'{:02.0f}'.format(pt)+'$ GeV', markersize=12)
    plt.axis([qt_max[0]-5.0, qt_max[-1]+5.0, 1e-05, 1e-02 ])
    legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
    plt.yscale('log')
    plt.xlabel('$q_{T}^{max}$ GeV', fontsize=20)
    plt.ylabel('$\Delta p_{T}/p_{T}$', fontsize=20)
    plt.title('bias', fontsize=20)
    plt.grid(True)
    plt.show()
    plt.savefig('test.png')
    plt.close()


def plot_closure_test(charge='Wplus', var='Wdress', coeff_eval='val', min_val = -999., save_2D=True, save_summary=True):

    print "*** plot_closure_test ***"

    from pylab import rcParams
    rcParams['figure.figsize'] = 14,8

    from tree_utils import np_bins_qt, np_bins_y

    #np_bins_qt = np.array([0.0, 2.0])
    np_bins_y_p0 = np.linspace(-4.0, -2.5,  4)
    np_bins_y_p1 = np.linspace(-2.0, +2.0, 21)
    np_bins_y_p2 = np.linspace(+2.5, +4.0,  4)
    np_bins_y = np.append( np.append(np_bins_y_p0, np_bins_y_p1), np_bins_y_p2)

    f = ROOT.TFile('../root/tree_histos2_CC.root')

    counter = 0
    for iqt in range(np_bins_qt.size-1):
        for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
            counter += 1
    res_mean      = np.zeros(counter)
    res_mean_err  = np.zeros(counter)
    res_sigma     = np.zeros(counter)
    res_sigma_err = np.zeros(counter)
    p_values = np.zeros(counter)
    print counter, 'histograms'

    counter = 0
    for iqt in range(np_bins_qt.size-1):
        for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
            bin_y = 'y{:03.2f}'.format(np_bins_y[iy])+'_'+'y{:03.2f}'.format(np_bins_y[iy+1])
            bin_qt = 'qt{:03.1f}'.format(np_bins_qt[iqt])+'_'+'qt{:03.1f}'.format(np_bins_qt[iqt+1])
            name = charge+'_'+var+'_'+coeff_eval+'_'+bin_y+'_'+bin_qt
            canvas = ROOT.TCanvas("c", "canvas", 600, 600)            
            h = f.Get(charge+'/'+var+'/'+coeff_eval+'/'+name)
            h_norm = f.Get(charge+'/'+var+'/'+coeff_eval+'/'+name+'_norm')
            norm2 = 0.0
            err2  = 0.0
            ndof = 0
            for ix in range(1, h.GetNbinsX()+1):
                for iy in range(1, h.GetNbinsY()+1):
                    if (h_norm.GetBinContent(ix,iy) < min_val or h_norm.GetBinContent(ix,iy)==0.):
                        continue
                    ndof += 1
                    err = h.GetBinError(ix,iy)
                    norm2 +=  h.GetBinContent(ix,iy)/err/err
                    err2  +=  1./err/err
            norm = norm2/err2 if err2>0.0 else 0.0

            chi2 = 0.0
            h_pull = ROOT.TH1F('h_pull_byInvPDF'+name, '', 81, -4.0, 4.0)                        
            for ix in range(1, h.GetNbinsX()+1):
                for iy in range(1, h.GetNbinsY()+1):
                    val = h.GetBinContent(ix,iy)
                    err = h.GetBinError(ix,iy)
                    if (h_norm.GetBinContent(ix,iy) < min_val or h_norm.GetBinContent(ix,iy)==0.):
                        continue
                    pull = (val - norm)/err
                    chi2 +=  pull*pull
                    h_pull.Fill(pull)
                    h.SetBinContent(ix,iy, h.GetBinContent(ix,iy)-norm )

            p_values[counter] =  (1-stats.chi2.cdf(chi2, ndof-1 ) - 0.5)
            res_mean[counter] = h_pull.GetMean()
            res_mean_err[counter] = h_pull.GetMeanError()
            res_sigma[counter] = h_pull.GetRMS()-1.0
            res_sigma_err[counter] = h_pull.GetRMSError()
            print name, 'mean[', counter , ']=', res_mean[counter] 

            if save_2D:
                h_pull.Draw("HIST")
                canvas.SaveAs('plots/pull_1D_byInvPDF_'+name+'.png')
                h.Draw("COLZ") 
                canvas.SaveAs('plots/pull_2Dflat_byInvPDF_COLZ_'+name+'.png')
                h.Draw("SURF") 
                canvas.SaveAs('plots/pull_2Dflat_byInvPDF_SURF_'+name+'.png')
            canvas.IsA().Destructor( canvas )
            h_pull.IsA().Destructor( h_pull ) 
            counter += 1
            
    if save_summary:
        plt.figure()
        fig, ax = plt.subplots()
        colors = ['b', 'r', 'g']
        fmts = ['o', '^', '+']
        x = np.linspace(0,counter-1,counter)
        ax.errorbar(x, res_mean, xerr=0.0, yerr=0.0, fmt=fmts[0], color=colors[0], label='mean', markersize=6)
        ax.errorbar(x, res_sigma, xerr=0.0, yerr=0.0, fmt=fmts[1], color=colors[1], label='RMS-1', markersize=6)
        ax.errorbar(x, p_values, xerr=0.0, yerr=0.0, fmt=fmts[2], color=colors[2], label='p-0.5', markersize=8)
        plt.axis([-1, counter, -0.51, +0.51])
        legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
        plt.xlabel('Bin number', fontsize=20)
        plt.ylabel('Value', fontsize=20)
        plt.title('['+charge+', '+var+'], type = '+coeff_eval+'. Flattened by 1/pdf', fontsize=20)
        plt.grid(True)
        plt.show()
        plt.savefig('plots/summary_byInvPDF_'+charge+'_'+var+'_'+coeff_eval+'.png')
        plt.close()

    print 'Done'


def plot_closure_test_cum(charge='Wplus', var='Wdress', coeff_eval='val', min_val=-999., verbose=False, save_2D=True, save_summary=True):

    print "*** plot_closure_test_cum ***"

    from pylab import rcParams
    rcParams['figure.figsize'] = 14,8

    from tree_utils import np_bins_qt, np_bins_y, weight_coeff_cumulative
    #np_bins_qt = np.array([0.0, 2.0])
    #np_bins_y  = np.array([-0.2, 0.0, 0.2])
    np_bins_y_p0 = np.linspace(-4.0, -2.5,  4)
    np_bins_y_p1 = np.linspace(-2.0, +2.0, 21)
    np_bins_y_p2 = np.linspace(+2.5, +4.0,  4)
    np_bins_y = np.append( np.append(np_bins_y_p0, np_bins_y_p1), np_bins_y_p2)

    f = ROOT.TFile('../root/tree_histos2_CC.root')
    #res = pickle.load( open(os.environ['CMSSW_BASE']+'/src/Wmass/data/'+'fit_results_CC_'+charge+'_all.pkl') )
    res = pickle.load( open(os.environ['CMSSW_BASE']+'/src/Wmass/test/plots/'+'fit_results_CC_'+charge+'_all_A0-7.pkl') )

    counter = 0
    for iqt in range(np_bins_qt.size-1):
        for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
            counter += 1
    res_mean      = np.zeros(counter)
    res_mean_err  = np.zeros(counter)
    res_sigma     = np.zeros(counter)
    res_sigma_err = np.zeros(counter)
    p_values = np.zeros(counter)
    print counter, 'histograms'

    counter = 0
    for iqt in range(np_bins_qt.size-1):
        for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
            bin_y  = 'y{:03.2f}'.format(np_bins_y[iy])+'_'+'y{:03.2f}'.format(np_bins_y[iy+1])
            bin_qt = 'qt{:03.1f}'.format(np_bins_qt[iqt])+'_'+'qt{:03.1f}'.format(np_bins_qt[iqt+1])
            name = charge+'_'+var+'_'+coeff_eval+'_'+bin_y+'_'+bin_qt
            canvas = ROOT.TCanvas("c", "canvas", 600, 600)            
            h_pull = ROOT.TH1F('h_pull_byPDF_'+name, '', 81, -4.0, 4.0)
            h = f.Get(charge+'/'+var+'/'+coeff_eval+'/'+name+'_norm')
            chi2 = 0.0
            integ = h.Integral()
            ndof = 0
            for ix in range(1, h.GetNbinsX()+1):
                bin_cos = (h.GetXaxis().GetBinLowEdge(ix), 
                           h.GetXaxis().GetBinLowEdge(ix)+h.GetXaxis().GetBinWidth(ix))
                for iy in range(1, h.GetNbinsY()+1):
                    bin_phi = (h.GetYaxis().GetBinLowEdge(iy), 
                               h.GetYaxis().GetBinLowEdge(iy)+h.GetYaxis().GetBinWidth(iy))                    
                    norm = weight_coeff_cumulative(res=res, 
                                                   coeff_eval=coeff_eval, 
                                                   bin_y=bin_y, 
                                                   qt=0.5*( np_bins_qt[iqt]+ np_bins_qt[iqt+1]), 
                                                   bin_cos=bin_cos, bin_phi=bin_phi, 
                                                   coeff=['A0', 'A1', 'A2', 'A3', 'A4'])
                    norm *= integ
                    if verbose:
                        print '\tBin (', ix, ',', iy , ') = ' , h.GetBinContent(ix,iy), ' +/- ', h.GetBinError(ix,iy)
                    if (h.GetBinContent(ix,iy) < min_val or h.GetBinContent(ix,iy)==0.):
                        continue
                    pull = (h.GetBinContent(ix,iy) - norm)/h.GetBinError(ix,iy)
                    h_pull.Fill( pull )
                    chi2 += pull*pull
                    ndof += 1
                    h.SetBinContent(ix,iy, h.GetBinContent(ix,iy) - norm )

            p_values[counter] =  (1-stats.chi2.cdf(chi2, ( ndof - 1) ) - 0.5)
            res_mean[counter] = h_pull.GetMean()
            print name, 'mean[', counter , ']=', res_mean[counter] 
            res_mean_err[counter] = h_pull.GetMeanError()
            res_sigma[counter] = h_pull.GetRMS()-1.0
            res_sigma_err[counter] = h_pull.GetRMSError()
            counter += 1
            if save_2D:
                h_pull.Draw("HIST")
                canvas.SaveAs('plots/pull_1D_byPDF_'+name+'.png')
                h.SetStats(0)                    
                h.Draw("COLZ") 
                canvas.SaveAs('plots/pull_2Dflat_byPDF_COLZ_'+name+'.png')
                h.Draw("SURF") 
                canvas.SaveAs('plots/pull_2Dflat_byPDF_SURF_'+name+'.png')
            h_pull.IsA().Destructor( h_pull )
            canvas.IsA().Destructor( canvas )

    if save_summary:
        plt.figure()
        fig, ax = plt.subplots()
        colors = ['b', 'r', 'g']
        fmts = ['o', '^', '+']
        x = np.linspace(0,counter-1,counter)
        ax.errorbar(x, res_mean, xerr=0.0, yerr=0.0, fmt=fmts[0], color=colors[0], label='mean', markersize=6)
        ax.errorbar(x, res_sigma, xerr=0.0, yerr=0.0, fmt=fmts[1], color=colors[1], label='RMS-1', markersize=6)
        ax.errorbar(x, p_values, xerr=0.0, yerr=0.0, fmt=fmts[2], color=colors[2], label='p-0.5', markersize=8)
        plt.axis([-1, counter, -0.51, +0.51])
        legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
        plt.xlabel('Bin number', fontsize=20)
        plt.ylabel('Value', fontsize=20)
        plt.title('['+charge+', '+var+'], type = '+coeff_eval+'. $\chi^2$ from pdf', fontsize=20)
        plt.grid(True)
        plt.show()
        plt.savefig('plots/summary_byPDF_'+charge+'_'+var+'_'+coeff_eval+'.png')
        plt.close()

    print 'Done'
