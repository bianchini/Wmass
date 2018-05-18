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

# closure test of the angular coefficients
def plot_closure_test(charge='Wplus', DY='CC', var='Wdress', coeff_eval='val', 
                      min_val=-999., 
                      coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7'],
                      byInvPdf=False,
                      verbose=False, 
                      save_2D=True, save_pdf=True, save_summary=True, 
                      do_toy=False, extra_variance_toy=0.07,
                      do_pt_eta=False):
    
    print "Running plot_closure_test"
    postfix = ''
    if byInvPdf:
        postfix = 'byInvPDF'
        print '> Compare MC distribution in (cos,phi) weighted by 1/pdf with flat distribution'
    else:
        postfix = 'byPDF'
        print '> Compare MC distribution in (cos,phi) with N*pdf'        

    if do_toy:
        print 'Set rnd seed to 0'
        np.random.seed(0)

    from pylab import rcParams
    rcParams['figure.figsize'] = 14,8

    from tree_utils import np_bins_qt, np_bins_y, weight_coeff_cumulative, get_coeff_vals
    np_bins_y_p0 = np.linspace(-4.0, -2.5,  4)
    np_bins_y_p1 = np.linspace(-2.0, +2.0, 21)
    np_bins_y_p2 = np.linspace(+2.5, +4.0,  4)
    np_bins_y = np.append( np.append(np_bins_y_p0, np_bins_y_p1), np_bins_y_p2)

    # for fast check:
    np_bins_qt = np.array([0.0, 2.0])    
    #np_bins_y  = np.array([0.8, 1.0])  
    print 'The following (qt,y) bins will be considered:'
    print '> qt:', np_bins_qt
    print '>  y:', np_bins_y

    fin_name = '../root/tree_histos'+('2' if not do_pt_eta else '3')+'_'+DY+'.root'
    fin = ROOT.TFile(fin_name, 'READ')
    print 'Read histograms from '+fin_name

    res_name = os.environ['CMSSW_BASE']+'/src/Wmass/data/'+'fit_results_'+DY+'_'+charge+'_'+var+'_all_A0-7.pkl'
    res = pickle.load( open(res_name) )
    print 'Read coefficients from '+res_name

    counter = 0
    for iqt in range(np_bins_qt.size-1):
        for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
            counter += 1
    res_mean      = np.zeros(counter)
    res_mean_err  = np.zeros(counter)
    res_sigma     = np.zeros(counter)
    res_sigma_err = np.zeros(counter)
    p_values = np.zeros(counter)
    print 'Total number of histograms:', counter

    x_labels = []
    counter = 0
    for iqt in range(np_bins_qt.size-1):
        bin_qt = 'qt{:03.1f}'.format(np_bins_qt[iqt])+'_'+'qt{:03.1f}'.format(np_bins_qt[iqt+1])
        x_labels.append( '['+'{:0.0f}'.format(np_bins_qt[iqt])+','+'{:0.0f}'.format(np_bins_qt[iqt+1])+']' )
        for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
            bin_y  = 'y{:03.2f}'.format(np_bins_y[iy])+'_'+'y{:03.2f}'.format(np_bins_y[iy+1])
            name = charge+'_'+var+'_'+coeff_eval+'_'+bin_y+'_'+bin_qt
            canvas = ROOT.TCanvas("c", "canvas", 600, 600)            
            h_pull = ROOT.TH1F('h_pull_byPDF_'+name, '', 81, -4.0, 4.0)

            (h, h_norm) = (None, None)
            if not do_pt_eta:
                (h, h_norm) = (fin.Get(charge+'/'+var+'/'+coeff_eval+'/'+name),
                               fin.Get(charge+'/'+var+'/'+coeff_eval+'/'+name+'_norm'))
            else:
                coeff_vals = get_coeff_vals(res=res, coeff_eval=coeff_eval, bin_y=bin_y, qt=0.5*(np_bins_qt[iqt]+np_bins_qt[iqt+1]), coeff=coeff)
                h = fin.Get(charge+'/'+var+'/'+coeff_eval+'/'+'M'+'{:05.3f}'.format(80.419)+'/'+bin_qt+'_'+bin_y+'/'+charge+'_'+var+'_'+coeff_eval+'_'+'M'+'{:05.3f}'.format(80.419)+'_'+bin_qt+'_'+bin_y+'_')
                h_norm_UL = fin.Get(charge+'/'+var+'/'+coeff_eval+'/'+'M'+'{:05.3f}'.format(80.419)+'/'+bin_qt+'_'+bin_y+'/'+charge+'_'+var+'_'+coeff_eval+'_'+'M'+'{:05.3f}'.format(80.419)+'_'+bin_qt+'_'+bin_y+'_UL')
                scale_UL = 1.0
                for ic,c in enumerate(coeff):
                    scale_UL -= coeff_vals[ic]
                    (dir_name, h_name) = (charge+'/'+var+'/'+coeff_eval+'/'+'M'+'{:05.3f}'.format(80.419)+'/'+bin_qt+'_'+bin_y, 
                                          charge+'_'+var+'_'+coeff_eval+'_'+'M'+'{:05.3f}'.format(80.419)+'_'+bin_qt+'_'+bin_y+'_'+c)
                    if h_norm==None:
                        h_norm = fin.Get(dir_name+'/'+h_name).Clone(h_name)
                        h_norm.Scale( coeff_vals[ic] )
                    else:
                        h_norm.Add( fin.Get(dir_name+'/'+h_name), coeff_vals[ic] )
                h_norm.Add(h_norm_UL, scale_UL)

            h_pdf = h_norm.Clone(h_norm.GetName()+'_pdf')
            h_pdf.Reset()

            if do_toy:
                from tree_utils import make_grid_CS
                ntoys=int(h_norm.Integral())
                print 'Generating ', ntoys, 'toys'
                h_norm.Reset()
                make_grid_CS(res=res, h=h_norm, coeff_eval=coeff_eval, 
                             bin_y=bin_y, qt=0.5*( np_bins_qt[iqt]+ np_bins_qt[iqt+1] ), 
                             coeff=coeff, ntoys=ntoys)
                             
            integ = h_norm.Integral() if not byInvPdf else h.Integral()
            chi2 = 0.0
            ndof = 0
            for ix in range(1, h_norm.GetNbinsX()+1):
                bin_cos = (h_norm.GetXaxis().GetBinLowEdge(ix), 
                           h_norm.GetXaxis().GetBinLowEdge(ix)+h_norm.GetXaxis().GetBinWidth(ix))
                for iy in range(1, h_norm.GetNbinsY()+1):
                    bin_phi = (h_norm.GetYaxis().GetBinLowEdge(iy), 
                               h_norm.GetYaxis().GetBinLowEdge(iy)+h_norm.GetYaxis().GetBinWidth(iy))                    

                    n = h_norm.GetBinContent(ix,iy) if not byInvPdf else h.GetBinContent(ix,iy)
                    err = h_norm.GetBinError(ix,iy) if not byInvPdf else h.GetBinError(ix,iy)

                    if do_toy:
                        err_weight = -1.0
                        while err_weight<0.:
                            err_weight = err*(1.0 + np.random.normal(0.0, extra_variance_toy))
                        err = err_weight

                    mu = 0.0
                    if byInvPdf:
                        mu = 1.0/(h.GetNbinsX()*h.GetNbinsY())
                    else:
                        if not do_pt_eta:
                            mu = weight_coeff_cumulative(res=res, 
                                                         coeff_eval=coeff_eval, 
                                                         bin_y=bin_y, 
                                                         qt=0.5*( np_bins_qt[iqt]+ np_bins_qt[iqt+1]), 
                                                         bin_cos=bin_cos, bin_phi=bin_phi, 
                                                         coeff=coeff)
                        else:
                            mu = h.GetBinContent(ix,iy)/integ

                    if verbose:
                        print '\tBin (', ix, ',', iy , ') = ' , n, ' +/- ', err
                    if (n < min_val or n==0.):
                        continue

                    pull = (n - mu*integ)/err
                    h_pull.Fill( pull )
                    chi2 += pull*pull
                    ndof += 1
                    h_norm.SetBinContent(ix,iy, pull )
                    h_pdf.SetBinContent(ix,iy, mu )

            p_values[counter] =  (1-stats.chi2.cdf(chi2, ( ndof - 1) ) - 0.5)
            res_mean[counter] = h_pull.GetMean()
            print name, 'mean[', counter , ']=', res_mean[counter] 
            res_mean_err[counter]  = h_pull.GetMeanError()
            res_sigma[counter]     = h_pull.GetRMS()-1.0
            res_sigma_err[counter] = h_pull.GetRMSError()
            counter += 1
            if save_2D:
                h_pull.Draw("HIST")
                canvas.SaveAs('plots/pull_1D_'+postfix+'_'+DY+'_'+name+('_toy' if do_toy else '')+'.png')                
                h_norm.SetMaximum(+4.0)
                h_norm.SetMinimum(-4.0)
                h_norm.SetStats(0)                    
                h_norm.Draw("COLZ") 
                canvas.SaveAs('plots/pull_2Dflat_'+postfix+'_COLZ_'+DY+'_'+name+('_toy' if do_toy else '')+'.png')
                h_norm.Draw("LEGO") 
                canvas.SaveAs('plots/pull_2Dflat_'+postfix+'_LEGO_'+DY+'_'+name+('_toy' if do_toy else '')+'.png')
            if save_pdf:
                h_pdf.Draw("COLZ") 
                canvas.SaveAs('plots/pdf_2Dflat_'+postfix+'_COLZ_'+DY+'_'+name+('_toy' if do_toy else '')+'.png')
                h_pdf.Draw("SURF") 
                canvas.SaveAs('plots/pdf_2Dflat_'+postfix+'_SURF_'+DY+'_'+name+('_toy' if do_toy else '')+'.png')
            h_pull.IsA().Destructor( h_pull )
            canvas.IsA().Destructor( canvas )

    if save_summary:
        print 'Summary plot...'
        plt.figure()
        fig, ax = plt.subplots()
        colors = ['b', 'r', 'g']
        fmts = ['o', '^', '+']
        x = np.linspace(0,counter-1,counter)
        ax.errorbar(x, res_mean, xerr=0.0, yerr=0.0, fmt=fmts[0], color=colors[0], label='mean', markersize=6)
        ax.errorbar(x, res_sigma, xerr=0.0, yerr=res_sigma_err, fmt=fmts[1], color=colors[1], label='RMS-1', markersize=6)
        ax.errorbar(x, p_values, xerr=0.0, yerr=0.0, fmt=fmts[2], color=colors[2], label='p-0.5', markersize=8)
        plt.axis([-1, counter, -0.51, +0.51])        
        ax.set_xticks( x[ np.arange(0, x.size-1, np_bins_y.size/2) ] )
        ax.set_xticklabels( x_labels, rotation=45, ha='center', fontsize=15 )
        legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
        plt.ylabel('Value', fontsize=20)
        plt.title('['+charge+', '+var+'], type = '+coeff_eval+('. Flattened by 1/pdf' if byInvPdf else '. Compared to N*pdf.')+' Data = '+('toys' if do_toy else 'MC'), fontsize=20)
        plt.grid(True)
        plt.show()
        plt.savefig('plots/summary_'+postfix+'_'+DY+'_'+charge+'_'+var+'_'+coeff_eval+('_toy' if do_toy else '')+'.png')
        plt.close()

    fin.Close()
    print 'Done!'

def merge_templates(charges=['Wplus'], var=['WpreFSR'], coeff_eval=['val'], masses=[80.419], coeff=['A0'], 
                    np_bins_template_qt=np.array([]), np_bins_template_y=np.array([]), rebin=(),
                    postfix=''):

    from tree_utils import np_bins_y, np_bins_qt, np_bins_eta, np_bins_pt, angular_pdf_string

    # available templates
    np_bins_y_extL = np.insert(np_bins_y,   0, [-10.])
    np_bins_y_ext  = np.append(np_bins_y_extL, [+10.])
    np_bins_qt_ext = np.append(np_bins_qt, 999.)
    nbins_y  = np_bins_y_ext.size - 1 
    nbins_qt = np_bins_qt_ext.size - 1 

    # target templates
    np_bins_template_y_extL = np.insert(np_bins_template_y,   0, [-10.])
    np_bins_template_y_ext  = np.append(np_bins_template_y_extL, [+10.])
    np_bins_template_qt_ext = np.append(np_bins_template_qt, 999.)
    nbins_template_y  = np_bins_template_y_ext.size - 1 
    nbins_template_qt = np_bins_template_qt_ext.size - 1 

    (np_bins_rebin_eta, np_bins_rebin_pt) = (np_bins_eta,np_bins_pt)    
    if rebin!=():
        np_bins_rebin_eta = np.linspace(np_bins_eta[0], np_bins_eta[-1], (np_bins_eta.size-1)/rebin[0]+1 )
        np_bins_rebin_pt  = np.linspace(np_bins_pt[0],  np_bins_pt[-1],  (np_bins_pt.size-1)/rebin[1]+1 )

    xx, yy = np.meshgrid(np_bins_rebin_pt, np_bins_rebin_eta)        
    coeff_ext = coeff+['UL', '']
    np_coeff_ext = np.chararray(len(coeff_ext), itemsize=2, unicode=True)
    np_coeff_ext[:] = coeff_ext

    fin = ROOT.TFile('../root/tree_histos3_CC_FxFx_masses.root', 'READ')

    # create TH2D
    for q in charges:
        for v in var:
            for ceval in coeff_eval:
                template = np.zeros( (len(masses), nbins_template_qt, nbins_template_y/2, len(coeff_ext), np_bins_rebin_eta.size-1, np_bins_rebin_pt.size-1) )    
                mc_acceptances = np.zeros( (len(masses), nbins_template_qt, nbins_template_y/2) )
                bin_template_qt_index = -1
                for qt in range(1, nbins_template_qt+1):
                    bin_template_qt_index += 1
                    qt_template_bin = 'qt{:03.1f}'.format(np_bins_template_qt_ext[qt-1])+'_'+'qt{:03.1f}'.format(np_bins_template_qt_ext[qt]) if qt<nbins_template_qt else 'OF'
                    (qt_template_low, qt_template_high) = (np_bins_template_qt_ext[qt-1], np_bins_template_qt_ext[qt])

                    bin_template_y_index = -1
                    for y in range(nbins_template_y/2+1, nbins_template_y+1):
                        bin_template_y_index += 1
                        y_template_bin = 'y{:03.2f}'.format(np_bins_template_y_ext[y-1])+'_'+'y{:03.2f}'.format(np_bins_template_y_ext[y]) if y<nbins_template_y else 'OF'
                        (y_template_low, y_template_high) = (np_bins_template_y_ext[y-1], np_bins_template_y_ext[y])
                                                  
                        template_bin_name = qt_template_bin+'_'+y_template_bin
                        print 'Doing....'+template_bin_name

                        iqts = []
                        for iqt in range(1, nbins_qt+1):
                            (qt_low, qt_high) = (np_bins_qt_ext[iqt-1], np_bins_qt_ext[iqt])
                            if (qt_low>qt_template_low or np.isclose(qt_low,qt_template_low)) and (qt_high<qt_template_high or np.isclose(qt_high,qt_template_high)):
                                print '\t\t[',qt_low,',',qt_high, '] < [', qt_template_low,',',qt_template_high, ']'
                                iqts.append(iqt)

                        iys = []
                        for iy in range(nbins_y/2+1, nbins_y+1):
                            (y_low, y_high) = (np_bins_y_ext[iy-1], np_bins_y_ext[iy])
                            if (y_low>y_template_low or np.isclose(y_low,y_template_low)) and (y_high<y_template_high or np.isclose(y_high,y_template_high)):
                                print '\t\t[',y_low,',',y_high, '] < [', y_template_low,',',y_template_high, ']'
                                iys.append(iy)

                        print '\tMerging qt bins with index:', [ iqt-1 for iqt in iqts ]
                        print '\tMerging y  bins with index:', [ iy-1  for iy  in iys ]

                        for im,m in enumerate(masses):
                            mass_str = 'M'+'{:05.3f}'.format(m)
                            for ic,c in enumerate(coeff_ext):                                
                                norm_c = 0.0 
                                print '\t\tCoefficient...'+c
                                template_name = q+'_'+v+'_'+ceval+'_'+mass_str+'_'+qt_template_bin+'_'+y_template_bin+'_'+c
                                h = None
                                for iqt in iqts:
                                    qt_bin = 'qt{:03.1f}'.format(np_bins_qt_ext[iqt-1])+'_'+'qt{:03.1f}'.format(np_bins_qt_ext[iqt]) if iqt<nbins_qt else 'OF'
                                    for iy in iys:
                                        y_bin = 'y{:03.2f}'.format(np_bins_y_ext[iy-1])+'_'+'y{:03.2f}'.format(np_bins_y_ext[iy]) if iy<nbins_y else 'OF'
                                        (dir_name, h_name) = (q+'/'+v+'/'+ceval+'/'+mass_str+'/'+qt_bin+'_'+y_bin, q+'_'+v+'_'+ceval+'_'+mass_str+'_'+qt_bin+'_'+y_bin+'_'+c)
                                        h_tmp = fin.Get(dir_name+'/'+h_name)
                                        norm_c += h_tmp.Integral(0, h_tmp.GetNbinsX()+1,0, h_tmp.GetNbinsY()+1)
                                        if h==None:
                                            h = h_tmp.Clone(template_name)
                                        else:
                                            h.Add( h_tmp )                                        

                                if rebin!=():
                                    h.Rebin2D(rebin[0], rebin[1])

                                # templates are normalised in acceptance
                                for ipt in range(np_bins_rebin_pt.size-1):
                                    for ieta in range(np_bins_rebin_eta.size-1):
                                        template[im][bin_template_qt_index][bin_template_y_index][ic][ieta][ipt] = h.GetBinContent(ieta+1,ipt+1)/(norm_c if c!='' else 1.0)
                                        
                                # save a map of MC acceptances (for closure-test)
                                if c=='':
                                    mc_acceptances[im][bin_template_qt_index][bin_template_y_index] = h.Integral()/norm_c 

                                plt.pcolormesh(yy, xx, template[im][bin_template_qt_index][bin_template_y_index][ic])
                                plt.colorbar()

                                words = template_name.split('_')
                                title = r'$'+words[0][0]+'^{'+('+' if 'plus' in words[0] else '-')+'}$, '
                                title += r'$q_{T}\in['+'{:0.0f}'.format(np_bins_template_qt_ext[qt-1])+', '+('{:0.0f}'.format(np_bins_template_qt_ext[qt]) if qt<nbins_template_qt else '\infty')+']$ GeV'+', '
                                title += r'$|y|\in['+'{:0.1f}'.format(np_bins_template_y_ext[y-1])+', '+('{:0.1f}'.format(np_bins_template_y_ext[y]) if y<nbins_template_y else '\infty')+']$'
                                plt.title(title, fontsize=20)
                                coeff_vals = np.zeros(8)
                                if c not in ['UL', '']:
                                    coeff_vals[int(c[1])] = 1.0
                                pdf_c = angular_pdf_string(coeff_vals=coeff_vals) if c!='' else 'MC'
                                plt.axis([np_bins_rebin_eta[0], np_bins_rebin_eta[-1], np_bins_rebin_pt[0], np_bins_rebin_pt[-1]])        
                                plt.figtext(0.15, 0.86, r'$M_{W}$ = '+words[3][1:]+' GeV', color='white')
                                acc = template[im][bin_template_qt_index][bin_template_y_index][ic].sum() / (norm_c if c=='' else 1.0)
                                acc_err = math.sqrt(acc*(1-acc)/norm_c)
                                plt.figtext(0.15, 0.81, r'$\frac{d\sigma}{dq_{T}d|y|}$ = '+'{:0.1f}'.format(norm_c)+', '+
                                            '$\epsilon_{A}$ = '+'{:0.4f}'.format(acc)+
                                            '$\pm$'+'{:0.4f}'.format(acc_err), color='white')
                                plt.figtext(0.15, 0.76, r'$\frac{1}{\sigma}\frac{d\sigma}{d\Omega}$ = '+pdf_c, color='white')
                                plt.xlabel('$\eta$', fontsize=20)
                                plt.ylabel('$p_{T}$ (GeV)', fontsize=20)
                                plt.show()
                                plt.savefig('plots/template_'+template_name+'.png')
                                plt.close('all')
                                #return

                outname = 'plots/template_'+q+'_'+v+'_'+ceval+postfix
                np.savez(outname,
                         template, 
                         np.array(masses), 
                         np_bins_template_qt_ext, 
                         np_bins_template_y_ext[(nbins_template_y/2):], 
                         np_coeff_ext, 
                         np_bins_rebin_eta, np_bins_rebin_pt,
                         mc_acceptances)

                # for validation
                print 'Content of file '+outname+'.npz'
                res = np.load(outname+'.npz')
                for f in res.files:
                    print 'File: '+f, 'with size = ', res[f].size
                    print res[f]
                    
# draw derivative of template against one coefficient at the time
def derivative_templates(charge='Wplus', var='WpreFSR', coeff_eval='val', bin=(0,0)):

    rcParams['figure.figsize'] = 8,8

    outname = 'plots/template_'+charge+'_'+var+'_'+coeff_eval
    print 'Content of file '+outname+'.npz'
    res = np.load(outname+'.npz')
    for f in res.files:
        print 'File: '+f, 'with size = ', res[f].shape

    (iqt,iy) = bin

    templates = []
    for i in range(8):
        templates.append( (res['arr_0'][0][iqt][iy][i]-res['arr_0'][0][iqt][iy][-2]) )

    np_bins_qt = res['arr_2']
    np_bins_y  = res['arr_3']
    np_bins_rebin_eta = res['arr_5']
    np_bins_rebin_pt  = res['arr_6']
    xx,yy = np.meshgrid(np_bins_rebin_pt, np_bins_rebin_eta)        
    
    for i in range(8):
        plt.subplot(4, 2, i+1)
        plt.pcolormesh(yy, xx, templates[i])
        plt.colorbar(format='%.0e')
        plt.title(r'$A_{'+str(i)+r'}$')
        plt.axis([np_bins_rebin_eta.min(), np_bins_rebin_eta.max(), np_bins_rebin_pt.min(), np_bins_rebin_pt.max()])        
        if i%2==0:
            plt.ylabel('$p_{T}$ (GeV)', fontsize=12)
        if i>=6:
            plt.xlabel('$\eta$', fontsize=12)

    plt.suptitle(r'$q_{T} \in ['+'{:0.1f}'.format(np_bins_qt[iqt])+', '+'{:0.1f}'.format(np_bins_qt[iqt+1])+']$ GeV, '+'$|y| \in ['+'{:0.1f}'.format(np_bins_y[iy])+', '+'{:0.1f}'.format(np_bins_y[iy+1])+']$', fontsize=20)
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.show()
    plt.savefig('plots/derivative_'+charge+'_'+var+'_'+coeff_eval+'.png')
    plt.close('all')



