import math
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 8,7

from scipy import stats
from scipy.stats import f
from array import array
import ROOT


def polynomial(x, coeff, order):
    val = np.zeros(x.size)
    for p in range(order+1):
        val +=  coeff[p]*np.power(x,p)
    return val

# Perform Fisher-test based on chi2 of fit to pol_order: [0] + [1]*x + [2]*x*x + ... 
# Stop when p-value of chi2>0.05
# If force_order>=0, force order; fix_to_zero requires [0] := 0.0
def Fisher_test(h=None, coeff='A0', fix_to_zero={}, fit_range=[0.0,50.0], force_order=-1, start_order=0, threshold=0.05, threshold_chi2=0.03):

    # a large initial value
    chi2 = 99999.

    # start from pol1
    order = start_order if force_order<0 else force_order

    # initial values
    r = None
    fit = None
    pval = 0.0
    chi2_prob = 0.0

    while( pval<threshold ):

        formula = ''
        for p in range(order+1):
            formula += ('['+str(p)+']*TMath::Power(x,'+str(p)+')'+('+' if p<order else ''))
        #print formula
        this_fit = ROOT.TF1('fit_order'+str(order), formula,  fit_range[0], fit_range[1]) 
        this_fit.SetNpx(10000)
        for p in range(order+1):
            if p in fix_to_zero[coeff]:
                this_fit.FixParameter(p, 0.0 )
                continue
            this_fit.SetParameter(p, math.pow(10,-p))

        this_r = h.Fit('fit_order'+str(order), 'SRQ')
        ndof = this_r.Ndf()

        # protection against impossible fits
        if ndof<=1:
            break

        this_chi2_prob = 1-stats.chi2.cdf( this_r.Chi2(), this_r.Ndf() )

        if force_order>=0:
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()
            chi2_prob = this_chi2_prob
            break

        if order==start_order:
            print '\tOrder', order, 'chi2=', this_r.Chi2(), '/', this_r.Ndf(), ' [p =', this_chi2_prob, ']: go to next order' 
            order += 1 
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()    
            chi2_prob = this_chi2_prob
            continue

        if this_chi2_prob<threshold_chi2:
            print '\tOrder', order, 'chi2=', this_r.Chi2(), '/', ndof,  ' [p =', this_chi2_prob, ']: this fit is bad, continue'
            order += 1         
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()        
            chi2_prob = this_chi2_prob            
            continue

        if chi2_prob<threshold_chi2:
            print '\tOrder', order, 'chi2=', this_r.Chi2(), '/', ndof,  ' [p =', this_chi2_prob, ']: this fit is ok, previous is bad, continue'
            order += 1         
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()        
            chi2_prob = this_chi2_prob            
            continue

        pval = f.cdf(this_r.Chi2()/chi2*(ndof)/(ndof-1) , ndof, ndof-1)
        if pval<threshold:
            print '\tOrder', order, 'chi2=', this_r.Chi2(), '/', ndof,  \
                ' [p =', '{:0.3f}'.format(this_chi2_prob), ']',  ': pval=', '{:0.3f}'.format(pval), \
                ' small, continue' 
            order += 1         
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()        
            chi2_prob = this_chi2_prob            
        else:
            print '\tOrder', order, 'chi2=', '{:0.2f}'.format(this_r.Chi2()), '/', ndof, \
                ' [p =', '{:0.3f}'.format(this_chi2_prob), ']', ': pval=', '{:0.3f}'.format(pval), '>=', threshold, \
                ', stop!'        

    if force_order<0:
        # go back by one order: the match was at the order before!
        order -= 1
        print '\tP-value threshold at order', order, \
            'with chi2=', '{:0.2f}'.format(r.Chi2()), '/', r.Ndf()
    else:
        print '\tFit performed at order ', force_order, \
            'with chi2=', '{:0.2f}'.format(r.Chi2()), '/', r.Ndf() 

    # nice formatting of results
    res = ''
    for p in range(order+1):
        exp = '{:.1E}'.format(r.Parameter(p))[-3]+'{:.1E}'.format(r.Parameter(p))[-1] 
        val = '{:.1E}'.format(r.Parameter(p))[:-4]
        err = '{:.1E}'.format(r.ParError(p))[:-4]
        res += 'c_{'+str(p)+'}=('+val+' #pm '+err+')#times10^{'+exp+'}, '

    # remove 1 from order to go to the previous one
    return (r,order,fit,res,chi2_prob)

# plot covariance matrix provided by cov as a .png, .npy, and .C
def plot_cov_matrix(n_vars=0, cov=None, label='', plot='cov', save_all=False):
    c = ROOT.TCanvas("canvas", "canvas", 800, 800) 
    h2 = ROOT.TH2D('cov', '', n_vars, 0, n_vars, n_vars, 0, n_vars)
    h2.SetStats(0) 
    for i in range(n_vars):
        for j in range(n_vars):            
            rho_ij = 0.0
            if plot=='corr':
                rho_ij = cov[i][j]/math.sqrt(cov[i,i]*cov[j,j]) if cov[i][i]>0.0 and cov[j][j]>0.0 else 0.0
            else:
                rho_ij = cov[i][j]
            h2.SetBinContent(i+1, j+1, rho_ij )
    h2.Draw("COLZ")
    c.SaveAs('plots/'+label+'_'+plot+'.png')
    if save_all:
        c.SaveAs('plots/'+label+'_'+plot+'.C')
        np.save('plots/'+label+'_'+plot, cov)
    c.IsA().Destructor( c )


# axis range when plotting the coefficients
def ranges_for_coeff():
    ranges = {}
    ranges['A0'] = (-0.2, 1.5)
    ranges['A1'] = (-0.8, 0.8)
    ranges['A2'] = (-0.2, 1.5)
    ranges['A3'] = (-2.0, 2.0)
    ranges['A4'] = (-2.5, 2.5)
    ranges['A5'] = (-0.5, 0.5)
    ranges['A6'] = (-0.5, 0.5)
    ranges['A7'] = (-0.5, 0.5)
    return ranges

def ranges_for_coeff_zoom(q='Wplus'):
    ranges = {}
    ranges['A0'] = [-0.1, 0.7]
    ranges['A1'] = [-0.1, 0.6]
    ranges['A2'] = [-0.1, 0.8]
    ranges['A3'] = [-1.5, 0.1] if q=='Wplus' else [-0.1, +1.5]
    ranges['A4'] = [-2.5, 0.1] if q=='Wplus' else [-0.1, +2.5]
    ranges['A5'] = [-0.5, 0.5]
    ranges['A6'] = [-0.5, 0.5]
    ranges['A7'] = [-0.5, 0.5]
    return ranges


# project the (y,qt) plot along the qt axis for all bins of y
def draw_y_slice(fname='./tree.root', var='Wdress', coeff='A0', weight_name=0, do_fit=True, threshold_chi2=0.02):

    fin = ROOT.TFile.Open(fname,'READ')

    ranges = ranges_for_coeff()

    histos = {}
    nbins_y = 0
    nbins_qt = 0
    for q in ['Wplus', 'Wminus']:
        histos[q] = fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name))
        histos[q+'_norm'] = fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name)+'_norm')
        nbins_y = histos[q].GetNbinsX()
        nbins_qt = histos[q].GetNbinsY()

    # merge bins corresponding to y and -y
    for y in range(nbins_y/2+1, nbins_y+1):
        c = ROOT.TCanvas('canvas_'+'y_'+str(y), "canvas", 800, 600) 
        leg = ROOT.TLegend(0.10,0.70,0.55,0.88, "","brNDC");
        leg.SetFillStyle(0);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.02);
        leg.SetFillColor(10);
        hslices = {} 
        for q in ['Wplus', 'Wminus']:
            h = histos[q] 
            h_norm = histos[q+'_norm']
            y_bin = 'y{:03.2f}'.format(h.GetXaxis().GetBinLowEdge(y))+'_'+'y{:03.2f}'.format(h.GetXaxis().GetBinLowEdge(y)+h.GetXaxis().GetBinWidth(y))
            print 'Bin Y: ', y_bin
            hslice_minus = h.ProjectionY(str(y)+'_'+q+'_minus_py', nbins_y+1-y, nbins_y+1-y)
            hslice_plus  = h.ProjectionY(str(y)+'_'+q+'_plus_py', y, y)
            hnorm_minus  = h_norm.ProjectionY(str(y)+'_'+q+'_minus_norm_py', nbins_y+1-y, nbins_y+1-y)
            hnorm_plus   = h_norm.ProjectionY(str(y)+'_'+q+'_plus_norm_py', y, y)
            print  '(', nbins_y+1-y, ' + ', y, ') /', nbins_y
            hslice_minus.Add(hslice_plus)
            hnorm_minus.Add(hnorm_plus)
            hslice_minus.Divide(hnorm_minus)
            hslice_minus.SetTitle(coeff[0]+'_{'+coeff[1]+'} for |y| #in ['+y_bin[1:]+']')
            hslice_minus.SetMinimum(ranges[coeff][0])
            hslice_minus.SetMaximum(ranges[coeff][1])
            hslice_minus.SetStats(0)
            hslice_minus.SetXTitle('q_{T} (GeV)')

            if do_fit:
                (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice_minus, coeff=coeff, 
                                                          fix_to_zero={'A0': [0], 'A1': [0], 'A2': [0], 'A3': [0], 'A4': [], 'A5':[0], 'A6': [0], 'A7': [0]}, 
                                                          fit_range=[0.0,50.0])
                if chi2_prob<threshold_chi2:
                    print 'Iterative fit did not converge: try with looser threshold on chi2_prob'
                    (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice_minus, coeff=coeff, 
                                                              fix_to_zero={'A0': [0], 'A1': [0], 'A2': [0], 'A3': [0], 'A4': [], 'A5': [0], 'A6': [0], 'A7': [0]},  
                                                              fit_range=[0.0,50.0], threshold_chi2=1e-04)
                    
                hslice_minus.GetFunction('fit_order'+str(order+1)).SetBit(ROOT.TF1.kNotDraw)

            hslice_minus.SetLineWidth(2)
            hslice_minus.SetLineColor(ROOT.kRed if q=='Wplus' else ROOT.kBlue)
            hslice_minus.SetFillColor(ROOT.kRed if q=='Wplus' else ROOT.kBlue)
            hslice_minus.SetFillStyle(3002 if q=='Wplus' else 3003)
            leg_text = 'W^{'+('+' if q=='Wplus' else '-')+'}'

            if do_fit:
                fit.SetLineWidth(2)
                fit.SetLineStyle(ROOT.kDashed)
                fit.SetLineColor(ROOT.kRed if q=='Wplus' else ROOT.kBlue)
                hslices[q+'_fit'] = fit
                leg_text = '#splitline{W^{'+('+' if q=='Wplus' else '-')+'}, #chi^{2}/ndof = '+'{:02.1f}'.format(r.Chi2()/r.Ndf())+' (prob='+'{:03.2f}'.format(chi2_prob)+ ')}{'+res+'}'
                leg.AddEntry(fit, 'Fit', "L")

            leg.AddEntry(hslice_minus, leg_text, "f")
            hslices[q] = hslice_minus

        hslices['Wplus'].Draw('E3')
        hslices['Wminus'].Draw('E3SAME')
        if do_fit:
            hslices['Wplus_fit'].Draw('SAME')
            hslices['Wminus_fit'].Draw('SAME')

        leg.Draw()
        c.SaveAs('plots/coefficient_'+var+'_'+coeff+'_'+str(weight_name)+'_'+y_bin+'.png')
        c.IsA().Destructor( c )        
        leg.IsA().Destructor( leg )        

    fin.Close()

# project the (y,qt) plot along the y axis for all bins of qt
def draw_qt_slice(fname='./tree.root', var='Wdress', coeff='A0', weight_name=0, do_fit=True):

    ranges = ranges_for_coeff()

    fin = ROOT.TFile.Open(fname,'READ')

    histos = {}
    nbins_y = 0
    nbins_qt = 0
    for q in ['Wplus', 'Wminus']:
        histos[q] = fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name))
        histos[q+'_norm'] = fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name)+'_norm')
        nbins_y = histos[q].GetNbinsX()
        nbins_qt = histos[q].GetNbinsY()

    for qt in range(1, nbins_qt+1):
        c = ROOT.TCanvas('canvas_'+'qt_'+str(qt), "canvas", 800, 600) 
        leg = ROOT.TLegend(0.10,0.70,0.55,0.88, "","brNDC");
        leg.SetFillStyle(0);
        leg.SetBorderSize(0);
        leg.SetTextSize(0.02);
        leg.SetFillColor(10);
        hslices = {} 
        for q in ['Wplus', 'Wminus']:
            h = histos[q] 
            h_norm = histos[q+'_norm']
            qt_bin = 'qt_{:03.0f}'.format(h.GetYaxis().GetBinLowEdge(qt))+'-'+'{:03.0f}'.format(h.GetYaxis().GetBinLowEdge(qt)+h.GetYaxis().GetBinWidth(qt))
            print 'Bin qt: ', qt_bin
            hslice  = h.ProjectionX(str(qt)+'_'+q+'_plus_px', qt, qt)
            hnorm   = h_norm.ProjectionX(str(qt)+'_'+q+'_plus_norm_px', qt, qt)

            # merge bins corresponding to y and -y
            for y in range(nbins_y/2+1, nbins_y+1):
                hslice[y] += hslice[nbins_y+1-y]
                hnorm[y] += hnorm[nbins_y+1-y]

            hslice.Divide(hnorm)
            hslice.SetTitle(coeff[0]+'_{'+coeff[1]+'} for q_{T} #in ['+qt_bin[3:]+'] GeV')
            hslice.SetMinimum(ranges[coeff][0])
            hslice.SetMaximum(ranges[coeff][1])
            hslice.SetStats(0)
            hslice.SetXTitle('|y|')

            if do_fit:
                (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero={'A0': [], 'A1': [0], 'A2': [], 'A3': [0], 'A4': [], 'A5':[], 'A6': [0], 'A7': [0]}, fit_range=[0.0, 3.0])
                if chi2_prob<threshold_chi2:
                    print 'Iterative fit did not converge: try with looser threshold on chi2_prob'
                    (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice_minus, coeff=coeff, 
                                                              fix_to_zero={'A0': [], 'A1': [0], 'A2': [], 'A3': [0], 'A4': [], 'A5':[], 'A6': [0], 'A7': [0]}, fit_range=[0.0, 3.0], threshold_chi2=1e-04)
                hslice.GetFunction('fit_order'+str(order+1)).SetBit(ROOT.TF1.kNotDraw)

            hslice.SetLineWidth(2)
            hslice.SetLineColor(ROOT.kRed if q=='Wplus' else ROOT.kBlue)
            hslice.SetFillColor(ROOT.kRed if q=='Wplus' else ROOT.kBlue)
            hslice.SetFillStyle(3002 if q=='Wplus' else 3003)
            leg_text = 'W^{'+('+' if q=='Wplus' else '-')+'}'

            if do_fit:
                fit.SetLineWidth(2)
                fit.SetLineStyle(ROOT.kDashed)
                fit.SetLineColor(ROOT.kRed if q=='Wplus' else ROOT.kBlue)
                hslices[q+'_fit'] = fit
                leg_text = '#splitline{W^{'+('+' if q=='Wplus' else '-')+'}, #chi^{2}/ndof = '+'{:02.1f}'.format(r.Chi2()/r.Ndf())+'}{'+res+'}'
                leg.AddEntry(fit, 'Fit', "L")

            leg.AddEntry(hslice, leg_text, "f")
            hslices[q] = hslice

        hslices['Wplus'].GetXaxis().SetRangeUser(0.0, 5.0)
        hslices['Wplus'].Draw('E3')
        hslices['Wminus'].Draw('E3SAME')
        if do_fit:
            hslices['Wplus_fit'].Draw('SAME')
            hslices['Wminus_fit'].Draw('SAME')

        leg.Draw()
        c.SaveAs('plots/coefficient_'+var+'_'+coeff+'_'+str(weight_name)+'_'+qt_bin+'.png')
        c.IsA().Destructor( c )        
        leg.IsA().Destructor( leg )        

    fin.Close()


def rebin(h=None, np_bins_template_qt=np.array([]), np_bins_template_y=np.array([]), verbose=False):
    
    np_bins_template_y_extL = np.insert(np_bins_template_y,   0, [h.GetXaxis().GetXmin()])
    np_bins_template_y_ext  = np.append(np_bins_template_y_extL, [h.GetXaxis().GetXmax()])
    np_bins_template_qt_ext = np.append(np_bins_template_qt,     h.GetYaxis().GetXmax())
    nbins_template_y  = np_bins_template_y_ext.size - 1 
    nbins_template_qt = np_bins_template_qt_ext.size - 1 

    h_rebin = ROOT.TH2D(h.GetName()+'_rebin', h.GetTitle(), nbins_template_y, array( 'f',  np_bins_template_y_ext ), nbins_template_qt, array( 'f',  np_bins_template_qt_ext ))

    for iy in range(1, h_rebin.GetNbinsX()+1):
        (y_low, y_high) = (h_rebin.GetXaxis().GetBinLowEdge(iy), h_rebin.GetXaxis().GetBinLowEdge(iy)+h_rebin.GetXaxis().GetBinWidth(iy))

        iys = []
        for iyy in range(1, h.GetNbinsX()+1):
            (yy_low, yy_high) = (h.GetXaxis().GetBinLowEdge(iyy), h.GetXaxis().GetBinLowEdge(iyy)+h.GetXaxis().GetBinWidth(iyy))
            if (yy_low>y_low or np.isclose(yy_low,y_low)) and (yy_high<y_high or np.isclose(yy_high,y_high)):
                if verbose:
                    print '\t\t[',yy_low,',',yy_high, '] < [', y_low,',',y_high, ']'
                iys.append(iyy)

        for iqt in range(1, h_rebin.GetNbinsY()+1):
            (qt_low, qt_high) = (h_rebin.GetYaxis().GetBinLowEdge(iqt), h_rebin.GetYaxis().GetBinLowEdge(iqt)+h_rebin.GetYaxis().GetBinWidth(iqt))

            iqts = []
            for iqqt in range(1, h.GetNbinsY()+1):
                (qqt_low, qqt_high) = (h.GetYaxis().GetBinLowEdge(iqqt), h.GetYaxis().GetBinLowEdge(iqqt)+h.GetYaxis().GetBinWidth(iqqt))
                if (qqt_low>qt_low or np.isclose(qqt_low,qt_low)) and (qqt_high<qt_high or np.isclose(qqt_high,qt_high)):
                    if verbose:
                        print '\t\t[',qqt_low,',',qqt_high, '] < [', qt_low,',',qt_high, ']'
                    iqts.append(iqqt)

            val = 0.
            err2 = 0.
            for iyy in iys:
                for iqqt in iqts:
                    val += h.GetBinContent(iyy,iqqt)
                    err2 += math.pow(h.GetBinError(iyy,iqqt), 2.0)
            h_rebin.SetBinContent(iy,iqt, val)
            h_rebin.SetBinError(iy,iqt, math.sqrt(err2))
            
    return (h_rebin, np_bins_template_qt_ext, np_bins_template_y_ext)
    


def get_covariance(fname='./tree.root', DY='CC', q='Wplus', var='Wdress', 
                   coefficients=['A0'], weights={}, add_stat_uncert=False, postfix='',
                   fix_to_zero={}, fit_range=[0.0, 50.0], threshold_chi2=0.02, verbose=False, 
                   save_corr=True, save_coeff=True, save_tree=True, save_pkl=True,
                   forced_orders={},
                   np_bins_template_qt=np.array([]), np_bins_template_y=np.array([]),
                   plot_updown=False):

    # bins used for the (y,qt) plots
    if np_bins_template_qt.size==0 or np_bins_template_y.size==0:
        print 'Use default binning'
        from tree_utils import np_bins_qt, np_bins_y

    # tree saving the results of the fit for the PDF replicas and scales
    fout = ROOT.TFile.Open('plots/'+'covariance_'+DY+'_'+q+'_'+var+'_'+postfix+'.root', 'RECREATE')
    tree = ROOT.TTree('cov','cov')    
    variables = {}

    # fit & bin-by-bin results (to be saved in a pkl file)
    results = {}
    results['fit_range'] = fit_range

    # covariance matrix from the fit (stat. only)
    fit_covariances = {}

    # map between coeff_bin_y and order of polynomial
    orders = {}
    pvalues = {}

    # input file
    fin = ROOT.TFile.Open(fname,'READ')

    # first loop: determine order of polynomials and book tree branches
    # n_vars = total number of variables = [order+1]*len(coefficients)*n_bins_y
    n_vars = 0
    for coeff in coefficients:

        (h, h_norm) = (None, None)
        if  (np_bins_template_qt.size==0 or np_bins_template_y.size==0):
            (h, h_norm) = (fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(0)), 
                           fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(0)+'_norm')) 
        else:
            (rebinned, rebinned_norm) = (rebin(h=fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(0)), 
                                               np_bins_template_qt=np_bins_template_qt, np_bins_template_y=np_bins_template_y ),
                                         rebin(h=fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(0)+'_norm'), 
                                               np_bins_template_qt=np_bins_template_qt, np_bins_template_y=np_bins_template_y ) )
            (h, h_norm) = (rebinned[0], rebinned_norm[0])
            (np_bins_qt, np_bins_y) = (rebinned[1], rebinned[2])

        np_bins_qt_width = np.array( [np_bins_qt[i+1]-np_bins_qt[i] for i in range(np_bins_qt.size-1)] )
        np_bins_qt_mid = np.array( [(np_bins_qt[i+1]+np_bins_qt[i])*0.5 for i in range(np_bins_qt.size-1)] )
        np_bins_qt_mid_from_zero = np.append(np.array([0.0]), np_bins_qt_mid)
        np_bins_y_width  = np.array( [np_bins_y[i+1]-np_bins_y[i] for i in range(np_bins_y.size-1)] )
        np_bins_y_mid    = np.array( [(np_bins_y[i+1]+np_bins_y[i])*0.5 for i in range(np_bins_y.size-1)] )        
        nbins_y  = np_bins_y.size - 1 
        nbins_qt = np_bins_qt.size - 1                        

        for y in range(nbins_y/2+1, nbins_y+1):
            y_bin = 'y{:03.2f}'.format(np_bins_y[y-1])+'_'+'y{:03.2f}'.format(np_bins_y[y])
            print 'Bin Y: ', y_bin
            name = q+'_'+str(0)
            (hslice, hslice_plus, hnorm, hnorm_plus) = (h.ProjectionY(str(y)+'_'+name+'_py', nbins_y+1-y, nbins_y+1-y),
                                                        h.ProjectionY(str(y)+'_'+name+'_plus_py', y, y),
                                                        h_norm.ProjectionY(str(y)+'_'+name+'_norm_py', nbins_y+1-y, nbins_y+1-y),
                                                        h_norm.ProjectionY(str(y)+'_'+name+'_plus_norm_py', y, y))
            print  '(', nbins_y+1-y, ' + ', y, ') /', nbins_y
            hslice.Add(hslice_plus)
            hnorm.Add(hnorm_plus)
            hslice.Divide(hnorm)

             # make the fit to decide the pol order
            if len(forced_orders)>0:
                (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=fix_to_zero, fit_range=fit_range, force_order=forced_orders[coeff])
            else:
                (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=fix_to_zero, fit_range=fit_range)
                if chi2_prob<threshold_chi2:
                    print 'Iterative fit did not converge: try with looser threshold on chi2_prob'
                    (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=fix_to_zero, fit_range=fit_range, threshold_chi2=1e-04)

            orders[coeff+'_'+y_bin] = order
            pvalues[coeff+'_'+y_bin] = chi2_prob
            # protection against fit with zero parameters (cov = ())
            fit_covariances[coeff+'_'+y_bin] = r.GetCovarianceMatrix() if r.NFreeParameters()>0 else ROOT.TMatrixDSym(1)

            # save the point values with errors
            results[coeff+'_'+y_bin+'_val'] = [hslice.GetBinContent(i+1) for i in range(hslice.GetNbinsX())] 
            results[coeff+'_'+y_bin+'_val_err'] = [hslice.GetBinError(i+1) for i in range(hslice.GetNbinsX())] 
            results[coeff+'_'+y_bin+'_fit'] = [] 

            # create tree branches
            for o in range(order+1):
                nuis_name = coeff+'_'+y_bin+'_p'+str(o) 
                results[coeff+'_'+y_bin+'_fit'].append( r.Parameter(o) ) 
                variables[nuis_name] = array( 'f', [ 0.0 ] )
                variables[nuis_name+'_id'] = array( 'i', [ 0 ] )
                tree.Branch(nuis_name, variables[nuis_name], nuis_name+'/F')
                tree.Branch(nuis_name+'_id', variables[nuis_name+'_id'], nuis_name+'_id'+'/I')
                n_vars += 1

        if  (np_bins_template_qt.size!=0 or np_bins_template_y.size!=0):
            h.IsA().Destructor( h )
            h_norm.IsA().Destructor( h_norm )


    # map of np matrix of all fit-coefficients for all weights
    data = {}
    for syst in ['scale', 'pdf']:
        ws =  weights[syst]
        data[syst] = np.zeros( (n_vars, len(ws)) )
        for iw,w in enumerate(ws):
            vars_count = 0         
            for coeff in coefficients:

                (h, h_norm) = (None, None)
                if  (np_bins_template_qt.size==0 or np_bins_template_y.size==0):
                    (h, h_norm) = (fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)), 
                                   fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)+'_norm')) 
                else:
                    (rebinned, rebinned_norm) = (rebin(h=fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)), 
                                                       np_bins_template_qt=np_bins_template_qt, np_bins_template_y=np_bins_template_y ),
                                                 rebin(h=fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)+'_norm'), 
                                                       np_bins_template_qt=np_bins_template_qt, np_bins_template_y=np_bins_template_y ) )
                    (h, h_norm) = (rebinned[0], rebinned_norm[0])
                    (np_bins_qt, np_bins_y) = (rebinned[1], rebinned[2])
                    
                np_bins_qt_width = np.array( [np_bins_qt[i+1]-np_bins_qt[i] for i in range(np_bins_qt.size-1)] )
                np_bins_qt_mid = np.array( [(np_bins_qt[i+1]+np_bins_qt[i])*0.5 for i in range(np_bins_qt.size-1)] )
                np_bins_qt_mid_from_zero = np.append(np.array([0.0]), np_bins_qt_mid)
                np_bins_y_width  = np.array( [np_bins_y[i+1]-np_bins_y[i] for i in range(np_bins_y.size-1)] )
                np_bins_y_mid    = np.array( [(np_bins_y[i+1]+np_bins_y[i])*0.5 for i in range(np_bins_y.size-1)] )        
                nbins_y  = np_bins_y.size - 1 
                nbins_qt = np_bins_qt.size - 1                        

                for y in range(nbins_y/2+1, nbins_y+1):
                    y_bin = 'y{:03.2f}'.format(np_bins_y[y-1])+'_'+'y{:03.2f}'.format(np_bins_y[y])
                    name = q+'_'+str(w)
                    #(h, h_norm) = (fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)),
                    #               fin.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)+'_norm'))
                    (hslice, hslice_plus, hnorm, hnorm_plus) = (h.ProjectionY(str(y)+'_'+name+'_py', nbins_y+1-y, nbins_y+1-y),
                                                                h.ProjectionY(str(y)+'_'+name+'_plus_py', y, y),
                                                                h_norm.ProjectionY(str(y)+'_'+name+'_norm_py', nbins_y+1-y, nbins_y+1-y),
                                                                h_norm.ProjectionY(str(y)+'_'+name+'_plus_norm_py', y, y))
                    hslice.Add(hslice_plus)
                    hnorm.Add(hnorm_plus)
                    hslice.Divide(hnorm)
                    print 'Syst:', syst, 'weight:', str(w), 'coeff:', coeff, 'y_bin:', y_bin
                    (r,order,fit,res,chi2_prob) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=fix_to_zero, fit_range=fit_range, 
                                                              force_order=orders[coeff+'_'+y_bin])                        
                    for o in range(order+1):
                        nuis_name = coeff+'_'+y_bin+'_p'+str(o) 
                        variables[nuis_name][0] = r.Parameter(o)
                        variables[nuis_name+'_id'][0] = int(w)
                        data[syst][vars_count][iw] = r.Parameter(o)
                        vars_count += 1

                if  (np_bins_template_qt.size!=0 or np_bins_template_y.size!=0):
                    h.IsA().Destructor( h )
                    h_norm.IsA().Destructor( h_norm )                    

            # fill the tree
            tree.Fill()    

    # covariance matrix for n_vars patrameters
    print 'Making the covariance matrix...'
    cov_map = {}
    dict_cov_map = {}
    for syst in ['pdf', 'scale', 'stat']:
        # statistical one is formed from 
        # the cov matrix of the fit
        print '\t> '+syst
        if syst=='stat':
            cov_map['stat']  = np.zeros((n_vars,n_vars))
            vars_count1 = 0 
            for coeff1 in coefficients:
                for y1 in range(nbins_y/2+1, nbins_y+1):
                    y_bin1 = 'y{:03.2f}'.format(np_bins_y[y1-1])+'_'+'y{:03.2f}'.format(np_bins_y[y1])
                    for o1 in range(orders[coeff1+'_'+y_bin1]+1):
                        nuis_name1 = coeff1+'_'+y_bin1
                        if not dict_cov_map.has_key(y_bin1+'_'+coeff1+'_'+'pol'+str(orders[coeff1+'_'+y_bin1])+'_p'+str(o1)):
                            dict_cov_map[y_bin1+'_'+coeff1+'_'+'pol'+str(orders[coeff1+'_'+y_bin1])+'_p'+str(o1)] = vars_count1
                        vars_count2 = 0 
                        for coeff2 in coefficients:
                            for y2 in range(nbins_y/2+1, nbins_y+1):
                                y_bin2 = 'y{:03.2f}'.format(np_bins_y[y2-1])+'_'+'y{:03.2f}'.format(np_bins_y[y2])
                                for o2 in range(orders[coeff2+'_'+y_bin2]+1):
                                    nuis_name2 = coeff2+'_'+y_bin2
                                    if  nuis_name1==nuis_name2:
                                        cov_map['stat'][vars_count1,vars_count2] += fit_covariances[coeff1+'_'+y_bin1][o1][o2]
                                    vars_count2 += 1
                        vars_count1 += 1

        # systematics formed from 
        # the pdf and scale variations
        else:                
            cov_map[syst] = np.cov(data[syst])

        # make a snapshot of the matrix
        if save_corr:
            syst_label = (DY+'_'+q+'_'+var+'_'+syst)+'_'+postfix
            plot_cov_matrix(n_vars=n_vars, cov=cov_map[syst], label='correlation_'+syst_label, plot='corr')

    # total covariance matrix
    cov_map['sum']  = np.zeros((n_vars,n_vars))
    cov_map['syst'] = np.zeros((n_vars,n_vars))
    for syst in ['pdf', 'scale', 'stat']:
        cov_map['sum'] += cov_map[syst]
        if syst in ['pdf', 'scale']:
            cov_map['syst'] += cov_map[syst]
    # make a snapshot of the matrix
    if save_corr:
        cov_label = (DY+'_'+q+'_'+var+'_stat_plus_syst')+'_'+postfix
        plot_cov_matrix(n_vars=n_vars, cov=cov_map['sum'], label='correlation_'+cov_label, plot='corr')

    # save firt results as a pkl file
    if save_pkl:
        pkl_name = 'fit_results_'+DY+'_'+q+'_'+var+'_'+postfix+'.pkl'
        print 'Save results to pickle file....'+'plots/'+pkl_name
        import pickle
        pickle.dump(results, open('plots/'+pkl_name,'wb') )

    bin_count = 0 
    last_bin = np.where(np_bins_qt_mid_from_zero<=fit_range[1])[0][-1] 

    print 'Making plots per coefficient/bin_y...'
    for coeff in coefficients:
        if not save_coeff:
            continue
        print '\t> '+coeff
        for y in range(nbins_y/2+1, nbins_y+1):
            y_bin = 'y{:03.2f}'.format(np_bins_y[y-1])+'_'+'y{:03.2f}'.format(np_bins_y[y])
            bin_name = coeff+'_'+y_bin
            order = orders[bin_name]
            pvalue = pvalues[bin_name]
            print 'Taking sub-matrix: [', bin_count, ',' , bin_count+order , ']' 
            y     = results[bin_name+'_val']
            y_err = results[bin_name+'_val_err']
            p     = results[bin_name+'_fit']
            print 'Save fit:', bin_name, 'at order :', order

            if verbose:
                print '\t Val         :', y
                print '\t Err         :', y_err
                print '\t Coefficients:', p            
                print cov_map['scale'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] 
                print cov_map['pdf'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] 
                print cov_map['stat'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] 

            plt.figure()
            fig, ax = plt.subplots()
            ntoys = 100
            x = np_bins_qt_mid_from_zero[0:last_bin+1]

            y_rnd = np.zeros( (x.size, ntoys) )

            # 1sigma CL from STAT+PDF+SCALE
            for itoy in range(ntoys):
                p_rnd_sum = np.random.multivariate_normal(p, cov_map['sum'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] )
                y_rnd[:, itoy] = polynomial(x=x, coeff=p_rnd_sum, order=order)
            ax.fill_between(x,  polynomial(x=x, coeff=p, order=order)-np.std(y_rnd, axis=1), polynomial(x=x, coeff=p, order=order)+np.std(y_rnd, axis=1), 
                            color='y', linestyle='-', label=r'PDF $\oplus$ scale $\oplus$ stats.')

            # 1sigma CL from scale
            if not plot_updown:
                for itoy in range(ntoys):
                    p_rnd_scale = np.random.multivariate_normal(p, cov_map['scale'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] )
                    y_rnd[:, itoy] = polynomial(x=x, coeff=p_rnd_scale, order=order)
                ax.fill_between(x,  polynomial(x=x, coeff=p, order=order)-np.std(y_rnd, axis=1), polynomial(x=x, coeff=p, order=order)+np.std(y_rnd, axis=1), 
                                color='b', linestyle='-', label=r'Scale ($\mu_R$, $\mu_F$)' )

            # 1sigma CL from PDF
            for itoy in range(ntoys):
                p_rnd_pdf = np.random.multivariate_normal(p, cov_map['pdf'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] )
                y_rnd[:, itoy] = polynomial(x=x, coeff=p_rnd_pdf, order=order)
            ax.fill_between(x,  polynomial(x=x, coeff=p, order=order)-np.std(y_rnd, axis=1), polynomial(x=x, coeff=p, order=order)+np.std(y_rnd, axis=1), 
                            color='g', linestyle='-', label=r'PDF (replicas)')

            # central fit
            ax.plot(x, polynomial(x=x, coeff=p, order=order), 'r--', 
                    label=r'Fit ($\mathrm{pol}_{'+str(order)+'}$), $p$-value: '+'{:0.2f}'.format(pvalue), linewidth=3.0)

            # histogram
            ax.errorbar(np_bins_qt_mid[0:last_bin], y[0:last_bin], 
                        xerr=np_bins_qt_width[0:last_bin]/2, 
                        yerr=y_err[0:last_bin], fmt='o', color='black', label='$'+coeff[0]+'_{'+coeff[1]+'}$')            

            if plot_updown:
                scale_up   = data['scale'][bin_count:(bin_count+order+1), 0]
                scale_down = data['scale'][bin_count:(bin_count+order+1), 1]            
                ax.fill_between(x,  polynomial(x=x, coeff=scale_down, order=order), polynomial(x=x, coeff=scale_up, order=order))
                print '********************************************'
                for iscale,scale in enumerate(weights['scale']):  
                    print  scale, '=>', polynomial(x=x[1:6], coeff=data['scale'][bin_count:(bin_count+order+1), iscale] , order=order)/ \
                        polynomial(x=x[1:6], coeff=p, order=order)

            plt.axis( [0.0, np_bins_qt[last_bin]] + ranges_for_coeff_zoom(q=q)[coeff] )
            plt.grid(True)

            legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
            plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
            plt.ylabel('$'+coeff[0]+'_{'+coeff[1]+'}$', fontsize=20)
            plt.title(DY+', charge='+q[1:]+', var='+var[1:]+', $|y| \in ['+y_bin[1:5]+','+y_bin[7:11]+']$', fontsize=20)
            plt.show()
            plt.savefig('plots/coefficient_'+DY+'_'+q+'_'+var+'_'+coeff+'_'+y_bin+'_fit.png')
            plt.close('all')            
            bin_count += (order+1)

    # save cov matric as np array
    for key,p in cov_map.items():
        print 'Save covariance matrix for', key
        cov_label = (DY+'_'+q+'_'+var+'_'+key)+'_'+postfix
        np.save('plots/covariance_'+cov_label, p)
    pickle.dump(dict_cov_map, open('plots/covariance_dict_'+(DY+'_'+q+'_'+var)+'_'+postfix+'.pkl','wb') )

    # save output tree
    fout.cd()
    if save_tree:
        tree.Write("cov", ROOT.TObject.kOverwrite)
    fout.Close()
    fin.Close()
    


# compare fit results
def compare_fit_results(DYs=['CC_MG5', 'CC_FxFx'], charges=['Wplus', 'Wminus'], variables=['Wdress', 'Wbare', 'WpreFSR'], coefficients=['A0'], fit_range=[0.0, 50.0], postfix='', compare=''):

    import pickle
    from tree_utils import np_bins_qt, np_bins_y, np_bins_qt_width, np_bins_y_width, np_bins_qt_mid, np_bins_y_mid
    np_bins_qt_mid_from_zero = np.append(np.array([0.0]), np_bins_qt_mid)
    nbins_y  = np_bins_y.size - 1 
    nbins_qt = np_bins_qt.size - 1 

    last_bin = np.where(np_bins_qt_mid_from_zero<=fit_range[1])[0][-1]
    x = np_bins_qt_mid_from_zero[0:last_bin+1]

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    fmts   = ['o', 'v', '^', '>', '<', '.', '+', '-']

    for coeff in coefficients:
        print '\>Doing coef '+coeff
        for y in range(nbins_y/2+1, nbins_y+1):
            y_bin = 'y{:03.2f}'.format(np_bins_y[y-1])+'_'+'y{:03.2f}'.format(np_bins_y[y])
            bin_name = coeff+'_'+y_bin        
            print '\t>Doing bin '+y_bin
            plt.figure()
            fig, ax = plt.subplots()
            counter = 0
            for iDY,DY in enumerate(DYs):
                for iq,q in enumerate(charges):
                    plt.axis( [0.0, np_bins_qt[last_bin]] + ranges_for_coeff_zoom(q=q)[coeff] )
                    for ivar,var in enumerate(variables):
                        pkl_name = 'fit_results_'+DY+'_'+q+'_'+var+'_'+postfix+'.pkl'
                        fin = open('plots/'+pkl_name,'r')
                        results = pickle.load(fin)
                        y     = results[bin_name+'_val']
                        y_err = results[bin_name+'_val_err']
                        #p     = results[bin_name+'_fit']
                        #order = len(p)-1
                        label = DY[3:]+' '+('$W' if 'CC' in DY else '$Z')+('^{+}$' if q[1:]=='plus' else '^{-}$')+' '+var[1:]
                        ax.errorbar(np_bins_qt_mid[0:last_bin], y[0:last_bin], xerr=np_bins_qt_width[0:last_bin]/2, yerr=y_err[0:last_bin], fmt=fmts[counter], color=colors[counter], label=label)            
                        fin.close()
                        counter += 1

            plt.grid(True)
            legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
            plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
            plt.ylabel('$'+coeff[0]+'_{'+coeff[1]+'}$', fontsize=20)
            plt.title('$|y| \in ['+y_bin[1:5]+','+y_bin[7:11]+']$', fontsize=20)
            plt.show()
            plt.savefig('plots/coefficient_compare_'+coeff+'_'+y_bin+'_'+compare+'.png')
            plt.close('all')






def fit_BreitWigner(tree=ROOT.TTree(), running=0, var='', cut='', tag=''):

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
        c.SaveAs('plots/fit_BW_'+('run-width_' if running else 'nonrun-width_')+var+'_'+tag+'.png')
    c.IsA().Destructor( c )
    leg.IsA().Destructor( leg )
    return (r.Parameter(2), r.ParError(2), r.Parameter(1), r.ParError(1))
