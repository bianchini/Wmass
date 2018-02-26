import math
import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 8,7

from scipy import stats
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
def Fisher_test(h=None, coeff='A0', fix_to_zero=['A4'], fit_range=[0.0,50.0], force_order=-1):

    # a large initial value
    chi2 = 99999.

    # start from pol1
    order = 1
    if force_order>=0:
        order = force_order

    r = None
    fit = None
    pval = 0.0
    while( pval<0.05 ):

        # standard polynomial
        formula = ''
        for p in range(order+1):
            formula += ('['+str(p)+']*TMath::Power(x,'+str(p)+')'+('+' if p<order else ''))

        #print formula
        this_fit = ROOT.TF1('fit_order'+str(order), formula,  fit_range[0], fit_range[1]) 
        this_fit.SetNpx(10000)
        for p in range(order+1):
            this_fit.SetParameter(p, math.pow(10,-p))
        if coeff in fix_to_zero:
            this_fit.FixParameter(0, 0.0 )

        this_r = h.Fit('fit_order'+str(order), 'SRQ')
        delta_chi2 = chi2 - this_r.Chi2()
        pval = 1-stats.chi2.cdf(delta_chi2, 1) 

        if force_order<0:
            print '\tOrder', order, 'pval=', pval, 'chi2=', this_r.Chi2(), '/', this_r.Ndf()
        else:
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()        
            order += 1 
            break

        if pval<0.05:            
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()        
            order += 1

    if force_order<0:
        print 'P-value threshold at order', order-1, 'with chi2=', r.Chi2(), '/', r.Ndf()
    else:
        print 'Fit performed at order ', order-1, 'with chi2=', r.Chi2(), '/', r.Ndf() 

    # nice formatting of results
    res = ''
    for p in range(order):
        exp = '{:.1E}'.format(r.Parameter(p))[-3]+'{:.1E}'.format(r.Parameter(p))[-1] 
        val = '{:.1E}'.format(r.Parameter(p))[:-4]
        err = '{:.1E}'.format(r.ParError(p))[:-4]
        res += 'c_{'+str(p)+'}=('+val+' #pm '+err+')#times10^{'+exp+'}, '

    # remove 1 from order to go to the previous one
    return (r,order-1,fit,res)

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
def draw_y_slice(fname='./tree.root', var='Wdress', coeff='A0', weight_name=0, do_fit=True):

    f = ROOT.TFile.Open(fname,'READ')

    ranges = ranges_for_coeff()

    histos = {}
    nbins_y = 0
    nbins_qt = 0
    for q in ['Wplus', 'Wminus']:
        histos[q] = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name))
        histos[q+'_norm'] = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name)+'_norm')
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
            y_bin = 'y_{:03.1f}'.format(h.GetXaxis().GetBinLowEdge(y))+'-'+'{:03.1f}'.format(h.GetXaxis().GetBinLowEdge(y)+h.GetXaxis().GetBinWidth(y))
            print 'Bin Y: ', y_bin
            hslice_minus = h.ProjectionY(str(y)+'_'+q+'_minus_py', nbins_y+1-y, nbins_y+1-y)
            hslice_plus  = h.ProjectionY(str(y)+'_'+q+'_plus_py', y, y)
            hnorm_minus  = h_norm.ProjectionY(str(y)+'_'+q+'_minus_norm_py', nbins_y+1-y, nbins_y+1-y)
            hnorm_plus   = h_norm.ProjectionY(str(y)+'_'+q+'_plus_norm_py', y, y)
            print  '(', nbins_y+1-y, ' + ', y, ') /', nbins_y
            hslice_minus.Add(hslice_plus)
            hnorm_minus.Add(hnorm_plus)
            hslice_minus.Divide(hnorm_minus)
            hslice_minus.SetTitle(coeff[0]+'_{'+coeff[1]+'} for |y| #in ['+y_bin[2:]+']')
            hslice_minus.SetMinimum(ranges[coeff][0])
            hslice_minus.SetMaximum(ranges[coeff][1])
            hslice_minus.SetStats(0)
            hslice_minus.SetXTitle('q_{T} (GeV)')

            if do_fit:
                (r,order,fit,res) = Fisher_test(h=hslice_minus, coeff=coeff, fix_to_zero=['A0','A1','A2','A3','A5','A6','A7'], fit_range=[0.0,50.0])
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
                leg_text = '#splitline{W^{'+('+' if q=='Wplus' else '-')+'}, #chi^{2}/ndof = '+'{:02.1f}'.format(r.Chi2()/r.Ndf())+'}{'+res+'}'
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

    f.Close()

# project the (y,qt) plot along the y axis for all bins of qt
def draw_qt_slice(fname='./tree.root', var='Wdress', coeff='A0', weight_name=0, do_fit=True):

    ranges = ranges_for_coeff()

    f = ROOT.TFile.Open(fname,'READ')

    histos = {}
    nbins_y = 0
    nbins_qt = 0
    for q in ['Wplus', 'Wminus']:
        histos[q] = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name))
        histos[q+'_norm'] = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name)+'_norm')
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
                (r,order,fit,res) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=['A1', 'A3', 'A6', 'A7'], fit_range=[0.0, 3.0])
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

    f.Close()


def get_covariance(fname='./tree.root', DY='CC', var='Wdress', q='Wplus', coefficients=['A0'], weights={}, add_stat_uncert=False, postfix='',
                   fix_to_zero=['A0','A1','A2','A3','A5','A6','A7'], fit_range=[0.0, 50.0]):

    # bins used for the (y,qt) plots
    from tree_utils import np_bins_qt, np_bins_y, np_bins_qt_width, np_bins_y_width, np_bins_qt_mid, np_bins_y_mid
    #np_bins_y = np_bins_y[15:18]
    np_bins_qt_mid_from_zero = np.append(np.array([0.0]), np_bins_qt_mid)
    nbins_y = np_bins_y.size - 1 
    nbins_qt = np_bins_qt.size - 1 

    out = ROOT.TFile.Open('plots/'+'covariance_'+q+postfix+'.root', 'RECREATE')
    tree = ROOT.TTree('cov','cov')
    variables = {}

    import pickle
    results = {}
    covariances = {}
    orders = {}

    # input file
    f = ROOT.TFile.Open(fname,'READ')

    # first loop: determine order of polynomials and book tree branches
    n_vars = 0
    for coeff in coefficients:
        for y in range(nbins_y/2+1, nbins_y+1):
            y_bin = 'y_{:03.1f}'.format(np_bins_y[y-1])+'_'+'{:03.1f}'.format(np_bins_y[y])
            print 'Bin Y: ', y_bin
            name = q+'_'+str(0)
            h = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(0))
            h_norm = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(0)+'_norm')
            hslice = h.ProjectionY(str(y)+'_'+name+'_py', nbins_y+1-y, nbins_y+1-y)
            hslice_plus  = h.ProjectionY(str(y)+'_'+name+'_plus_py', y, y)
            hnorm  = h_norm.ProjectionY(str(y)+'_'+name+'_norm_py', nbins_y+1-y, nbins_y+1-y)
            hnorm_plus   = h_norm.ProjectionY(str(y)+'_'+name+'_plus_norm_py', y, y)
            print  '(', nbins_y+1-y, ' + ', y, ') /', nbins_y
            hslice.Add(hslice_plus)
            hnorm.Add(hnorm_plus)
            hslice.Divide(hnorm)
            (r,order,fit,res) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=fix_to_zero, fit_range=fit_range)
            orders[coeff+'_'+y_bin] = order
            covariances[coeff+'_'+y_bin] = r.GetCovarianceMatrix()

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

    # arrays that contain the coefficients for all pdf and scale variations
    data = {}
    for syst in ['pdf', 'scale']:
        ws =  weights[syst]
        data[syst] = np.zeros( (n_vars, len(ws)) )
        for iw,w in enumerate(ws):
            vars_count = 0         
            for coeff in coefficients:
                for y in range(nbins_y/2+1, nbins_y+1):
                    y_bin = 'y_{:03.1f}'.format(np_bins_y[y-1])+'_'+'{:03.1f}'.format(np_bins_y[y])
                    name = q+'_'+str(w)
                    h = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w))
                    h_norm = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(w)+'_norm')
                    hslice = h.ProjectionY(str(y)+'_'+name+'_py', nbins_y+1-y, nbins_y+1-y)
                    hslice_plus  = h.ProjectionY(str(y)+'_'+name+'_plus_py', y, y)
                    hnorm  = h_norm.ProjectionY(str(y)+'_'+name+'_norm_py', nbins_y+1-y, nbins_y+1-y)
                    hnorm_plus   = h_norm.ProjectionY(str(y)+'_'+name+'_plus_norm_py', y, y)
                    hslice.Add(hslice_plus)
                    hnorm.Add(hnorm_plus)
                    hslice.Divide(hnorm)
                    (r,order,fit,res) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=fix_to_zero, fit_range=fit_range, force_order=orders[coeff+'_'+y_bin])
                    for o in range(order+1):
                        nuis_name = coeff+'_'+y_bin+'_p'+str(o) 
                        variables[nuis_name][0] = r.Parameter(o)
                        variables[nuis_name+'_id'][0] = int(w)
                        data[syst][vars_count][iw] = r.Parameter(o)
                        vars_count += 1
            tree.Fill()    

    # covariance matrix
    cov_map = {}
    for syst in ['pdf', 'scale', 'stat']:
        # statistical one is formed from the cov matrix of the fits
        if syst=='stat':
            cov_map['stat']  = np.zeros((n_vars,n_vars))
            vars_count1 = 0 
            for coeff1 in coefficients:
                for y1 in range(nbins_y/2+1, nbins_y+1):
                    y_bin1 = 'y_{:03.1f}'.format(np_bins_y[y1-1])+'_'+'{:03.1f}'.format(np_bins_y[y1])
                    for o1 in range(orders[coeff1+'_'+y_bin1]+1):
                        nuis_name1 = coeff1+'_'+y_bin1
                        vars_count2 = 0 
                        for coeff2 in coefficients:
                            for y2 in range(nbins_y/2+1, nbins_y+1):
                                y_bin2 = 'y_{:03.1f}'.format(np_bins_y[y2-1])+'_'+'{:03.1f}'.format(np_bins_y[y2])
                                for o2 in range(orders[coeff2+'_'+y_bin2]+1):
                                    nuis_name2 = coeff2+'_'+y_bin2
                                    if  nuis_name1==nuis_name2:
                                        cov_map['stat'][vars_count1,vars_count2] += covariances[coeff1+'_'+y_bin1][o1][o2]
                                    vars_count2 += 1
                        vars_count1 += 1
        # systematics formed from the pdf and scal variations
        else:                
            cov_map[syst] = np.cov(data[syst])
        # plot the matrix
        syst_label = (q+'_'+syst)+postfix
        plot_cov_matrix(n_vars=n_vars, cov=cov_map[syst], label=syst_label, plot='corr')

    # total covariance matrix
    cov_map['sum'] = np.zeros((n_vars,n_vars))
    cov_map['syst'] = np.zeros((n_vars,n_vars))
    for syst in ['pdf', 'scale', 'stat']:
        cov_map['sum'] += cov_map[syst]
        if syst in ['pdf', 'scale']:
            cov_map['syst'] += cov_map[syst]
    cov_label = (q+'_stat_plus_syst')+postfix
    plot_cov_matrix(n_vars=n_vars, cov=cov_map['sum'], label=cov_label, plot='corr')

    # save firt results as a pkl file
    pickle.dump(results, open('plots/'+q+'_results'+postfix+'.pkl','wb') )

    bin_count = 0 
    last_bin = 16
    for coeff in coefficients:
        for y in range(nbins_y/2+1, nbins_y+1):
            y_bin = 'y_{:03.1f}'.format(np_bins_y[y-1])+'_'+'{:03.1f}'.format(np_bins_y[y])
            bin_name = coeff+'_'+y_bin
            order = orders[bin_name]
            #print 'Taking sub-matrix: [', bin_count, ',' , bin_count+(order) , ']' 

            y = results[bin_name+'_val']
            y_err = results[bin_name+'_val_err']
            p = results[bin_name+'_fit']
            print 'Save fit:', bin_name, ':', order
            #print '\t Val         :', y
            #print '\t Err         :', y_err
            #print '\t Coefficients:', p

            plt.figure()
            fig, ax = plt.subplots()
            ntoys = 100
            x = np_bins_qt_mid_from_zero[0:last_bin+1]
            for itoy in range(ntoys):
                p_rnd_sum = np.random.multivariate_normal(p, cov_map['sum'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] )
                ax.plot(x, polynomial(x=x, coeff=p_rnd_sum, order=order), 'y-', label=('PDF $\otimes$ scale $\otimes$ stats.'  if itoy==0 else None) )
            for itoy in range(ntoys):
                p_rnd_scale = np.random.multivariate_normal(p, cov_map['scale'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] )
                ax.plot(x, polynomial(x=x, coeff=p_rnd_scale, order=order), 'b-', label=('Scale ($\mu_R$, $\mu_F$)' if itoy==0 else None) )
            for itoy in range(ntoys):
                p_rnd_pdf = np.random.multivariate_normal(p, cov_map['pdf'][bin_count:(bin_count+order+1), bin_count:(bin_count+order+1)] )
                ax.plot(x, polynomial(x=x, coeff=p_rnd_pdf, order=order), 'g-', label=('PDF (replicas)' if itoy==0 else None) )

            ax.plot(x, polynomial(x=x, coeff=p, order=order), 'r--', label='Fit', linewidth=3.0)
            ax.errorbar(np_bins_qt_mid[0:last_bin], y[0:last_bin], xerr=np_bins_qt_width[0:last_bin]/2, yerr=y_err[0:last_bin], fmt='o', color='black', label='$'+coeff[0]+'_{'+coeff[1]+'}$')
            
            plt.axis( [0.0, np_bins_qt[last_bin]] + ranges_for_coeff_zoom(q=q)[coeff] )
            plt.grid(True)

            legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
            plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
            plt.ylabel('$'+coeff[0]+'_{'+coeff[1]+'}$', fontsize=20)
            plt.title(DY+', charge='+q[1:]+', $|y| \in ['+y_bin[2:6]+','+y_bin[6:11]+']$', fontsize=20)
            plt.show()
            plt.savefig('plots/coefficient_'+q+'_'+var+'_'+coeff+'_'+y_bin+'_fit.png')
            plt.close()            
            bin_count += (order+1)

    # save output tree
    out.cd()
    tree.Write("cov", ROOT.TObject.kOverwrite)
    out.Close()
    f.Close()
    


