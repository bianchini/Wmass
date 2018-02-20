import math
import numpy as np 
from scipy import stats
from array import array
import ROOT

# perform Fisher-test based on chi2 of fit. Stop when p-value of chi2>0.05
def Fisher_test(h=None, coeff='A0', fix_to_zero=['A4'], fit_range=[0.0,50.0], force_order=-1):
    
    chi2 = 99999.
    order = 1
    if force_order>=0:
        order = force_order
    r = None
    fit = None
    pval = 0.0
    while( pval<0.05 ):
        formula = ''
        for p in range(order+1):
            formula += ('['+str(p)+']*TMath::Power(x,'+str(p)+')'+('+' if p<order else ''))
        #print formula
        this_fit = ROOT.TF1('fit_order'+str(order), formula,  fit_range[0], fit_range[1]) 
        for p in range(order+1):
            this_fit.SetParameter(p, math.pow(10,-p))
        if coeff in fix_to_zero:
            this_fit.FixParameter(0, 0.0 )

        this_fit.SetNpx(10000)
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

    res = ''
    for p in range(order):
        exp = '{:.1E}'.format(r.Parameter(p))[-3]+'{:.1E}'.format(r.Parameter(p))[-1] 
        val = '{:.1E}'.format(r.Parameter(p))[:-4]
        err = '{:.1E}'.format(r.ParError(p))[:-4]
        res += 'c_{'+str(p)+'}=('+val+' #pm '+err+')#times10^{'+exp+'}, '

    return (r,order-1,fit,res)


def plot_cov_matrix(n_vars=0, cov=None, label='', plot='cov'):
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
    c.SaveAs(label+'_'+plot+'.png')
    c.SaveAs(label+'_'+plot+'.C')
    np.save(label+'_'+plot, cov)
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


def get_covariance(fname='./tree.root', var='Wdress', q='Wplus', coefficients=['A0'], weights=[0], add_stat_uncert=False, postfix=''):

    from tree_utils import np_bins_qt, np_bins_y
    nbins_y = np_bins_y.size - 1 
    nbins_qt = np_bins_qt.size - 1 

    out = ROOT.TFile.Open('covariance_'+q+postfix+'.root', 'RECREATE')
    tree = ROOT.TTree('cov','cov')
    variables = {}
    covariances = {}
    orders = {}

    f = ROOT.TFile.Open(fname,'READ')

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
            (r,order,fit,res) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=['A0','A1','A2','A3','A5','A6','A7'], fit_range=[0.0, 50.0])
            orders[coeff+'_'+y_bin] = order
            covariances[coeff+'_'+y_bin] = r.GetCovarianceMatrix()
            for o in range(order+1):
                nuis_name = coeff+'_'+y_bin+'_p'+str(o) 
                variables[nuis_name] = array( 'f', [ 0.0 ] )
                tree.Branch(nuis_name, variables[nuis_name], nuis_name+'/F')
                n_vars += 1

    data = np.zeros((n_vars, len(weights)))
    for iw,w in enumerate(weights):
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
                (r,order,fit,res) = Fisher_test(h=hslice, coeff=coeff, fix_to_zero=['A0','A1','A2','A3','A5','A6','A7'], fit_range=[0.0, 50.0], force_order=orders[coeff+'_'+y_bin])
                for o in range(order+1):
                    nuis_name = coeff+'_'+y_bin+'_p'+str(o) 
                    variables[nuis_name][0] = r.Parameter(o)
                    data[vars_count][iw] = r.Parameter(o)
                    vars_count += 1
        tree.Fill()

    #print data
    cov_syst = np.zeros((n_vars,n_vars))
    cov_stat = np.zeros((n_vars,n_vars))

    cov_syst += np.cov(data)
    cov_label = (q+'_syst_only')+postfix
    plot_cov_matrix(n_vars=n_vars, cov=cov_syst, label=cov_label, plot='cov')
    plot_cov_matrix(n_vars=n_vars, cov=cov_syst, label=cov_label, plot='corr')

    if add_stat_uncert:
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
                                    stat = covariances[coeff1+'_'+y_bin1][o1][o2]
                                    print nuis_name1, o1, o2, ' ==> adding stat. uncert. ', stat, 'to syst. uncert. ', cov_syst[vars_count1,vars_count2]  
                                    cov_stat[vars_count1,vars_count2] += stat
                                vars_count2 += 1
                    vars_count1 += 1

        cov_label = (q+'_stat_only')+postfix
        plot_cov_matrix(n_vars=n_vars, cov=cov_stat, label=cov_label, plot='cov')
        plot_cov_matrix(n_vars=n_vars, cov=cov_stat, label=cov_label, plot='corr')

        cov = np.zeros((n_vars,n_vars))
        cov += cov_stat
        cov += cov_syst
        cov_label = (q+'_stat_plus_syst')+postfix
        plot_cov_matrix(n_vars=n_vars, cov=cov, label=cov_label, plot='cov')
        plot_cov_matrix(n_vars=n_vars, cov=cov, label=cov_label, plot='corr')

    out.cd()
    tree.Write("cov", ROOT.TObject.kOverwrite)
    out.Close()
    f.Close()
    


