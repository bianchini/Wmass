import copy
import math
from array import array
import numpy as np 
from scipy import stats
import ROOT
from ROOT import TLorentzVector

def isFromW(p):
    mother = p
    while(mother.numberOfMothers()>0):
        if abs(mother.pdgId())==24:
            return True
        mother = mother.mother(0)
    return False

def isFromZ(p):
    mother = p
    while(mother.numberOfMothers()>0):
        if abs(mother.pdgId())==23:
            return True
        mother = mother.mother(0)
    return False

def isMuon(p):
    if not (p.isPromptFinalState() and abs(p.pdgId())==13 and p.pt()>=0.0):
        return False
    return True
    #return isFromW(p)
    
def isNeutrino(p):
    if not (p.isPromptFinalState() and abs(p.pdgId())==14):
        return False
    return True
    #return isFromW(p)

def isPhoton(p):
    return (p.isPromptFinalState() and p.pdgId()==22 and p.pt()>=0.0)

def deltaR(a,b):
    return math.sqrt( math.pow(a.eta()-b.eta(),2) + math.pow( math.acos( math.cos(a.phi()-b.phi())),2) )

def azimuth(phi):
    if phi<0.0:
        phi += 2*math.pi
    return phi

def boost_to_CS_matrix(Lp4=ROOT.TLorentzVector(0,0,0,0), 
                       Wp4=ROOT.TLorentzVector(0,0,0,0) ):

    Lp4_rot = copy.deepcopy(Lp4)

    # rotate so that W is aligned along x
    Lp4_rot.RotateZ( -Wp4.Phi() )     

    x_lab = np.array([Lp4_rot.E(), Lp4_rot.Px(), Lp4_rot.Py(), Lp4_rot.Pz()])
    E  = Wp4.E()
    qt = Wp4.Pt()
    pz = Wp4.Pz()
    M  = Wp4.M()
    Xt = math.sqrt(M*M + qt*qt)
    boost = np.array([[ E/M, -qt/M, 0, -pz/M ],
                      [ -qt*E/M/Xt, Xt/M, 0, qt*pz/M/Xt],
                      [0., 0., 1., 0.],
                      [-pz/Xt, 0, 0, E/Xt]
                      ] )
    xCS = np.linalg.multi_dot([boost, x_lab])
    flip_z = -1 if Wp4.Rapidity()<0.0 else +1
    ps = (xCS[0], xCS[3]/math.sqrt(xCS[1]*xCS[1] + xCS[2]*xCS[2] + xCS[3]*xCS[3])*flip_z, azimuth(math.atan2(xCS[2],xCS[1])*flip_z) )
    return ps

def boost_to_CS_root(Lp4=ROOT.TLorentzVector(0,0,0,0), 
                     Wp4=ROOT.TLorentzVector(0,0,0,0) ):
    
    Wp4_rot = copy.deepcopy(Wp4)
    Lp4_rot = copy.deepcopy(Lp4)

    # align W/L along x axis
    Wp4_rot.RotateZ( -Wp4.Phi() )
    Lp4_rot.RotateZ( -Wp4.Phi() )

    # first boost
    boostL = Wp4_rot.BoostVector()
    boostL.SetX(0.0)
    boostL.SetY(0.0)
    Lp4_rot.Boost( -boostL )
    Wp4_rot.Boost( -boostL )

    # second boost
    boostT = Wp4_rot.BoostVector()
    Lp4_rot.Boost( -boostT )

    # the CS frame defines the z-axis according to the W pz in the lab 
    flip_z = -1 if Wp4.Rapidity()<0.0 else +1

    # compute PS point
    ps = (Lp4_rot.E(), Lp4_rot.CosTheta()*flip_z, azimuth(Lp4_rot.Phi()*flip_z) )
    return ps

def printp(tag, p, other):
    print tag+' ['+str(p.pdgId())+']: ('+'{:04.3f}'.format(p.pt())+',{:04.3f}'.format(p.eta())+',{:04.3f}'.format(p.phi())+') .... '+other

def add_vars(tree=None, debug=False):
    names = ['nuLost', 'muLost', 'mu_charge', 'nu_charge', 'isFromW', 'alphaQED', 'alphaQCD', 'scale', 'id1', 'id2', 'x1', 'x2']
    for v in ['qt', 'y', 'mass', 'phi']:
        names.append('lhe_'+v)
    for t in ['WpreFSR', 'Wdress', 'Wbare']:
        for v in ['qt', 'y', 'mass', 'phi', 'ECS','cosCS', 'phiCS', 'mu_pt', 'mu_eta', 'mu_phi', 'nu_pt', 'nu_eta', 'nu_phi' ]:
            names.append(t+'_'+v)    
    if debug and tree!=None:
        nevents = 3
        scan = 'id1:x1:id2:x2:alphaQCD:alphaQED:scale'
        tree.Scan(scan, "", "", nevents)
        scan = 'nuLost:muLost:isFromW:mu_charge:nu_charge'
        tree.Scan(scan, "", "", nevents)
        scan = 'lhe_mass:lhe_qt:lhe_y:lhe_phi:Wdress_mass:Wdress_qt:Wdress_y:Wdress_phi'
        tree.Scan(scan, "", "", nevents)
        scan = 'Wdress_mu_pt:Wdress_mu_eta:Wdress_nu_pt:Wdress_nu_eta:Wdress_ECS:Wdress_cosCS:Wdress_phiCS'
        tree.Scan(scan, "", "", nevents)
        return

    variables = {}
    for name in names:        
        variables[name] = np.zeros(1, dtype=float)
        if tree!=None:
            tree.Branch(name, variables[name], name+'/D')
    variables['weights'] = np.zeros(109, dtype=float)
    if tree!=None:
        tree.Branch('weights', variables['weights'], 'weights[109]/D')

    return variables
    
def fill_default(variables):
    for key,var in variables.items():
        var[0] = 0.0


def add_histo2D(charges=['Wminus','Wplus'], var=['Wdress'], coeff=['A0','A1','A2','A3','A4','A5', 'A6', 'A7'], weights=[0]):

    # binning
    np_bins_qt = np.append(np.append(np.linspace(0.0, 10.0, 6), np.linspace(12.0, 20.0, 5)), np.append(np.linspace(24.0, 40.0, 5), np.array([60,80,100,150, 200]) ) )
    np_bins_y  = np.append(np.append(np.linspace(-5.0, -2.5, 6), np.linspace(-2.0, +2.0, 21)), np.linspace(2.5, 5.0, 6))
    bins_qt = array( 'f',  np_bins_qt )
    bins_y  = array( 'f',  np_bins_y )
    
    histos = {}

    # create TH2D
    for q in charges:
        histos[q] = {}
        for v in var:
            histos[q][v] = {}
            for c in coeff:
                histos[q][v][c] = {}
                for w in weights:
                    name = q+'_'+v+'_'+c+'_'+str(w)
                    histos[q][v][c][str(w)] = (ROOT.TH2D(name, name, len(bins_y)-1, bins_y, len(bins_qt)-1, bins_qt),
                                               ROOT.TH2D(name+'_norm', name+'_norm', len(bins_y)-1, bins_y, len(bins_qt)-1, bins_qt))
                    # create sum w2
                    histos[q][v][c][str(w)][0].Sumw2()
                    histos[q][v][c][str(w)][1].Sumw2()
        
    return histos

def test_A(coeff='A0', ps=(0.0, 0.0)):

    val = 0.0
    if coeff=='A0':
        val = 20./3.*(0.5*(1-3*ps[0]*ps[0])) + 2./3.
    elif coeff=='A1':
        val = 5.*2.*math.sqrt(1-ps[0]*ps[0])*ps[0]*math.cos(ps[1])
    elif coeff=='A2':
        val = 10.*(1-ps[0]*ps[0])*math.cos(2*ps[1]) 
    elif coeff=='A3':
        val = 4.*math.sqrt(1-ps[0]*ps[0])*math.cos(ps[1])
    elif coeff=='A4':
        val = 4.*ps[0]
    elif coeff=='A5':
        val = 5.*(1-ps[0]*ps[0])*math.sin(2*ps[1])
    elif coeff=='A6':
        val = 5.*2.*math.sqrt(1-ps[0]*ps[0])*ps[0]*math.sin(ps[1])
    elif coeff=='A7':
        val = 4.*math.sqrt(1-ps[0]*ps[0])*math.sin(ps[1])
    return val


def fill_coefficients(histos={}, charge=0, var='Wdress', weight_name=0, ps_W=(), ps_CS=(), weight=1.0):
    q = 'Wplus' if charge==-13 else 'Wminus'
    for coeff in ['A0','A1','A2','A3','A4','A5', 'A6', 'A7']:
        (h,h_norm) = histos[q][var][coeff][str(weight_name)]
        h.Fill(ps_W[0],ps_W[1], weight*test_A(coeff=coeff,ps=ps_CS) )
        h_norm.Fill(ps_W[0],ps_W[1], weight )

def find_polynomial(h=None, coeff='A0'):
    
    chi2 = 99999.
    order = 1
    r = None
    fit = None
    pval = 0.0
    while( pval<0.05 ):
        formula = ''
        for p in range(order+1):
            formula += ('['+str(p)+']*TMath::Power(x,'+str(p)+')'+('+' if p<order else ''))
        print formula
        this_fit = ROOT.TF1('fit_order'+str(order), formula,  0, 50) 
        for p in range(order+1):
            this_fit.SetParameter(p, math.pow(10,-p))
        if coeff!='A4':
            this_fit.FixParameter(0, 0.0 )
        this_fit.SetNpx(10000)
        this_r = h.Fit('fit_order'+str(order), 'SRQ')
        delta_chi2 = chi2 - this_r.Chi2()
        pval = 1-stats.chi2.cdf(delta_chi2, 1) 
        print '\tOrder', order, 'pval=', pval, 'chi2=', this_r.Chi2(), '/', this_r.Ndf()
        if pval<0.05:            
            r = this_r
            fit = this_fit
            chi2 = this_r.Chi2()        
            order += 1

    print 'P-value threshold at order', order-1, 'with chi2=', r.Chi2(), '/', r.Ndf()
    return (r,order-1,fit)
        

def draw_slice(fname='./tree.root', var='Wdress', coeff='A0', weight_name=0):

    ranges = {}
    ranges['A0'] = (-0.2, 1.5)
    ranges['A1'] = (-0.8, 0.8)
    ranges['A2'] = (-0.2, 1.5)
    ranges['A3'] = (-2.0, 2.0)
    ranges['A4'] = (-2.5, 2.5)
    ranges['A5'] = (-0.5, 0.5)
    ranges['A6'] = (-0.5, 0.5)
    ranges['A7'] = (-0.5, 0.5)


    f = ROOT.TFile.Open(fname,'READ')

    histos = {}
    nbins_y = 0
    nbins_qt = 0
    for q in ['Wplus', 'Wminus']:
        histos[q] = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name))
        histos[q+'_norm'] = f.Get(q+'/'+var+'/'+coeff+'/'+q+'_'+var+'_'+coeff+'_'+str(weight_name)+'_norm')
        nbins_y = histos[q].GetNbinsX()
        nbins_qt = histos[q].GetNbinsY()

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
            hslice_minus.SetTitle(coeff[0]+'_{'+coeff[1]+'} for y_{W} #in ['+y_bin[2:]+']')
            hslice_minus.SetMinimum(ranges[coeff][0])
            hslice_minus.SetMaximum(ranges[coeff][1])
            hslice_minus.SetStats(0)
            hslice_minus.SetXTitle('q_{T} (GeV)')

            (r,order,fit) = find_polynomial(h=hslice_minus, coeff=coeff)
            hslice_minus.GetFunction('fit_order'+str(order+1)).SetBit(ROOT.TF1.kNotDraw)

            #fit = ROOT.TF1('fit_'+q, '[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x',  0, 50)
            #fit.SetParameter(0, 0.0)
            #fit.SetParameter(1, 0.1)
            #fit.SetParameter(2, 0.01 )
            #fit.SetParameter(3, 0.001 )
            #fit.SetParameter(3, 0.0001 )
            #if coeff!='A4':
            #    fit.FixParameter(0, 0.0 )
            #fit.SetNpx(10000)
            #r = hslice_minus.Fit('fit_'+q, 'SR')
            
            res = ''
            for p in range(order+1):
                exp = '{:.1E}'.format(r.Parameter(p))[-3]+'{:.1E}'.format(r.Parameter(p))[-1] 
                val = '{:.1E}'.format(r.Parameter(p))[:-4]
                err = '{:.1E}'.format(r.ParError(p))[:-4]
                res += 'c_{'+str(p)+'}=('+val+' #pm '+err+')#times10^{'+exp+'}, '
            if q=='Wplus':                
                hslice_minus.SetLineWidth(2)
                hslice_minus.SetLineColor(ROOT.kRed)
                hslice_minus.SetFillColor(ROOT.kRed)
                hslice_minus.SetFillStyle(3002)
                fit.SetLineWidth(2)
                fit.SetLineStyle(ROOT.kDashed)
                fit.SetLineColor(ROOT.kRed)
                leg.AddEntry(hslice_minus, '#splitline{W^{+}, #chi^{2}/ndof = '+'{:02.1f}'.format(r.Chi2()/r.Ndf())+'}{'+res+'}', "LP")
                hslices[q+'_fit'] = fit
            else:
                hslice_minus.SetLineWidth(2)
                hslice_minus.SetFillColor(ROOT.kBlue)
                hslice_minus.SetFillStyle(3003)
                hslice_minus.SetLineColor(ROOT.kBlue)
                fit.SetLineWidth(2)
                fit.SetLineStyle(ROOT.kDashed)
                fit.SetLineColor(ROOT.kBlue)
                leg.AddEntry(hslice_minus, '#splitline{W^{-}, #chi^{2}/ndof = '+'{:02.1f}'.format(r.Chi2()/r.Ndf())+'}{'+res+'}', "LP")            
                hslices[q+'_fit'] = fit

            hslices[q] = hslice_minus

        hslices['Wplus'].Draw('E3')
        hslices['Wminus'].Draw('E3SAME')
        hslices['Wplus_fit'].Draw('SAME')
        hslices['Wminus_fit'].Draw('SAME')

        leg.Draw()
        c.SaveAs('plots/coefficient_'+var+'_'+coeff+'_'+str(weight_name)+'_'+y_bin+'.png')
        #raw_input()
        c.IsA().Destructor( c )        
        leg.IsA().Destructor( leg )        
    f.Close()
