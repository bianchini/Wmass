import copy
import math
from array import array
import numpy as np 
import ROOT
from ROOT import TLorentzVector

# binning for templates
np_bins_qt_p0 = np.linspace( 0.0, 10.0, 6)
np_bins_qt_p1 = np.linspace(12.0, 20.0, 5)
np_bins_qt_p2 = np.linspace(24.0, 40.0, 5)
np_bins_qt_p3 = np.array([60, 80, 100, 150, 200]) 
#np_bins_qt    = np.append( np.append(np_bins_qt_p0, np_bins_qt_p1), np.append( np_bins_qt_p2, np_bins_qt_p3))
np_bins_qt    = np.array([0.0, 32.0])

np_bins_qt_width = np.array( [np_bins_qt[i+1]-np_bins_qt[i] for i in range(np_bins_qt.size-1)] )
np_bins_qt_mid   = np.array( [(np_bins_qt[i+1]+np_bins_qt[i])*0.5 for i in range(np_bins_qt.size-1)] )

np_bins_y_p0 = np.linspace(-5.0, -2.5,  6)
np_bins_y_p1 = np.linspace(-2.0, +2.0, 21)
np_bins_y_p2 = np.linspace(+2.5, +5.0,  6)
#np_bins_y    = np.append( np.append(np_bins_y_p0, np_bins_y_p1), np_bins_y_p2)
np_bins_y = np.array([-3.5, -2.5, 0.0, 2.5, 3.5])

np_bins_y_width = np.array( [np_bins_y[i+1]-np_bins_y[i] for i in range(np_bins_y.size-1)] )
np_bins_y_mid   = np.array( [(np_bins_y[i+1]+np_bins_y[i])*0.5 for i in range(np_bins_y.size-1)] )

np_bins_pt  = np.linspace( 25.0, 65.0, 81  )
np_bins_ptinv = np.linspace( 1./65.0, 1./25.0, 81  )
np_bins_eta = np.linspace( -2.5, 2.5,  101 )


# find bin corresponding to a given ps=(|y|,qt) point
# if an over/underflow is found, return 'OF' 
def find_y_qt_bin( ps=(), verbose=False ):

    y = abs(ps[0])
    iy_low  = np.where(np_bins_y<=y)[0][-1] if np.where(np_bins_y<=y)[0].size>0 else -1
    iy_high = np.where(np_bins_y>y)[0][0]   if np.where(np_bins_y>y)[0].size>0  else -1
    bin_y='OF'
    if iy_low==-1 or iy_high==-1:
        if verbose:
            print 'y=', y, 'yields (', iy_low, iy_high, ') => return'
    else:
        bin_y = 'y{:03.2f}'.format(np_bins_y[iy_low])+'_'+'y{:03.2f}'.format(np_bins_y[iy_high])

    qt = ps[1]
    iqt_low = np.where(np_bins_qt<=qt)[0][-1] if np.where(np_bins_qt<=qt)[0].size>0 else -1
    iqt_high = np.where(np_bins_qt>qt)[0][0] if np.where(np_bins_qt>qt)[0].size>0 else -1
    bin_qt='OF'
    if iqt_low==-1 or iqt_high==-1:
        if verbose:
            print 'qt=', qt, 'yields (', iqt_low, iqt_high, ') => return'
    else:
        bin_qt = 'qt{:03.1f}'.format(np_bins_qt[iqt_low])+'_'+'qt{:03.1f}'.format(np_bins_qt[iqt_high])

    if verbose:
        print 'y=', y, 'yields bins (', iy_low, iy_high, ') => Ok'
        print 'qt=', qt, 'yields bins (', iqt_low, iqt_high, ') => Ok'
    
    return (bin_y, iy_low, bin_qt, iqt_low)

def get_weight_meaning(wid=0):
    if wid==0:
        return r'$[\mu_{R},\mu_{F}]=[1.0,1.0]$'
    elif wid==1:
        return r'$[\mu_{R},\mu_{F}]=[1.0,2.0]$'
    elif wid==2:
        return r'$[\mu_{R},\mu_{F}]=[1.0,0.5]$'
    elif wid==3:
        return r'$[\mu_{R},\mu_{F}]=[2.0,1.0]$'
    elif wid==4:
        return r'$[\mu_{R},\mu_{F}]=[2.0,2.0]$'
    elif wid==5:
        return r'$[\mu_{R},\mu_{F}]=[2.0,0.5]$'
    elif wid==6:
        return r'$[\mu_{R},\mu_{F}]=[0.5,1.0]$'
    elif wid==7:
        return r'$[\mu_{R},\mu_{F}]=[0.5,2.0]$'
    elif wid==8:
        return r'$[\mu_{R},\mu_{F}]=[0.5,0.5]$'
    else:
        return r'NNPDF replica '+str(wid-9)


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

# Boost lab --> CS using a parametric representation of the boost matrix
# flip_z = the direction of the z-axis in the CS frame is determined by the sign of pz in the lab
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

# Boost lab --> CS using the TLorentzVector class
# flip_z = the direction of the z-axis in the CS frame is determined by the sign of pz in the lab
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
        variables[name] = array( 'f', [ 0.0 ] ) 
        if tree!=None:
            tree.Branch(name, variables[name], name+'/F')
    variables['weights'] = array( 'f', 109*[ 0.0 ] )
    if tree!=None:
        tree.Branch('weights', variables['weights'], 'weights[109]/F')

    return variables
    
def fill_default(variables):
    for key,var in variables.items():
        var[0] = 0.0

# Create 2D maps for (qT,y) filled with the mean of the test functions
# The binning is taken from global variables 
def add_histo2D(charges=['Wminus','Wplus'], var=['Wdress'], coeff=['A0','A1','A2','A3','A4','A5', 'A6', 'A7', 'A8'], weights=[0]):
    # binning
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

# Create 2D maps for the (cos*,phi*) distribution in bins of (qT,y)
# The binning is decided inside
def add_histo2D_CS(charges=['Wminus','Wplus'], var=['Wdress'], coeff_eval=['fit','val']):
    # binning
    bins_cos = array( 'f',  np.linspace(-1.0, 1.0, 21) )
    bins_phi = array( 'f',  np.linspace(0.0, 2*math.pi, 21) )    
    histos = {}
    # create TH2D
    for q in charges:
        histos[q] = {}
        for v in var:
            histos[q][v] = {}
            for c in coeff_eval:
                histos[q][v][c] = {}
                for iqt in range(np_bins_qt.size-1):
                    for iy in range((np_bins_y.size-1)/2, np_bins_y.size-1):
                        bin_y = 'y{:03.2f}'.format(np_bins_y[iy])+'_'+'y{:03.2f}'.format(np_bins_y[iy+1])
                        bin_qt = 'qt{:03.1f}'.format(np_bins_qt[iqt])+'_'+'qt{:03.1f}'.format(np_bins_qt[iqt+1])
                        name = q+'_'+v+'_'+c+'_'+bin_y+'_'+bin_qt
                        histos[q][v][c][name] = (ROOT.TH2F(name, name, len(bins_cos)-1, bins_cos, len(bins_phi)-1, bins_phi), 
                                               ROOT.TH2F(name+'_norm', name+'_norm', len(bins_cos)-1, bins_cos, len(bins_phi)-1, bins_phi) )
                        histos[q][v][c][name][0].Sumw2()
                        histos[q][v][c][name][1].Sumw2()
        
    return histos

def add_histo2D_lepton(charges=['Wminus','Wplus'], var=['Wdress'], coeff_eval=['val'], masses=[80.000], coeff=['A0'], plot_vars=['eta','pt']):
    # binning
    np_bins_xx = np_bins_eta
    np_bins_yy = np_bins_pt
    if plot_vars[1]=='1/pt':
        np_bins_yy = np_bins_ptinv
        
    bins_yy = array( 'f',  np_bins_yy )
    bins_xx = array( 'f',  np_bins_xx )

    np_bins_y_extL = np.insert(np_bins_y,   0, [-10.])
    np_bins_y_ext  = np.append(np_bins_y_extL, [+10.])
    np_bins_qt_ext = np.append(np_bins_qt, 999.)

    nbins_y = np_bins_y_ext.size - 1 
    nbins_qt = np_bins_qt_ext.size - 1 

    histos = {}

    # create TH2D
    for q in charges:
        histos[q] = {}
        for v in var:
            histos[q][v] = {}
            for ceval in coeff_eval:
                histos[q][v][ceval] = {}
                for m in masses:
                    mass_str = 'M'+'{:05.3f}'.format(m)
                    histos[q][v][ceval][mass_str] = {}
                    for qt in range(1, nbins_qt+1):
                        qt_bin = 'qt{:03.1f}'.format(np_bins_qt_ext[qt-1])+'_'+'qt{:03.1f}'.format(np_bins_qt_ext[qt]) if qt<nbins_qt else 'OF'
                        for y in range(nbins_y/2+1, nbins_y+1):
                            y_bin = 'y{:03.2f}'.format(np_bins_y_ext[y-1])+'_'+'y{:03.2f}'.format(np_bins_y_ext[y]) if y<nbins_y else 'OF'
                            bin_name = qt_bin+'_'+y_bin
                            histos[q][v][ceval][mass_str][bin_name] = {}                                
                            for c in (coeff+['', 'UL']):
                                name = q+'_'+v+'_'+ceval+'_'+mass_str+'_'+bin_name+'_'+c
                                h2 = ROOT.TH2F(name, name, len(bins_xx)-1, bins_xx, len(bins_yy)-1, bins_yy) 
                                h2.Sumw2()
                                histos[q][v][ceval][mass_str][bin_name][c]= h2
    return histos

def add_histo2D_lepton_mc_weights(charges=['Wminus','Wplus'],  masses=[80.000], weights=[], plot_vars=['eta','pt']):
    # binning
    np_bins_xx = np_bins_eta
    np_bins_yy = np_bins_pt
    if plot_vars[1]=='1/pt':
        np_bins_yy = np_bins_ptinv
        
    bins_yy = array( 'f',  np_bins_yy )
    bins_xx = array( 'f',  np_bins_xx )

    np_bins_y_extL = np.insert(np_bins_y,   0, [-10.])
    np_bins_y_ext  = np.append(np_bins_y_extL, [+10.])
    np_bins_qt_ext = np.append(np_bins_qt, 999.)

    nbins_y = np_bins_y_ext.size - 1 
    nbins_qt = np_bins_qt_ext.size - 1 

    histos = {}

    # create TH2D
    for q in charges:
        histos[q] = {}
        for m in masses:
            mass_str = 'M'+'{:05.3f}'.format(m)
            histos[q][mass_str] = {}
            for qt in range(1, nbins_qt+1):
                qt_bin = 'qt{:03.1f}'.format(np_bins_qt_ext[qt-1])+'_'+'qt{:03.1f}'.format(np_bins_qt_ext[qt]) if qt<nbins_qt else 'OF'
                for y in range(nbins_y/2+1, nbins_y+1):
                    y_bin = 'y{:03.2f}'.format(np_bins_y_ext[y-1])+'_'+'y{:03.2f}'.format(np_bins_y_ext[y]) if y<nbins_y else 'OF'
                    bin_name = qt_bin+'_'+y_bin
                    histos[q][mass_str][bin_name] = {}                                
                    for w in weights:
                        name = q+'_'+mass_str+'_'+bin_name+'_'+str(w)
                        h2 = ROOT.TH2F(name, name, len(bins_xx)-1, bins_xx, len(bins_yy)-1, bins_yy) 
                        h2.Sumw2()
                        histos[q][mass_str][bin_name][str(w)]= h2
    return histos

def add_histo2D_lepton_pt_scale(charges=['Wminus','Wplus'],  masses=[80.000], pt_scales=[], plot_vars=['eta','pt']):
    # binning
    np_bins_xx = np_bins_eta
    np_bins_yy = np_bins_pt
    if plot_vars[1]=='1/pt':
        np_bins_yy = np_bins_ptinv
        
    bins_yy = array( 'f',  np_bins_yy )
    bins_xx = array( 'f',  np_bins_xx )

    np_bins_y_extL = np.insert(np_bins_y,   0, [-10.])
    np_bins_y_ext  = np.append(np_bins_y_extL, [+10.])
    np_bins_qt_ext = np.append(np_bins_qt, 999.)

    nbins_y = np_bins_y_ext.size - 1 
    nbins_qt = np_bins_qt_ext.size - 1 

    histos = {}

    # create TH2D
    for q in charges:
        histos[q] = {}
        for m in masses:
            mass_str = 'M'+'{:05.3f}'.format(m)
            histos[q][mass_str] = {}
            for qt in range(1, nbins_qt+1):
                qt_bin = 'qt{:03.1f}'.format(np_bins_qt_ext[qt-1])+'_'+'qt{:03.1f}'.format(np_bins_qt_ext[qt]) if qt<nbins_qt else 'OF'
                for y in range(nbins_y/2+1, nbins_y+1):
                    y_bin = 'y{:03.2f}'.format(np_bins_y_ext[y-1])+'_'+'y{:03.2f}'.format(np_bins_y_ext[y]) if y<nbins_y else 'OF'
                    bin_name = qt_bin+'_'+y_bin
                    histos[q][mass_str][bin_name] = {}                                
                    for s in pt_scales:
                        name = q+'_'+mass_str+'_'+bin_name+'_'+'{:2.1f}'.format(s)
                        h2 = ROOT.TH2F(name, name, len(bins_xx)-1, bins_xx, len(bins_yy)-1, bins_yy) 
                        h2.Sumw2()
                        histos[q][mass_str][bin_name]['{:2.1f}'.format(s)]= h2
    return histos



# The test function for projecting-out one harmonic at the time
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
    elif coeff=='A8':
        val = 14.-10.*(1+ps[0]*ps[0])
    return val

# Fill a 2D map with the value of the test function needed to project-out one harmonic at the time.
def fill_coefficients(histos={}, q='', coefficients_for_histos=['A0'], var='Wdress', weight_name=0, ps_W=(), ps_CS=(), weight=1.0):
    for coeff in coefficients_for_histos:
        (h,h_norm) = histos[q][var][coeff][str(weight_name)]
        h.Fill(ps_W[0],ps_W[1], weight*test_A(coeff=coeff,ps=ps_CS) )
        h_norm.Fill(ps_W[0],ps_W[1], weight )

# mass reweighting
def mass_weight(mass=80.000, mass_target=80.419, mass_mc=80.419, width_target=2.047, width_mc=2.047):
    weight = (math.pow(mass*mass - mass_mc*mass_mc, 2.0) + math.pow(mass_mc*width_mc, 2.0))/(math.pow(mass*mass - mass_target*mass_target, 2.0) + math.pow(mass_target*width_target, 2.0))
    #print 'mass=', mass, ' => ', weight
    return weight

# Evaluate the angular pdf given a PS point (cos*,phi*) and the value of the 8 parameters
def angular_pdf(ps=(), coeff_vals=[], verbose=False):
    (x,y) = ps
    UL = (1.0 + x*x)
    L = 0.5*(1-3*x*x)
    T = 2.0*x*math.sqrt(1-x*x)*math.cos(y)
    I = 0.5*(1-x*x)*math.cos(2*y)
    A = math.sqrt(1-x*x)*math.cos(y)
    P = x
    p7 = (1-x*x)*math.sin(2*y)
    p8 = 2.0*x*math.sqrt(1-x*x)*math.sin(y)
    p9 = math.sqrt(1-x*x)*math.sin(y)
    if verbose:
        print ('\t pdf = 3./16./math.pi * ( UL + %s*L + %s*T + %s*I + %s*A + %s*P + %s*p7 + %s*p8 + %s*p9)' % (coeff_vals[0], coeff_vals[1], coeff_vals[2], coeff_vals[3], coeff_vals[4], coeff_vals[5], coeff_vals[6], coeff_vals[7]) )
    return 3./16./math.pi * ( UL + coeff_vals[0]*L + coeff_vals[1]*T + coeff_vals[2]*I + coeff_vals[3]*A + coeff_vals[4]*P + coeff_vals[5]*p7 + coeff_vals[6]*p8 + coeff_vals[7]*p9)

def angular_pdf_string(coeff_vals=[]):
    title = r'$\frac{3}{16\pi}['
    UL = r'(1 + \cos^2\theta)'
    pols = [r'\frac{1}{2}(1-3\cos^2\theta)',   
            r'\sin2\theta \cos\phi',
            r'\frac{1}{2}\sin^2\theta\cos2\phi',
            r'\sin\theta\cos\phi',
            r'\cos\theta',
            r'\sin^2\theta\sin2\phi',
            r'\sin2\theta\sin\phi',
            r'\sin\theta\sin\phi'
            ]
    title += UL
    for ic,c in enumerate(coeff_vals):
        if c!=0.0:
            title += ('+'+pols[ic])
    title += r']$'
    #print title
    return title

# Evaluate the cumulative angular pdf given a PS point (cos*,phi*) and the value of the 8 parameters
def cumulative_angular_pdf(bin_cos=(), bin_phi=(), coeff_vals=[], verbose=False):
    (a,b) = bin_cos 
    (c,d) = bin_phi     
    
    P1 = lambda x : math.pow(x,3.)/3. + x
    UL = (P1(b)-P1(a))*(d-c)

    P2 = lambda x : x-math.pow(x,3.)
    L = 0.5*(P2(b)-P2(a))*(d-c)

    P3 = lambda x : -2./3.*math.pow(1-x*x, 3./2.)
    P4 = lambda x : math.sin(x)
    T = (P3(b)-P3(a))*(P4(d)-P4(c))

    P5 = lambda x : x - math.pow(x,3.)/3.
    P6 = lambda x : 0.5*math.sin(2*x)
    I = 0.5*(P5(b)-P5(a))*(P6(d)-P6(c))

    P7 = lambda x : 0.5*(math.sqrt(1-x*x)*x + math.asin(x))
    A =  (P7(b)-P7(a))*(P4(d)-P4(c))

    P8 = lambda x : x*x/2.
    P =  (P8(b)-P8(a))*(d-c)
    
    P9 = lambda x : -0.5*math.cos(2*x)
    p7 = (P2(b)-P2(a))*(P9(d)-P9(c))

    P10 = lambda x : -math.cos(x)
    p8 = (P3(b)-P3(a))*(P10(d)-P10(c))

    p9 = (P7(b)-P7(a))*(P10(d)-P10(c))

    if verbose:
        print ('\t pdf = 3./16./math.pi * ( UL + %s*L + %s*T + %s*I + %s*A + %s*P + %s*p7 + %s*p8 + %s*p9)' % (coeff_vals[0], coeff_vals[1], coeff_vals[2], coeff_vals[3], coeff_vals[4], coeff_vals[5], coeff_vals[6], coeff_vals[7]) )
    return 3./16./math.pi * ( UL + coeff_vals[0]*L + coeff_vals[1]*T + coeff_vals[2]*I + coeff_vals[3]*A + coeff_vals[4]*P + coeff_vals[5]*p7 + coeff_vals[6]*p8 + coeff_vals[7]*p9)
    #return (b-a)*(d-c)/(4.*math.pi)

# read coefficients from res
def get_coeff_vals(res={}, coeff_eval='fit', bin_y='', qt=0.0, coeff=['A0'], np_bins_template_qt_new=np.array([])):

    if np_bins_template_qt_new.size==0:
        np_bins_template_qt = np_bins_qt
    else:
        np_bins_template_qt = np_bins_template_qt_new

    coeff_vals = np.zeros(8)
    for ic,c in enumerate(coeff):
        coeff_val = 0.0
        if coeff_eval == 'fit':
            order = len(res[c+'_'+bin_y+'_fit'])
            for o in range(order):
                coeff_val += math.pow(qt,o)*res[c+'_'+bin_y+'_fit'][o]
        elif coeff_eval == 'val':
            iqt = np.where(np_bins_template_qt<=qt)[0][-1]
            coeff_val = res[c+'_'+bin_y+'_val'][iqt]
        coeff_vals[ic] = coeff_val
    return coeff_vals

# Determine the angular_pdf in a y-bin as a function of qt (read results 'res' from external file)
# Two modes:
#  - 'fit' : read the coefficients of a polynomial fit to A(qT). Rebuild the polynomnial from it.
#  - 'val' : read the bin-by-bin value of A
# Negative values for the pdf are possible, but expected only for qT outside the fit range in 'fit' mode 
def weight_coeff(res={}, coeff_eval='fit', bin_y='', qt=0.0, ps=(), coeff=['A0'], verbose=False):
    coeff_vals = get_coeff_vals(res=res, coeff_eval=coeff_eval, bin_y=bin_y, qt=qt, coeff=coeff)
    val = angular_pdf(ps=ps, coeff_vals=coeff_vals)
    if val > 0.0:        
        if verbose:
            print 'Fit type', coeff_eval, ': bin', bin_y, 'at qt =', qt, 'for ps=', ps, 'yields val = ', val
        return val
    else:
        if verbose:
            print ('Fit type', coeff_eval, ': bin', bin_y, 'at qt =', qt, 'for ps=', ps, 'yields pdf = ', val, '<=0.0. Return 1.0')
        return 1.0

# Determine the cumulative angular_pdf in a y-bin as a function of qt (read results 'res' from external file)
# Two modes:
#  - 'fit' : read the coefficients of a polynomial fit to A(qT). Rebuild the polynomnial from it.
#  - 'val' : read the bin-by-bin value of A
# Negative values for the pdf are possible, but expected only for qT outside the fit range in 'fit' mode 
def weight_coeff_cumulative(res={}, coeff_eval='fit', bin_y='', qt=0.0, bin_cos=(), bin_phi=(), coeff=['A0'], verbose=False):
    coeff_vals = get_coeff_vals(res=res, coeff_eval=coeff_eval, bin_y=bin_y, qt=qt, coeff=coeff)
    val = cumulative_angular_pdf(bin_cos=bin_cos, bin_phi=bin_phi, coeff_vals=coeff_vals)
    if val > 0.0:        
        if verbose:
            print 'Fit type', coeff_eval, ': bin', bin_y, 'at qt =', qt, 'for ps=', ps, 'yields val = ', val
        return val
    else:
        if verbose:
            print ('Fit type', coeff_eval, ': bin', bin_y, 'at qt =', qt, 'for ps=', ps, 'yields pdf = ', val, '<=0.0. Return 1.0')
        return 1.0


# Given a ps_W point, find the bin in np_bins_y/qt. Return if the point is not within the bins
# Fill a 2D map with ps_CS with weight 1/weight_coeff
def fill_weighted_CS(res={}, histos={}, q='', var='Wdress', coeff_eval=['fit'], ps_W=(), ps_CS=(), weight=1.0, coeff=['A0', 'A1', 'A2', 'A3', 'A4'], verbose=False):

    bin = find_y_qt_bin( ps=ps_W, verbose=verbose)
    if bin[0]=='OF' or bin[2]=='OF':
        return

    (bin_y, iy_low, bin_qt, iqt_low) = bin 

    for ceval in coeff_eval:
        name = q+'_'+var+'_'+ceval+'_'+bin_y+'_'+bin_qt
        (h,h_norm) = histos[q][var][ceval][name]
        h.Fill(ps_CS[0], ps_CS[1], weight/weight_coeff(res=res, coeff_eval=ceval, bin_y=bin_y, qt=ps_W[1], ps=ps_CS, coeff=coeff) )
        h_norm.Fill(ps_CS[0], ps_CS[1], weight)

# fill lepton kinematics (pt,eta) in the lab
def fill_lepton_lab(res={}, histos={}, q='', var='WpreFSR', coeff_eval=['val'], masses=[80.419], ps_W=(), ps_CS=(), ps_lep=(), weight=1.0, coeff=[]):

    # if (qt,y) not in the bins, use MC
    useMC = False
    (bin_y, iy_low, bin_qt, iqt_low) = find_y_qt_bin( ps=ps_W, verbose=False )
    if bin_y=='OF' or bin_qt=='OF':
        useMC = True

    # loop over evaluation methods ['val', 'fit']
    for ceval in coeff_eval:
        # the total pdf ~ [(1+cos^2) + Sum Ai*Pi]
        pdf = weight_coeff(res=res, coeff_eval=ceval, bin_y=bin_y, qt=ps_W[1], ps=ps_CS, coeff=coeff) if (not useMC) else 1.0
        if pdf<=0.0:
            print ('Fit type', coeff_eval, ': bin', bin_y, 'at qt =', qt, 'for ps=', ps_CS, 'yields pdf = ', pdf, '<=0.0. Return 1.0')
            return
        # loop over MW hypotheses
        for m in masses:
            wm = mass_weight(mass=ps_W[2], mass_target=m )
            mass_str = 'M'+'{:05.3f}'.format(m)
            # '' => MC, UL => (1+cos^2), others : [0,0,...,1,...,0,0]
            for c in (coeff+['', 'UL']):
                wc = 1.0
                if c=='':
                    wc = 1.0
                else:
                    coeff_vals = np.zeros(8)
                    if c!='UL':
                        coeff_vals[int(c[1])] = 1.0
                    wc = angular_pdf(ps=ps_CS, coeff_vals=coeff_vals)/pdf if (not useMC) else 1.0
                #print c, 'weight=', wc
                h2 = histos[q][var][ceval][mass_str][bin_qt+'_'+bin_y][c]
                h2.Fill( ps_lep[0], ps_lep[1], weight*wc*wm)

# fill lepton kinematics (pt,eta) in the lab
def fill_lepton_lab_mc_weights(histos={}, q='', masses=[80.419], ps_W=(), ps_CS=(), ps_lep=(), weight=1.0, weight_id=0):

    # if (qt,y) not in the bins, use MC
    useMC = False
    (bin_y, iy_low, bin_qt, iqt_low) = find_y_qt_bin( ps=ps_W, verbose=False )
    if bin_y=='OF' or bin_qt=='OF':
        useMC = True

    # loop over MW hypotheses
    for m in masses:
        wm = mass_weight(mass=ps_W[2], mass_target=m )
        mass_str = 'M'+'{:05.3f}'.format(m)
        h2 = histos[q][mass_str][bin_qt+'_'+bin_y][str(weight_id)]
        h2.Fill( ps_lep[0], ps_lep[1], weight*wm)

# fill lepton kinematics (pt,eta) in the lab
def fill_lepton_lab_pt_scale(histos={}, q='', masses=[80.419], ps_W=(), ps_CS=(), ps_lep=(), weight=1.0, pt_scale=1.0):

    # if (qt,y) not in the bins, use MC
    useMC = False
    (bin_y, iy_low, bin_qt, iqt_low) = find_y_qt_bin( ps=ps_W, verbose=False )
    if bin_y=='OF' or bin_qt=='OF':
        useMC = True

    # loop over MW hypotheses
    for m in masses:
        wm = mass_weight(mass=ps_W[2], mass_target=m )
        mass_str = 'M'+'{:05.3f}'.format(m)
        h2 = histos[q][mass_str][bin_qt+'_'+bin_y]['{:2.1f}'.format(pt_scale)]
        h2.Fill( ps_lep[0], ps_lep[1]*(1+pt_scale*1e-04), weight*wm)


# compute dsigma/dphidcos in the CS frame
def angular_pdf_CS(x, y, coeff=[]):
    UL = (1.0 + x*x)
    L = 0.5*(1-3*x*x)
    T = 2.0*x*np.sqrt(1-x*x)*np.cos(y)
    I = 0.5*(1-x*x)*np.cos(2*y)
    A = np.sqrt(1-x*x)*np.cos(y)
    P = x
    p7 = (1-x*x)*np.sin(2*y)
    p8 = 2.0*x*np.sqrt(1-x*x)*np.sin(y)
    p9 = np.sqrt(1-x*x)*np.sin(y)
    return 3./16./math.pi * ( UL + coeff[0]*L + coeff[1]*T + coeff[2]*I + coeff[3]*A + coeff[4]*P + coeff[5]*p7 + coeff[6]*p8 + coeff[7]*p9)
    #return (1.0 + 0.0*x)/(4.0*math.pi)

# fill (cos,phi) grid in CS frame
def make_grid_CS(res={}, coeff_eval='fit', bin_y='', qt=0.0, h=None, coeff=[], ntoys=1000):

    coeff_vals = get_coeff_vals(res=res, coeff_eval=coeff_eval, bin_y=bin_y, qt=qt, coeff=coeff)
    bins_cos = np.linspace(h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax()-h.GetXaxis().GetBinWidth(1), h.GetNbinsX())
    bins_phi = np.linspace(h.GetYaxis().GetXmin(), h.GetYaxis().GetXmax()-h.GetYaxis().GetBinWidth(1), h.GetNbinsY())
    xx, yy = np.meshgrid( bins_cos, bins_phi )        

    pdf_evals = np.zeros( bins_cos.size*bins_phi.size)
    #pdf_evals = angular_pdf_CS(xx,yy, coeff=coeff_vals).flatten()    
    for ix in range(10):
        for iy in range(10):
            delta_xx = ix*h.GetXaxis().GetBinWidth(1)/10.
            delta_yy = iy*h.GetYaxis().GetBinWidth(1)/10.
            pdf_evals += angular_pdf_CS(xx+delta_xx,yy+delta_yy, coeff=coeff_vals).flatten()    

    # prevent negative pdf values 
    pdf_evals = np.absolute(pdf_evals)

    pdf_evals_norm = np.sum(pdf_evals)    
    pdf_evals *= (1./pdf_evals_norm if pdf_evals_norm>0. else 1.0)        

    rnd = np.random.choice(  np.arange(0, len(bins_cos)*len(bins_phi), 1), size=ntoys, p=pdf_evals )

    for i,idx in enumerate(rnd):
        idx_cos =  idx%len(bins_cos)+1
        idx_phi =  idx/len(bins_cos)+1
        h.Fill( h.GetXaxis().GetBinCenter(idx_cos),  h.GetYaxis().GetBinCenter(idx_phi) )

###########################
