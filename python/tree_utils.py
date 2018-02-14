import copy
import math
import numpy as np 
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
    if debug:
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
        tree.Branch(name, variables[name], name+'/D')

    variables['weights'] = np.zeros(109, dtype=float)
    tree.Branch('weights', variables['weights'], name+'[109]/D')

    return variables
    

def fill_default(variables):
    for key,var in variables.items():
        var[0] = 0.0
