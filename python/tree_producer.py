from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )
from ROOT import TLorentzVector

import math
import numpy as np 

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Handle, Events

verbose = False

def isMuon(p):
    if not (p.isPromptFinalState() and abs(p.pdgId())==13 and p.pt()>20. and abs(p.eta())<2.4):
        return False
    mother = p
    while(mother.numberOfMothers()>0):
        if abs(mother.pdgId())==24:
            return True
        mother = mother.mother(0)
    return False

def isNeutrino(p):
    if not (p.isPromptFinalState() and abs(p.pdgId())==14):
        return False
    mother = p
    while(mother.numberOfMothers()>0):
        if abs(mother.pdgId())==24:
            return True
        mother = mother.mother(0)
    return False

def isPhoton(p):
    return (p.isPromptFinalState() and p.pdgId()==22 and p.pt()>0.0 and abs(p.eta())<2.5)

def deltaR(a,b):
    return math.sqrt( math.pow(a.eta()-b.eta(),2) + math.pow( math.acos( math.cos(a.phi()-b.phi())),2) )

def printp(tag, p, other):
    print tag+' ['+str(p.pdgId())+']: ('+'{:04.3f}'.format(p.pt())+',{:04.3f}'.format(p.eta())+',{:04.3f}'.format(p.phi())+') .... '+other

outfile = ROOT.TFile("tree.root", "RECREATE")
outtree = ROOT.TTree('tree', 'tree')

lhe_mW_     = np.zeros(1, dtype=float)
preFSR_mW_  = np.zeros(1, dtype=float)
recFSR_mW_  = np.zeros(1, dtype=float)
postFSR_mW_ = np.zeros(1, dtype=float)

outtree.Branch('lhe_mW', lhe_mW_, 'lhe_mW/D')
outtree.Branch('preFSR_mW',  preFSR_mW_,  'preFSR_mW/D')
outtree.Branch('recFSR_mW',  recFSR_mW_,  'recFSR_mW/D')
outtree.Branch('postFSR_mW', postFSR_mW_, 'postFSR_mW/D')


print "Opening file..."
events = Events('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0AF0207B-EFBE-E611-B4BE-0CC47A7FC858.root')
print "File opened..."


genH, genN = Handle("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
infoH, infoN = Handle("GenEventInfoProduct"), "generator"
lheH, lheN = Handle("LHEEventProduct"), "externalLHEProducer"

for i,event in enumerate(events):    

    if i%100==0:
        print "Event", i

    event.getByLabel(genN,genH)
    event.getByLabel(infoN,infoH)
    event.getByLabel(lheN,lheH)

    generator = infoH.product()
    weights = list(generator.weights())

    lhe = lheH.product()
    hepeup = lhe.hepeup()
    Wp4_lhe = [0.,0.,0.,0.,0.]
    for p in range(hepeup.NUP):
        if abs(hepeup.IDUP[p])==24:
            for x in range(5):
                Wp4_lhe[x] = hepeup.PUP[p][x]

    genParticles = list(genH.product())
    muons = [p for p in genParticles if isMuon(p)] 
    neutrinos = [p for p in genParticles if isNeutrino(p)] 
    gammas = [p for p in genParticles if isPhoton(p)] 
    
    if len(muons)==0:
        continue
    if len(neutrinos)==0:
        print "No muon neutrinos for event", i
        continue

    print '***********'
    muons.sort(reverse=True)
    mu = muons[0]
    if verbose:
        printp('muon', mu, '')

    neutrinos.sort(reverse=True)
    nu = neutrinos[0]
    if verbose:
        printp('neut', nu, '')

    mother = mu
    mu_prefsr = mu
    while(mother.numberOfMothers()>0):
        if verbose:
            printp('MOTH', mother, '')
        if abs(mother.pdgId()==13) and mother.statusFlags().isLastCopyBeforeFSR():
            mu_prefsr = mother
        mother = mother.mother(0)

    fsr = []
    gammas.sort(reverse=True)
    for ng,g in enumerate(gammas):
        dR = deltaR(mu,g)
        if dR<0.1:
            fsr.append(g)
            if verbose:
                printp('>gam', g, 'dR:{:03.2f}'.format(dR))

    nup4 = nu.p4() 
    mup4 = mu.p4() 
    mup4_prefsr = mu_prefsr.p4()
    mup4_recfsr = mup4
    for g in fsr:
        mup4_recfsr += g.p4()
        
    mass_lhe = Wp4_lhe[4]
    mass = (mup4+nup4).mass()
    mass_recfsr = (mup4_recfsr+nup4).mass()
    mass_prefsr = (mup4_prefsr+nup4).mass()    
    if verbose:
        print 'LHE:        '+'{:04.3f}'.format(mass_lhe)
        print 'pre-FSR mu: '+'{:04.3f}'.format(mass_prefsr)
        print 'dressed mu: '+'{:04.3f}'.format(mass_recfsr)
        print 'post-FSR  : '+'{:04.3f}'.format(mass)
    lhe_mW_[0] = mass_lhe
    preFSR_mW_[0] = mass_prefsr
    recFSR_mW_[0] = mass_recfsr
    postFSR_mW_[0] = mass

    outtree.Fill()

outfile.cd()
outtree.Write("tree", ROOT.TObject.kOverwrite)
outfile.Close()
