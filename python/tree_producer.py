from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )
from ROOT import TLorentzVector

import time
import math
import os

from tree_utils import *

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Handle, Events


DY      = 'NC'
verbose = False
debug   = True

# open output file
outfile = None
if debug:
    outfile = ROOT.TFile(os.environ['CMSSW_BASE']+'/src/Wmass/test/'+'tree_'+argv[1]+'.root', "RECREATE")
else:
    outfile = ROOT.TFile(os.environ['CMSSW_BASE']+'/src/Wmass/test/'+'tree.root', "RECREATE")
outtree = ROOT.TTree('tree', 'tree')

# add branches to tree
variables = add_vars(outtree)

print "Opening file..."
filename = []
if DY == 'NC':
    filename.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root')
else:
    filename.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0AF0207B-EFBE-E611-B4BE-0CC47A7FC858.root')
    #filename.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/002F2CE1-38BB-E611-AF9F-0242AC130005.root')
    #filename.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0A85AA82-45BB-E611-8ACD-001E674FB063.root')
    #filename.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0CA050B2-57BB-E611-8A7A-001E674FBA1D.root')

events = Events(filename)
print "File opened.... Tot. num of events:", events.size()

# handles
genH, genN = Handle("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
infoH, infoN = Handle("GenEventInfoProduct"), "generator"
lheH, lheN = Handle("LHEEventProduct"), "externalLHEProducer"

start = time.time()
# loop over the events
for i,event in enumerate(events):    

    if i%1000==0:
        print "Processing event", i, '/', events.size()

    if debug:
        if i%100==0:
            print "Processing event", i, '/', events.size()
        if i>2000:
            break

    # fill with default values
    fill_default(variables)
    
    # load LHE product
    event.getByLabel(lheN,lheH)
    lhe = lheH.product()
    hepeup = lhe.hepeup()

    isW = False
    isWToMuNu = False

    Wp4_lhe = [0.,0.,0.,0.,0.]
    for p in range(hepeup.NUP):
        if abs(hepeup.IDUP[p])==13:
            isWToMuNu = True
        if abs(hepeup.IDUP[p])== (24 if DY=='CC' else 23):
            isW = True
            for x in range(5):                
                Wp4_lhe[x] = hepeup.PUP[p][x]
            Wp4 = ROOT.TLorentzVector(Wp4_lhe[0], Wp4_lhe[1], Wp4_lhe[2], Wp4_lhe[3])     
            variables['lhe_y'][0]   =  Wp4.Rapidity()
            variables['lhe_qt'][0]  =  Wp4.Pt()
            variables['lhe_mass'][0]=  Wp4.M()
            variables['lhe_phi'][0] =  Wp4.Phi()
   
    # consider only W->munu events
    if not (isW and isWToMuNu):
        continue

    # LHE weights
    norm = abs(lhe.weights()[0].wgt)
    variables['weight'][0]   = lhe.weights()[0].wgt/norm # nominal (+/- 1)
    variables['weightNU'][0] = lhe.weights()[1].wgt/norm # muR *= 1 ,  muF *= 2 
    variables['weightND'][0] = lhe.weights()[2].wgt/norm # muR *= 1 ,  muF *= 0.5
    variables['weightUN'][0] = lhe.weights()[3].wgt/norm # muR *= 2 ,  muF *= 1
    variables['weightUU'][0] = lhe.weights()[4].wgt/norm # muR *= 2 ,  muF *= 2
    variables['weightDN'][0] = lhe.weights()[6].wgt/norm # muR *= 0.5, muF *= 1
    variables['weightDD'][0] = lhe.weights()[8].wgt/norm # muR *= 0.5, muF *= 0.5
   
    if verbose:
        print "Weights:"
        for w in ['', 'NU', 'ND', 'UN', 'UU', 'DN', 'DD']:
            print w, variables['weight'+w][0]

    # read gen particles
    event.getByLabel(genN,genH)
    genParticles = list(genH.product())

    # filter muons, neutrinos and gammas
    muons = [p for p in genParticles if isMuon(p)] 
    neutrinos = [p for p in genParticles if isNeutrino(p)] 
    gammas = [p for p in genParticles if isPhoton(p)] 

    # sort by pt
    muons.sort(reverse=True)
    neutrinos.sort(reverse=True)
    
    # CC: consider events with 0 muons to study acceptance 
    if DY=='CC' and len(muons)==0:
        if verbose:
            print "No muon in W>munu for event", i, ". Try to understand why:"
            print " > W rapidity:    ", variables['lhe_y'][0]
        variables['muLost'][0] += 1.0
        if len(neutrinos)>0:
            variables[t+'_nu_pt'][0]  = neutrinos[0].p4().Pt()
            variables[t+'_nu_eta'][0] = neutrinos[0].p4().Eta()
            variables[t+'_nu_phi'][0] = neutrinos[0].p4().Phi()

    # consider events with 0 neutrinos to study acceptance 
    if DY=='CC' and len(neutrinos)==0:
        if verbose:
            print "No muon neutrinos for event", i, ". Try to understand why:"
            print " > W rapidity:    ", variables['lhe_y'][0]
            print " > Muon rapidity: ", muons[0].p4().Rapidity()
        variables['nuLost'][0] += 1.0
        if len(muons)>0:
            variables[t+'_mu_pt'][0]  = muons[0].p4().Pt()
            variables[t+'_mu_eta'][0] = muons[0].p4().Eta()
            variables[t+'_mu_phi'][0] = muons[0].p4().Phi()

    # NC: consider events with 1 muons to study acceptance 
    if DY=='NC' and len(muons)<2:
        if verbose:
            print "No muons in Z>mumu for event", i, ". Try to understand why:"
            print " > Z rapidity:    ", variables['lhe_y'][0]
        variables['muLost'][0] += 1.0
        if len(muons)==1:
            variables[t+'_nu_pt'][0]  = muons[0].p4().Pt()
            variables[t+'_nu_eta'][0] = muons[0].p4().Eta()
            variables[t+'_nu_phi'][0] = muons[0].p4().Phi()

    # if no muons or no neutrinos, save the event and continue
    if DY=='CC':
        if (len(muons)==0 or len(neutrinos)==0):
            outtree.Fill()
            continue
    elif DY=='NC':
        if len(muons)<2:
            outtree.Fill()
            continue
        
    # the muon is the first ranked by pt if CC else the mu-
    mu = muons[0]
    if DY=='NC':
        mu = (muons[0] if muons[0].pdgId()==-13 else muons[1])
    variables['isFromW'][0] += int(isFromW(mu) if DY=='CC' else isFromZ(mu))
    variables['mu_charge'][0] = mu.pdgId()
    if verbose:
        printp('muon', mu, '')

    # the neutrino is the first ranked by pt if CC else the mu+
    nu = None
    if DY=='CC':
        nu = neutrinos[0]
    elif DY=='NC':
        nu = (muons[0] if muons[0].pdgId()==+13 else muons[1])
    variables['isFromW'][0] += int(isFromW(nu) if DY=='CC' else isFromZ(nu))
    variables['nu_charge'][0] = nu.pdgId()
    if verbose:
        printp('neut', nu, '')

    # pre-FSR
    mother = mu
    mu_prefsr = mu
    while(mother.numberOfMothers()>0):
        if verbose:
            printp('MOTH', mother, '')
        if abs(mother.pdgId())==13 and mother.statusFlags().isLastCopyBeforeFSR():
            mu_prefsr = mother
        mother = mother.mother(0)

    mother = nu
    nu_prefsr = nu
    while(mother.numberOfMothers()>0):
        if verbose:
            printp('MOTH', mother, '')
        if abs(mother.pdgId())==13 and mother.statusFlags().isLastCopyBeforeFSR():
            nu_prefsr = mother
        mother = mother.mother(0)

    # standard dressing alorithm
    mu_fsr,nu_fsr = [],[]
    gammas.sort(reverse=True)
    for ng,g in enumerate(gammas):
        dR_mu = deltaR(mu,g)
        dR_nu = deltaR(nu,g)
        dR = (dR_mu if DY=='CC' else min(dR_mu,dR_nu))
        if dR<0.1:
            if DY=='CC':
                mu_fsr.append(g)
            elif DY=='NC':
                if dR_mu<dR_nu:
                    mu_fsr.append(g)
                else:
                    nu_fsr.append(g)
            if verbose:
                printp('>gam', g, 'dR:{:03.2f}'.format(dR))

    # bare
    mup4 = ROOT.TLorentzVector(mu.p4().Px(), mu.p4().Py(), mu.p4().Pz(), mu.p4().E() )
    nup4 = ROOT.TLorentzVector(nu.p4().Px(), nu.p4().Py(), nu.p4().Pz(), nu.p4().E() )

    # pre-FSR
    mup4_prefsr = ROOT.TLorentzVector(mu_prefsr.p4().Px(), mu_prefsr.p4().Py(), mu_prefsr.p4().Pz(), mu_prefsr.p4().E() )
    nup4_prefsr = nup4 if DY=='CC' else ROOT.TLorentzVector(nu_prefsr.p4().Px(), nu_prefsr.p4().Py(), nu_prefsr.p4().Pz(), nu_prefsr.p4().E() )

    # dressed
    mup4_recfsr = mup4
    for g in mu_fsr:
        mup4_recfsr += ROOT.TLorentzVector(g.p4().Px(), g.p4().Py(), g.p4().Pz(), g.p4().E() )
    nup4_recfsr = nup4
    if DY=='NC':
        for g in nu_fsr:
            nup4_recfsr += ROOT.TLorentzVector(g.p4().Px(), g.p4().Py(), g.p4().Pz(), g.p4().E() )

    # list of p4 for W
    Wp4 = {}
    Wp4['bare']  = mup4 + nup4                # the mu+nu p4 after FSR
    Wp4['preFSR']= mup4_prefsr + nup4_prefsr  # the mu+nu p4 pre-FSR
    Wp4['dress'] = mup4_recfsr + nup4_recfsr  # the mu+nu p4 w/ dressed mu

    for t in ['bare', 'dress', 'preFSR']:
        variables[t+'_mass'][0] = Wp4[t].M() 
        variables[t+'_qt'][0]   = Wp4[t].Pt() 
        variables[t+'_y'][0]    = Wp4[t].Rapidity() 
        variables[t+'_phi'][0]  = Wp4[t].Phi() 
        mup4LV = mup4
        nup4LV = nup4
        if t=='dress':
            mup4LV = mup4_recfsr
            nup4LV = nup4_recfsr
        elif t=='preFSR':
            mup4LV = mup4_prefsr
            nup4LV = nup4_prefsr
        variables[t+'_mu_pt'][0]  = mup4LV.Pt()
        variables[t+'_mu_eta'][0] = mup4LV.Eta()
        variables[t+'_mu_phi'][0] = mup4LV.Phi()
        variables[t+'_nu_pt'][0]  = nup4LV.Pt()
        variables[t+'_nu_eta'][0] = nup4LV.Eta()
        variables[t+'_nu_phi'][0] = nup4LV.Phi()
        
        # the CS variables
        ps = boost_to_CS_root(Lp4=mup4LV, Wp4=Wp4[t] ) 
        variables[t+'_ECS'][0]   = ps[0] 
        variables[t+'_cosCS'][0] = ps[1] 
        variables[t+'_phiCS'][0] = ps[2] 

    # fill the tree
    outtree.Fill()

stop = time.time()
# save and close
outfile.cd()
outtree.Write("tree", ROOT.TObject.kOverwrite)
if debug:
    add_vars(tree=outtree, debug=True)
outfile.Close()
print "Output file closed. Processed ", i , "events in "+"{:03.0f}".format(stop-start)+" sec.("+"{:03.0f}".format(i/(stop-start))+" Hz)" 
