from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )
from ROOT import TLorentzVector

import math
import os

from tree_utils import *

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Handle, Events

verbose = False

#outfile = ROOT.TFile(os.environ['CMSSW_BASE']+'/src/Wmass/test/'+'tree_'+argv[1]+'.root', "RECREATE")
outfile = ROOT.TFile(os.environ['CMSSW_BASE']+'/src/Wmass/test/'+'tree.root', "RECREATE")
outtree = ROOT.TTree('tree', 'tree')
variables = add_vars(outtree)

print "Opening file..."
filename = [
    'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0AF0207B-EFBE-E611-B4BE-0CC47A7FC858.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/002F2CE1-38BB-E611-AF9F-0242AC130005.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/009CE684-45BB-E611-A261-001E67E6F8FA.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/044FF9CC-42BB-E611-ACB0-0CC47AD98BC2.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/06103109-48BB-E611-86BE-001E673968A6.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0843C79F-FCBD-E611-B38C-001E67A3F8A8.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0881BCD8-8FBE-E611-8796-002590FD5A72.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08E524F3-0ABC-E611-984F-141877639F59.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/08F5FD50-23BC-E611-A4C2-00259073E3DA.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0A85AA82-45BB-E611-8ACD-001E674FB063.root',
    #'root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0CA050B2-57BB-E611-8A7A-001E674FBA1D.root'
    ]
events = Events(filename)
print "File opened.... Tot. num of events:", events.size()


genH, genN = Handle("std::vector<pat::PackedGenParticle>"), "packedGenParticles"
infoH, infoN = Handle("GenEventInfoProduct"), "generator"
lheH, lheN = Handle("LHEEventProduct"), "externalLHEProducer"

for i,event in enumerate(events):    

    if i%1000==0:
        print "Processing event", i, '/', events.size()
    
    event.getByLabel(genN,genH)
    event.getByLabel(lheN,lheH)

    #event.getByLabel(infoN,infoH)
    #generator = infoH.product()
    #weights = list(generator.weights())
    #variables['weight'][0] = weights[0]/abs(weights[0])

    lhe = lheH.product()
    norm = abs(lhe.weights()[0].wgt)
    variables['weight'][0]   = lhe.weights()[0].wgt/norm
    variables['weightNU'][0] = lhe.weights()[1].wgt/norm
    variables['weightND'][0] = lhe.weights()[2].wgt/norm
    variables['weightUN'][0] = lhe.weights()[3].wgt/norm
    variables['weightUU'][0] = lhe.weights()[4].wgt/norm
    variables['weightDN'][0] = lhe.weights()[6].wgt/norm
    variables['weightDD'][0] = lhe.weights()[8].wgt/norm
   
    if verbose:
        print "Weights:"
        for w in ['', 'NU', 'ND', 'UN', 'UU', 'DN', 'DD']:
            print w, variables['weight'+w][0]


    hepeup = lhe.hepeup()
    Wp4_lhe = [0.,0.,0.,0.,0.]
    isMuThere = False
    for p in range(hepeup.NUP):
        if abs(hepeup.IDUP[p])==13:
            isMuThere = True
        if abs(hepeup.IDUP[p])==24:
            for x in range(5):                
                Wp4_lhe[x] = hepeup.PUP[p][x]

    genParticles = list(genH.product())
    muons = [p for p in genParticles if isMuon(p)] 
    neutrinos = [p for p in genParticles if isNeutrino(p)] 
    gammas = [p for p in genParticles if isPhoton(p)] 
    
    variables['muLost'][0] = 0.0
    variables['nuLost'][0] = 0.0

    if len(muons)==0:
        if isMuThere:
            if verbose:
                print "No muon in W>munu for event", i, ". Try to understand why:"
                print " > W rapidity:    ", 0.5*math.log((Wp4_lhe[3]+Wp4_lhe[2])/(Wp4_lhe[3]-Wp4_lhe[2]))
            fill_default(variables)
            variables['muLost'][0] = 1.0
            Wp4 = ROOT.TLorentzVector(Wp4_lhe[0], Wp4_lhe[1], Wp4_lhe[2], Wp4_lhe[3])    
            variables['lhe_y'][0] = Wp4.Rapidity()
            variables['lhe_qt'][0] =  Wp4.Pt()
            variables['lhe_mass'][0] =  Wp4.M()
            variables['lhe_phi'][0] =  Wp4.Phi()
            outtree.Fill()
        continue

    if len(neutrinos)==0:
        if verbose:
            print "No muon neutrinos for event", i, ". Try to understand why:"
            print " > W rapidity:    ", 0.5*math.log((Wp4_lhe[3]+Wp4_lhe[2])/(Wp4_lhe[3]-Wp4_lhe[2]))
            print " > Muon rapidity: ", muons[0].p4().Rapidity()
        fill_default(variables)
        variables['nuLost'][0] = 1.0
        Wp4 = ROOT.TLorentzVector(Wp4_lhe[0], Wp4_lhe[1], Wp4_lhe[2], Wp4_lhe[3])    
        variables['lhe_y'][0] = Wp4.Rapidity()
        variables['lhe_qt'][0] =  Wp4.Pt()
        variables['lhe_mass'][0] =  Wp4.M()
        variables['lhe_phi'][0] =  Wp4.Phi()
        outtree.Fill()
        continue

    variables['isW'][0] = 0.0

    muons.sort(reverse=True)
    mu = muons[0]
    variables['isW'][0] += int(isFromW(mu))

    variables['charge'][0] = mu.pdgId()
    if verbose:
        print '********************************'
        printp('muon', mu, '')

    neutrinos.sort(reverse=True)
    nu = neutrinos[0]
    variables['isW'][0] += int(isFromW(nu))
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

    Wp4 = {}
    nup4 = ROOT.TLorentzVector(nu.p4().Px(), nu.p4().Py(), nu.p4().Pz(), nu.p4().E() )
    mup4 = ROOT.TLorentzVector(mu.p4().Px(), mu.p4().Py(), mu.p4().Pz(), mu.p4().E() )
    mup4_prefsr = ROOT.TLorentzVector(mu_prefsr.p4().Px(), mu_prefsr.p4().Py(), mu_prefsr.p4().Pz(), mu_prefsr.p4().E() )
    mup4_recfsr = mup4
    for g in fsr:
        mup4_recfsr += ROOT.TLorentzVector(g.p4().Px(), g.p4().Py(), g.p4().Pz(), g.p4().E() )

    # the W p4 from LHE
    Wp4['lhe'] = ROOT.TLorentzVector(Wp4_lhe[0], Wp4_lhe[1], Wp4_lhe[2], Wp4_lhe[3])   
    # the mu+nu p4 after FSR
    Wp4['postFSR'] = mup4+nup4
    # the mu+nu p4 pre-FSR
    Wp4['preFSR'] = mup4_prefsr+nup4
    # the mu+nu p4 w/ dressed mu
    Wp4['dressFSR'] = mup4_recfsr+nup4

    for t in ['lhe', 'preFSR', 'dressFSR', 'postFSR']:
        variables[t+'_mass'][0] = Wp4[t].M() 
        variables[t+'_qt'][0]   = Wp4[t].Pt() 
        variables[t+'_y'][0]    = Wp4[t].Rapidity() 
        variables[t+'_phi'][0]  = Wp4[t].Phi() 
        mup4LV = mup4
        if t=='dressFSR':
            mup4LV = mup4_recfsr
        elif t=='preFSR':
            mup4LV = mup4_prefsr
        variables[t+'_mu_pt'][0]  = mup4LV.Pt()
        variables[t+'_mu_eta'][0] = mup4LV.Eta()
        variables[t+'_mu_phi'][0] = mup4LV.Phi()
        ps = boost_to_CS_root(Lp4=mup4LV, Wp4=Wp4[t] ) 
        variables[t+'_ECS'][0]   = ps[0] 
        variables[t+'_cosCS'][0] = ps[1] 
        variables[t+'_phiCS'][0] = ps[2] 

    outtree.Fill()

outfile.cd()
outtree.Write("tree", ROOT.TObject.kOverwrite)
outfile.Close()
