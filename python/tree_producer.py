from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )
from ROOT import TLorentzVector

import time
import math
import os
import pickle

from tree_utils import *

ROOT.gSystem.Load("libFWCoreFWLite.so")
ROOT.gSystem.Load("libDataFormatsFWLite.so")
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Handle, Events

class TreeProducer:

    def __init__(self, DY='CC', verbose=False, debug=True, filenames=[], postfix='test', save_tree=True, save_histo1=False, save_histo2=False, save_histo3=False, masses=[80.419]):
    
        print "****** TreeProducer *****"
        print 'Running for '+DY+' Drell-Yan'
        self.DY = DY
        self.verbose = verbose
        self.debug = debug

        # save the TTree in the output file
        self.save_tree = save_tree

        # open output file
        self.outfile = None
        if self.debug:
            self.outfile = ROOT.TFile(os.environ['CMSSW_BASE']+'/src/Wmass/test/'+'tree_'+postfix+'.root', "RECREATE")
        else:
            self.outfile = ROOT.TFile('tree.root', "RECREATE")

        # out tree
        self.outtree = None
        if self.save_tree:
            self.outtree = ROOT.TTree('tree', 'tree')

        # mW for reweighting
        self.masses = masses

        # add histos
        #self.weights_for_histos = [0]
        self.weights_for_histos = range(109)
        self.coefficients_for_histos = ['A0','A1','A2','A3','A4','A5','A6','A7']

        # save the histograms in the output file
        self.save_histo1 = save_histo1 
        if self.save_histo1:
            self.histos = add_histo2D(charges=['Wminus','Wplus'], var=['Wdress', 'Wbare', 'WpreFSR'], 
                                      coeff=self.coefficients_for_histos, 
                                      weights=self.weights_for_histos)

        # read result on fit parameters
        self.save_histo2 = save_histo2
        self.fit_result = {}
        if save_histo2:
            for q in ['Wplus', 'Wminus']:
                self.fit_result[q]  = {}
                for var in ['Wdress', 'Wbare', 'WpreFSR']:
                    self.fit_result[q][var] = pickle.load( open(os.environ['CMSSW_BASE']+'/src/Wmass/data/'+'fit_results_'+DY+'_'+q+'_'+var+'_all_A0-7.pkl') )
            self.histos = add_histo2D_CS( charges=['Wminus','Wplus'], var=['Wdress', 'Wbare', 'WpreFSR'], coeff_eval=['val'])

        # pt-eta
        self.save_histo3 = save_histo3
        if save_histo3:
            for q in ['Wplus', 'Wminus']:
                self.fit_result[q]  = {}
                for var in ['Wdress', 'Wbare', 'WpreFSR']:
                    self.fit_result[q][var] = pickle.load( open(os.environ['CMSSW_BASE']+'/src/Wmass/data/'+'fit_results_'+DY+'_'+q+'_'+var+'_all_A0-7.pkl') )
            self.histos = add_histo2D_lepton( charges=['Wminus','Wplus'], var=['WpreFSR'], coeff_eval=['val'], masses=masses, coeff=['A0','A1','A2','A3','A4','A5','A6','A7'] ) 

        # add branches to tree (needed even if self.save_tree=False)
        self.variables = add_vars(self.outtree)

        self.filenames = filenames
        if self.debug:
            if 'NC' in self.DY:
                self.filenames.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root')
            elif self.DY=='CC_FxFx':
                self.filenames.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/0AF0207B-EFBE-E611-B4BE-0CC47A7FC858.root')
            elif self.DY=='CC_MG5':
                self.filenames.append('root://xrootd-cms.infn.it//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/100000/00933E2A-A0D5-E611-B2CD-00266CF89130.root')

        print "Opening file..."
        self.events = Events(self.filenames)
        print "File opened.... Tot. num of events:", self.events.size()


    def run(self):

        # handles
        genH, genN = Handle("std::vector<reco::GenParticle>"), "prunedGenParticles"
        infoH, infoN = Handle("GenEventInfoProduct"), "generator"
        lheH, lheN = Handle("LHEEventProduct"), "externalLHEProducer"

        start = time.time()

        # loop over the events
        ###############################################
        for i,event in enumerate(self.events):    

            if self.debug:
                if i%100==0:
                    print "Processing event", i, '/', self.events.size()
                if i>1000:
                    break
            else:
                if i%1000==0:
                    print "Processing event", i, '/', self.events.size()


            # fill with default values
            fill_default(self.variables)
            
            # load LHE product
            event.getByLabel(lheN,lheH)
            lhe = lheH.product()
            hepeup = lhe.hepeup()

            self.variables['scale'][0]    = hepeup.SCALUP
            self.variables['alphaQCD'][0] = hepeup.AQCDUP
            self.variables['alphaQED'][0] = hepeup.AQEDUP

            isW = False
            isWToMuNu = False

            Wp4_lhe = [0.,0.,0.,0.,0.]
            for p in range(hepeup.NUP):
                if self.verbose:
                    print 'HEPEUP...', p, 'pdg:', hepeup.IDUP[p], '(', hepeup.PUP[p][0], hepeup.PUP[p][1],hepeup.PUP[p][2],hepeup.PUP[p][3],')'
                # first parton
                if p==0:
                    self.variables['id1'][0] = hepeup.IDUP[p]
                    self.variables['x1'][0] = hepeup.PUP[p][3]/6500.
                elif p==1:
                    self.variables['id2'][0] = hepeup.IDUP[p] 
                    self.variables['x2'][0] = hepeup.PUP[p][3]/6500.
                if abs(hepeup.IDUP[p])==13:
                    isWToMuNu = True
                if abs(hepeup.IDUP[p])== (24 if 'CC' in self.DY else 23):
                    isW = True
                    for x in range(5):                
                        Wp4_lhe[x] = hepeup.PUP[p][x]
                    Wp4 = ROOT.TLorentzVector(Wp4_lhe[0], Wp4_lhe[1], Wp4_lhe[2], Wp4_lhe[3])     
                    self.variables['lhe_y'][0]   =  Wp4.Rapidity()
                    self.variables['lhe_qt'][0]  =  Wp4.Pt()
                    self.variables['lhe_mass'][0]=  Wp4.M()
                    self.variables['lhe_phi'][0] =  Wp4.Phi()
   
            # consider only W->munu events
            if not (isW and isWToMuNu):
                continue

            # LHE weights
            norm = abs(lhe.originalXWGTUP())
            wid=0
            for w in lhe.weights():
                if wid>=109:
                    continue
                #if self.verbose:
                #    print w.id, w.wgt
                self.variables['weights'][wid] = w.wgt/norm
                wid += 1

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
            if 'CC' in self.DY and len(muons)==0:
                if self.verbose:
                    print "No muon in W>munu for event", i, ". Try to understand why:"
                    print " > W rapidity:    ", self.variables['lhe_y'][0]
                self.variables['muLost'][0] += 1.0
                if len(neutrinos)>0:
                    self.variables[t+'_nu_pt'][0]  = neutrinos[0].p4().Pt()
                    self.variables[t+'_nu_eta'][0] = neutrinos[0].p4().Eta()
                    self.variables[t+'_nu_phi'][0] = neutrinos[0].p4().Phi()

            # consider events with 0 neutrinos to study acceptance 
            if 'CC' in self.DY and len(neutrinos)==0:
                if self.verbose:
                    print "No muon neutrinos for event", i, ". Try to understand why:"
                    print " > W rapidity:    ", self.variables['lhe_y'][0]
                    print " > Muon rapidity: ", muons[0].p4().Rapidity()
                self.variables['nuLost'][0] += 1.0
                if len(muons)>0:
                    self.variables[t+'_mu_pt'][0]  = muons[0].p4().Pt()
                    self.variables[t+'_mu_eta'][0] = muons[0].p4().Eta()
                    self.variables[t+'_mu_phi'][0] = muons[0].p4().Phi()

            # NC: consider events with 1 muons to study acceptance 
            if 'NC' in self.DY and len(muons)<2:
                if self.verbose:
                    print "No muons in Z>mumu for event", i, ". Try to understand why:"
                    print " > Z rapidity:    ", self.variables['lhe_y'][0]
                    self.variables['muLost'][0] += 1.0
                if len(muons)==1:
                    self.variables[t+'_nu_pt'][0]  = muons[0].p4().Pt()
                    self.variables[t+'_nu_eta'][0] = muons[0].p4().Eta()
                    self.variables[t+'_nu_phi'][0] = muons[0].p4().Phi()

            # if no muons or no neutrinos, save the event and continue
            if 'CC' in self.DY:
                if (len(muons)==0 or len(neutrinos)==0):
                    if self.save_tree:
                        self.outtree.Fill()
                    continue
            elif 'NC' in self.DY:
                if len(muons)<2:
                    if self.save_tree:
                        self.outtree.Fill()
                    continue
        
            # the muon is the first ranked by pt if CC else the mu+
            mu = muons[0]
            if 'NC' in self.DY:
                mu = (muons[0] if muons[0].pdgId()==-13 else muons[1])
            self.variables['isFromW'][0] += int(isFromW(mu) if 'CC' in self.DY else isFromZ(mu))
            self.variables['mu_charge'][0] = mu.pdgId()
            if self.verbose:
                printp('muon', mu, '')

            # the neutrino is the first ranked by pt if CC else the mu-
            nu = None
            if 'CC' in self.DY:
                nu = neutrinos[0]
            elif 'NC' in self.DY:
                nu = (muons[0] if muons[0].pdgId()==+13 else muons[1])
            self.variables['isFromW'][0] += int(isFromW(nu) if 'CC' in self.DY else isFromZ(nu))
            self.variables['nu_charge'][0] = nu.pdgId()
            if self.verbose:
                printp('neut', nu, '')

            # pre-FSR
            mother = mu
            mu_prefsr = mu
            while(mother.numberOfMothers()>0):
                if self.verbose:
                    printp('MOTH', mother, '')
                if abs(mother.pdgId())==13 and mother.statusFlags().isLastCopyBeforeFSR():
                    mu_prefsr = mother
                mother = mother.mother(0)

            mother = nu
            nu_prefsr = nu
            while(mother.numberOfMothers()>0):
                if self.verbose:
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
                dR = (dR_mu if 'CC' in self.DY else min(dR_mu,dR_nu))
                if dR<0.1:
                    if 'CC' in self.DY:
                        mu_fsr.append(g.p4())
                        if self.verbose:
                            printp('>gam', g, 'dR:{:03.2f}'.format(dR)+' to muon')
                    elif 'NC' in self.DY:
                        if dR_mu<dR_nu:
                            mu_fsr.append(g.p4())
                            if self.verbose:
                                printp('>gam', g, 'dR:{:03.2f}'.format(dR)+' to muon-')
                        else:
                            nu_fsr.append(g.p4())
                            if self.verbose:
                                printp('>gam', g, 'dR:{:03.2f}'.format(dR)+' to mu+')

            # bare
            mup4 = ROOT.TLorentzVector(mu.p4().Px(), mu.p4().Py(), mu.p4().Pz(), mu.p4().E() )
            nup4 = ROOT.TLorentzVector(nu.p4().Px(), nu.p4().Py(), nu.p4().Pz(), nu.p4().E() )
            
            # pre-FSR
            mup4_prefsr = ROOT.TLorentzVector(mu_prefsr.p4().Px(), mu_prefsr.p4().Py(), mu_prefsr.p4().Pz(), mu_prefsr.p4().E() )
            nup4_prefsr = copy.deepcopy(nup4) if 'CC' in self.DY else ROOT.TLorentzVector(nu_prefsr.p4().Px(), nu_prefsr.p4().Py(), nu_prefsr.p4().Pz(), nu_prefsr.p4().E() )

            # dressed
            mup4_recfsr = copy.deepcopy(mup4)
            for g in mu_fsr:
                mup4_recfsr += ROOT.TLorentzVector(g.Px(), g.Py(), g.Pz(), g.E() )
            nup4_recfsr = copy.deepcopy(nup4)
            if 'NC' in self.DY:
                for g in nu_fsr:
                    nup4_recfsr += ROOT.TLorentzVector(g.Px(), g.Py(), g.Pz(), g.E() )

            # list of p4 for W
            Wp4 = {}
            Wp4['Wbare']  = mup4 + nup4                # the mu+nu p4 after FSR
            Wp4['WpreFSR']= mup4_prefsr + nup4_prefsr  # the mu+nu p4 pre-FSR
            Wp4['Wdress'] = mup4_recfsr + nup4_recfsr  # the mu+nu p4 w/ dressed mu

            for t in ['Wbare', 'Wdress', 'WpreFSR']:
                self.variables[t+'_mass'][0] = Wp4[t].M() 
                self.variables[t+'_qt'][0]   = Wp4[t].Pt() 
                self.variables[t+'_y'][0]    = Wp4[t].Rapidity() 
                self.variables[t+'_phi'][0]  = Wp4[t].Phi() 
                mup4LV = mup4
                nup4LV = nup4
                if t=='Wdress':
                    mup4LV = mup4_recfsr
                    nup4LV = nup4_recfsr
                elif t=='WpreFSR':
                    mup4LV = mup4_prefsr
                    nup4LV = nup4_prefsr
                self.variables[t+'_mu_pt'][0]  = mup4LV.Pt()
                self.variables[t+'_mu_eta'][0] = mup4LV.Eta()
                self.variables[t+'_mu_phi'][0] = mup4LV.Phi()
                self.variables[t+'_nu_pt'][0]  = nup4LV.Pt()
                self.variables[t+'_nu_eta'][0] = nup4LV.Eta()
                self.variables[t+'_nu_phi'][0] = nup4LV.Phi()
        
                # the CS variables
                ps = boost_to_CS_root(Lp4=mup4LV, Wp4=Wp4[t] ) 
                self.variables[t+'_ECS'][0]   = ps[0] 
                self.variables[t+'_cosCS'][0] = ps[1] 
                self.variables[t+'_phiCS'][0] = ps[2] 

                if self.save_histo1:
                    if t not in ['Wdress', 'Wbare', 'WpreFSR']:
                        continue
                    for w in self.weights_for_histos:
                        q1 = 'Wplus' if self.variables['mu_charge'][0]==-13 else 'Wminus'
                        fill_coefficients(histos=self.histos, 
                                          q=q1, 
                                          coefficients_for_histos=self.coefficients_for_histos,
                                          var=t,
                                          weight_name=w, 
                                          ps_W=(Wp4[t].Rapidity(), Wp4[t].Pt()), 
                                          ps_CS=(ps[1],ps[2]), 
                                          weight=self.variables['weights'][w] )
                        if 'NC' in self.DY:
                            q2 = 'Wplus' if self.variables['nu_charge'][0]==-13 else 'Wminus'
                            fill_coefficients(histos=self.histos, 
                                              q=q2,
                                              coefficients_for_histos=self.coefficients_for_histos,
                                              var=t,
                                              weight_name=w, 
                                              ps_W=(Wp4[t].Rapidity(), Wp4[t].Pt()), 
                                              ps_CS=(-ps[1], ps[2] + math.pi - (2*math.pi if (ps[2] + math.pi) > 2*math.pi else 0.0 ) ), 
                                              weight=self.variables['weights'][w] )

                if self.save_histo2:
                    if t not in ['Wdress', 'Wbare', 'WpreFSR']:
                        continue
                    q1 = 'Wplus' if self.variables['mu_charge'][0]==-13 else 'Wminus'
                    fill_weighted_CS(res=getattr(self, "fit_result")[q1][t],
                                     histos=self.histos, 
                                     q=q1,
                                     var=t,
                                     coeff_eval=['val'],
                                     ps_W=(Wp4[t].Rapidity(), Wp4[t].Pt()), 
                                     ps_CS=(ps[1],ps[2]),
                                     weight=self.variables['weights'][0],
                                     coeff=['A0', 'A1', 'A2', 'A3', 'A4']
                                     )

                if self.save_histo3:
                    if t not in ['WpreFSR']:
                        continue
                    q1 = 'Wplus' if self.variables['mu_charge'][0]==-13 else 'Wminus'
                    fill_lepton_lab(res=getattr(self, "fit_result")[q1][t],
                                    histos=self.histos, 
                                    q=q1,
                                    var=t,
                                    coeff_eval=['val'],
                                    masses=self.masses,
                                    ps_W=(Wp4[t].Rapidity(), Wp4[t].Pt(), Wp4[t].M()), 
                                    ps_CS=(ps[1],ps[2]),
                                    ps_lep=(self.variables['Wbare_mu_eta'][0], self.variables['Wbare_mu_pt'][0]),
                                    weight=self.variables['weights'][0],
                                    coeff=['A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7']
                                    )


            # fill the tree
            if self.save_tree:
                self.outtree.Fill()
        ###############################################
        
        stop = time.time()

        # save and close
        self.outfile.cd()
        if self.save_tree:
            self.outtree.Write("tree", ROOT.TObject.kOverwrite)

        if self.save_histo1 or self.save_histo2:
            for kq,q in self.histos.items():
                print 'Charge: '+kq
                self.outfile.mkdir(kq)
                for kv,v in q.items():
                    print '\tVar: '+kv
                    self.outfile.mkdir(kq+'/'+kv)
                    for kc,c in v.items():
                        print '\t\tVar: '+kc
                        self.outfile.mkdir(kq+'/'+kv+'/'+kc)
                        for kw,w in c.items():
                            print '\t\t\tHisto: '+kw+'.....', w[0].GetEntries(), 'entries'
                            self.outfile.cd(kq+'/'+kv+'/'+kc)
                            w[0].Write('', ROOT.TObject.kOverwrite)
                            w[1].Write('', ROOT.TObject.kOverwrite)
                self.outfile.cd()

        if self.save_histo3:
            for kq,q in self.histos.items():
                print 'Charge: '+kq
                self.outfile.mkdir(kq)
                for kv,v in q.items():
                    print '\tVar: '+kv
                    self.outfile.mkdir(kq+'/'+kv)
                    for kc,c in v.items():
                        print '\t\tEval: '+kc
                        self.outfile.mkdir(kq+'/'+kv+'/'+kc)
                        for km,m in c.items():
                            print '\t\t\tMass: '+km
                            self.outfile.mkdir(kq+'/'+kv+'/'+kc+'/'+km)
                            for kb,b in m.items():
                                print '\t\t\t\tBin: '+kb
                                self.outfile.mkdir(kq+'/'+kv+'/'+kc+'/'+km+'/'+kb)
                                for ka,a in b.items():
                                    print '\t\t\t\t\tCoeff: '+ka+'.....', a.GetEntries(), 'entries'                                    
                                    self.outfile.cd(kq+'/'+kv+'/'+kc+'/'+km+'/'+kb)
                                    a.Write('', ROOT.TObject.kOverwrite)
                self.outfile.cd()


        if self.debug and self.save_tree:
            add_vars(tree=self.outtree, debug=True)

        self.outfile.Close()
        print "Output file closed. Processed ", i , "events in "+"{:03.0f}".format(stop-start)+" sec.("+"{:03.0f}".format(i/(stop-start))+" Hz)" 
