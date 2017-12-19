import os.path

from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import pickle
from pprint import pprint

ROOT.TGaxis.SetMaxDigits(2)

jobs = [ 
    ['1e6_pt_y_A0_A4_prior03',  40],
    ['1e7_pt_y_A0_A4_prior02',  40],
    #['1e6_pt_y_A0_A4_mass_prior02',  40],
    #['1e6_pt_y_A0_A4_mass_prior03', 40],
    #['1e7_pt_y_A0_A4_mass_prior01', 40],
    #['1e7_pt_y_A0_A4_mass_prior02', 40],
    #['1e6_pt_y_A0_A4_mass_prior02',  40],
    #['1e6_pt_y_A0_A4_mass_prior02',  40],
    #['1e6_pt_y_A0_A4_mass_prior03',  40],
    #['5e6_pt_y_A0_A4_mass_prior03',  40],
    #['1e7_pt_y_A0_A4_mass_prior03',  40],
    #['1e6_pt_y_A0_A4_mass_prior02',  40],
    #['5e6_pt_y_A0_A4_mass_prior02',  40],
    #['1e7_pt_y_A0_A4_mass_prior02',  40],
    #['1e6_pt_y_A0_A4_mass_prior01',  40],
    #['5e6_pt_y_A0_A4_mass_prior01',  40],
    #['1e7_pt_y_A0_A4_mass_prior01',  40],
    #['1e5_pt_y', 50], 
    #['1e6_pt_y', 50], 
    #['1e7_pt_y', 50], 
    #['1e5_pt_y_A0_A4', 50], 
    #['1e6_pt_y_A0_A4', 50], 
    #['1e7_pt_y_A0_A4', 50], 
    #['1e5_pt_y_A0_A4_prior05', 50], 
    #['1e6_pt_y_A0_A4_prior05', 50], 
    #['1e7_pt_y_A0_A4_prior05', 50], 
    #['1e5_pt_y_A0_A4_mass', 50], 
    #['1e6_pt_y_A0_A4_mass', 50], 
    #['1e7_pt_y_A0_A4_mass', 50], 
    #['1e5_y_A0_A4', 50], 
    #['1e6_y_A0_A4', 50], 
    #['1e7_y_A0_A4', 50], 
    ]

dir_name = 'V2'

for ijob,job in enumerate(jobs):
    #print job
    histos = {} 
    first = 0
    for j in range(job[1]):
        if not os.path.isfile('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl'):
            continue
        f = open('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl', 'r')
        res = pickle.load(f)
        #pprint(a)
        for key,p in res.items():        

            if not (key=='mass' or 'pt' in key):
                continue

            if first==0: 
                #print "fill for ", key
                histos[key+'_pull'] = ROOT.TH1F(key+'_pull', key+'_pull', 21, -4, 4) 
                ranges = [0., 1.0]
                if key=='mass':
                    ranges = [0., 1e-03]
                elif '_A' in key:
                    ranges = [-1.0, 1.0]
                else:
                    ranges = [0.0, 0.50]
                histos[key+'_relerr'] = ROOT.TH1F(key+'_relerr', key+'_relerr', 100, ranges[0], ranges[1]) 

            for toy in range(res['ntoys']):        
                if res['status'][toy] != 0:
                    continue
                pull = (res[key]['fit'][toy]-res[key]['true'][toy])/res[key]['err'][toy] if res[key]['err'][toy]>0 else 0.0            
                relerr = res[key]['err'][toy]/res[key]['fit'][toy] if '_A' not in key and res[key]['fit'][toy]>0. else res[key]['err'][toy]
                histos[key+'_pull'].Fill(pull)
                histos[key+'_relerr'].Fill(relerr)
        f.close()
        first = 1

    for key, p in res.items():
        if not (histos.has_key(key+'_pull') or histos.has_key(key+'_relerr') ):
            continue
        if 'pt_y' in key and 'pt_y' not in job[0]:
            continue
        if 'A' in key and 'A' not in job[0]:
            continue
        if 'mass' in key and 'mass' not in job[0]:
            continue

        print "Plot", key
        c = ROOT.TCanvas("c_"+job[0]+"_"+key, "canvas for "+job[0]+"_"+key, 1200, 600) 
        c.Divide(2,1)
        p = None
        for var in ['pull', 'relerr']:
            if var=='pull':
                c.cd(1)
                p = histos[key+'_pull']
            elif var=='relerr':
                c.cd(2)
                p = histos[key+'_relerr']
            title = var.title()+' of '+key+', N_{ev}='+str(res['num_events'])+', fix: '
            for ic,k in enumerate( res['fix'] ):
                title += (k+',')
            p.SetTitle(title)
            p.GetXaxis().SetTitle(var.title())
            p.GetXaxis().SetTitleSize(25)
            p.GetXaxis().SetTitleFont(43)
            p.GetXaxis().SetTitleOffset(1.0)
            p.GetXaxis().SetLabelFont(43) 
            p.GetXaxis().SetLabelSize(20)
            p.GetYaxis().SetTitle('toys')
            p.GetYaxis().SetTitleSize(25)
            p.GetYaxis().SetTitleFont(43)
            p.GetYaxis().SetTitleOffset(1.0)
            p.GetYaxis().SetLabelFont(43) 
            p.GetYaxis().SetLabelSize(20)
            p.SetLineWidth(3)
            p.SetLineStyle(ROOT.kSolid)
            p.SetLineColor(ROOT.kRed)
            p.Draw("HIST")
            #raw_input()
        c.SaveAs('plots/'+dir_name+'/'+key+'_'+job[0]+'.png')
        c.IsA().Destructor( c )
