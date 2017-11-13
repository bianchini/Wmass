import ROOT
from ROOT import TFile, TTree, TH1F
from array import array

from minimize_deltas import Minimze_deltas
ROOT.gROOT.SetBatch(True)

def create_tree(n_bins=10, x_low=-2.0, x_high=+2.0, 
                n_events=500, n_toys=100, 
                pdf_x='Gaus', pdf_x_par=[0.0,1.0,-99.], 
                pdf_smear='Gaus', pdf_smear_par=[0.05,0.01,-99.], 
                postfix='lowstat'):
 
    f = TFile( 'deltas/gen_toy_'+str(n_bins)+'bin_'+str(n_events)+'ev_'+postfix+'.root', 'recreate' )
    t = TTree( 'tree', 'tree for toys' )

    ran = ROOT.TRandom3()

    h_n = TH1F( 'h_n', 'nominal', n_bins, x_low, x_high )

    n_deltas = array( 'i', [ 1 ] )
    n_bin = array( 'i', [ 1 ] )
    n_toy = array( 'i', [ 1 ] )
    deltas = array( 'i', (4*n_bins+4)*[ 0 ] )
    bin_n = array( 'f', (n_bins+2)*[ 0. ] )
    bin_u = array( 'f', (n_bins+2)*[ 0. ] )
    bin_d = array( 'f', (n_bins+2)*[ 0. ] )

    t.Branch( 'n_deltas', n_deltas, 'n_deltas/I' )
    t.Branch( 'deltas', deltas, 'deltas[n_deltas]/I' )

    t.Branch( 'n_bin', n_bin, 'n_bin/I' )
    t.Branch( 'n_toy', n_toy, 'n_toy/I' )
    t.Branch( 'bin_n', bin_n, 'bin_n[n_bin]/F' )
    t.Branch( 'bin_u', bin_u, 'bin_u[n_bin]/F' )
    t.Branch( 'bin_d', bin_d, 'bin_d[n_bin]/F' )

    for toy in range(n_toys):

        n_toy[0] = toy
        n_deltas[0] = 4*n_bins+4
        n_bin[0] = n_bins

        # clean
        for i in range( n_bin[0] ):
            bin_n[i] = 0.
            bin_u[i] = 0.
            bin_d[i] = 0.        
        for i in range( n_deltas[0] ):
            deltas[i]=0

        n_events_rnd = ran.Poisson(n_events)
        for i in range(n_events_rnd):
            x_n = 0.0
            if pdf_x=="Gaus":
                x_n = getattr(ran, pdf_x)(pdf_x_par[0], pdf_x_par[1])
            elif pdf_x=="Exp":
                x_n = getattr(ran, pdf_x)(pdf_x_par[0])
            smear = getattr(ran, pdf_smear)(pdf_smear_par[0], pdf_smear_par[1])    
            x_u = x_n + smear
            x_d = x_n - smear    
            ibin_n = h_n.FindBin(x_n)
            ibin_u = h_n.FindBin(x_u)
            ibin_d = h_n.FindBin(x_d)
            bin_n[ibin_n] += 1.0
            bin_u[ibin_u] += 1.0
            bin_d[ibin_d] += 1.0

            deltas[0] += (ibin_n==0 and ibin_u==1)
            deltas[1] += (ibin_n==0 and ibin_d==1)
            for b in range(1, n_bins+1):
                deltas[2+0*n_bins + b - 1] += (ibin_n==b and ibin_u==b+1 )
                deltas[2+1*n_bins + b - 1] += (ibin_n==b and ibin_d==b+1 )
                deltas[2+2*n_bins + b - 1] += (ibin_n==b and ibin_u==b-1 )
                deltas[2+3*n_bins + b - 1] += (ibin_n==b and ibin_d==b-1 )
            deltas[ (4*n_bins+4) - 2] += ibin_n==n_bins+1 and ibin_u==n_bins
            deltas[ (4*n_bins+4) - 1] += ibin_n==n_bins+1 and ibin_d==n_bins

        t.Fill()

    f.Write()
    f.Close()
    ran.IsA().Destructor(ran)

##############################################################


def test_minimizer(fin_name='gen_toy_10bin_highstat', fout_name='result'): 

    fout = TFile('deltas/'+fout_name+'.root', 'RECREATE')
    tout = TTree('tree', 'tree')
    bin =  array( 'i', [ 0 ] )
    n =  array( 'f', [ 0 ] )
    delta_n_u = array( 'f', [ 0. ] )
    delta_n_d = array( 'f', [ 0. ] )
    delta_m_u = array( 'f', [ 0. ] )
    delta_m_d = array( 'f', [ 0. ] )
    var_delta_m = array( 'f', [ 0. ] )

    tout.Branch( 'bin', bin, 'bin/I' )
    tout.Branch( 'n', n, 'n/F' )
    tout.Branch( 'delta_n_u', delta_n_u, 'delta_n_u/F' )
    tout.Branch( 'delta_n_d', delta_n_d, 'delta_n_d/F' )
    tout.Branch( 'delta_m_u', delta_m_u, 'delta_m_u/F' )
    tout.Branch( 'delta_m_d', delta_m_d, 'delta_m_d/F' )
    tout.Branch( 'var_delta_m', var_delta_m, 'var_delta_m/F' )

    f = TFile.Open('deltas/'+fin_name+".root")
    tree = f.Get("tree")
    tree.SetBranchStatus("*", True)

    for iev in range(tree.GetEntries()):
        tree.GetEntry(iev)
        ev = tree 
        n_toy = ev.n_toy
        n_bins = ev.n_bin    
        n_deltas = ev.n_deltas
        deltas = ev.deltas
        bin_n = ev.bin_n
        bin_u = ev.bin_u
        bin_d = ev.bin_d

        print deltas
        fix_params = []
        #fix_params.extend([1, 2+3*n_bins])
        #fix_params.extend(range(n_bins+2, 2*n_bins+2))
        minimize_deltas = Minimze_deltas(deltas=deltas, 
                                         nbins=n_bins, 
                                         step=0.01, 
                                         fix_params=fix_params)
        res = minimize_deltas.run()
        for b in range(n_bins):
            bin[0] = b+1
            n[0] = bin_n[b]
            delta_n_u[0] = res[0][b]
            delta_n_d[0] = res[1][b]
            delta_m_u[0] = res[2][b]
            delta_m_d[0] = res[3][b]
            var_delta_m[0] = res[4][b]
            #print "Bin %s: [%s, %s] ==> %s +/- %s" % ( b+1, delta_n_u[0], delta_n_d[0], abs(delta_m_u[0]), math.sqrt(var_delta_m[0]))
            tout.Fill()

    f.Close()
    fout.cd()
    fout.Write()
    fout.Close()

##############################################################

def draw(title='', x_title='', y_title='', hists=[], legend=''):
    c1 = ROOT.TCanvas("c1_"+title,"c1",600,600)
    leg = ROOT.TLegend( 0.12  ,0.65, 0.45 ,0.88, "","brNDC")
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.06)
    leg.SetFillColor(10) 
    c1.cd()
    for hi in range(len(hists)):
        if hi==0:
            hists[hi].SetMinimum(  hists[hi].GetMinimum()*2.0 if 'RMS' not in y_title else 0.0 )
            hists[hi].SetMaximum(  hists[hi].GetMaximum()*2.0 )
            #hists[hi].SetMaximum(  hists[hi].GetMaximum()*1.4 )
            hists[hi].GetXaxis().SetTitle(x_title)
            hists[hi].GetYaxis().SetTitle(y_title)
            hists[hi].GetYaxis().SetTitleSize(25)
            hists[hi].GetYaxis().SetTitleOffset(0.8)
            hists[hi].GetYaxis().SetTitleFont(43)
            hists[hi].GetXaxis().SetTitleSize(25)
            hists[hi].GetXaxis().SetTitleFont(43)
            leg.AddEntry( hists[hi], 'nominal', 'L')
            #hists[hi].GetYaxis().SetLabelSize(30)
            #hists[hi].GetYaxis().SetLableFont(43)
            #hists[hi].GetXaxis().SetLableSize(20)
            #hists[hi].GetXaxis().SetLableFont(43)
            hists[hi].SetTitle("")
            hists[hi].SetStats(0)
            hists[hi].Draw('HISTE')
        elif hi==1:
            leg.AddEntry( hists[hi],  '0.5#times(u-d)' if '_u' in title else  '0.5#times(d-u)', 'L')
            hists[hi].Draw('HISTESAME')
        elif hi==2:
            #leg.AddEntry( hists[hi],  '0.5#times(|u|+|d|)#times sign('+('u-d' if '_u' in title else 'd-u')+')', 'L')
            #hists[hi].Draw('HISTESAME')
        #else:
            leg.AddEntry( hists[hi],  'fit', 'L')
            hists[hi].Draw('HISTESAME')
            
    leg.SetHeader(legend)
    leg.Draw()
    c1.SaveAs('deltas/'+title+'.png')
    c1.IsA().Destructor( c1 )
    leg.IsA().Destructor( leg )
    return


def make_plots(fin_name='10bin_500ev_gaus_5pc', n_bins=10, legend=''):
    
    fin = TFile('deltas/result_'+fin_name+'.root')
    tree = fin.Get("tree")
    fin.cd()

    histos = {}
    for key in ['n_u', 'n_d', 'm_u', 'm_d', 'av1_u', 'av1_d',  
                'av2_u', 'av2_d'
                ]:
        histos['h2_'+key] = ROOT.TH2F('h2_'+key, 'h2_'+key, n_bins, 1, n_bins+1, 10000, -100, 100)
        if 'av' not in key:
            tree.Draw('delta_'+key+':bin>>'+'h2_'+key)
        else:
            if key=='av1_u':
                tree.Draw('0.5*(delta_n_u-delta_n_d):bin>>'+'h2_'+key, "")
            elif key=='av1_d':
                tree.Draw('0.5*(delta_n_d-delta_n_u):bin>>'+'h2_'+key, "")
            elif key=='av2_u':
                tree.Draw('((abs(delta_n_u)+abs(delta_n_d))/2 * (-2*int(TMath::SignBit(delta_n_u-delta_n_d))+1)):bin>>'+'h2_'+key, "")
            elif key=='av2_d':
                tree.Draw('((abs(delta_n_u)+abs(delta_n_d))/2 * (-2*int(TMath::SignBit(delta_n_d-delta_n_u))+1)):bin>>'+'h2_'+key, "")
                
        #histos['h2_'+key].Draw('box')
        #raw_input()
        histos['h1_mean_'+key] = ROOT.TH1F('h1_mean_'+key, 'h1_'+key, n_bins, 1, n_bins+1)
        histos['h1_RMS_'+key] = histos['h1_mean_'+key].Clone('h1_RMS_'+key)
        for b in range(1,n_bins+1):
            hslice = histos['h2_'+key].ProjectionY("_y", b, b) 
            histos['h1_mean_'+key].SetBinContent(b, hslice.GetMean())
            histos['h1_mean_'+key].SetBinError(b, hslice.GetMeanError())
            histos['h1_RMS_'+key].SetBinContent(b, hslice.GetRMS())
            histos['h1_RMS_'+key].SetBinError(b, hslice.GetRMSError())

    for key in ['u', 'd']:
        for obs in ['mean', 'RMS']:
            histos['h1_'+obs+'_n_'+key].SetLineWidth(3)
            histos['h1_'+obs+'_n_'+key].SetLineColor(ROOT.kBlue)
            histos['h1_'+obs+'_m_'+key].SetLineWidth(3)
            histos['h1_'+obs+'_m_'+key].SetLineColor(ROOT.kRed)
            histos['h1_'+obs+'_av1_'+key].SetLineWidth(3)
            histos['h1_'+obs+'_av1_'+key].SetLineColor(ROOT.kGreen)
            histos['h1_'+obs+'_av2_'+key].SetLineWidth(3)
            histos['h1_'+obs+'_av2_'+key].SetLineColor(ROOT.kYellow)        
            histos['h1_'+obs+'_n_'+key].Draw()
            histos['h1_'+obs+'_m_'+key].Draw('same')
            histos['h1_'+obs+'_av1_'+key].Draw('same')
            histos['h1_'+obs+'_av2_'+key].Draw('same')
            draw(title='plot_'+obs+'_'+fin_name+'_'+key,  x_title='bin number', y_title=obs, 
                 hists=[ 
                    histos['h1_'+obs+'_n_'+key],  
                    histos['h1_'+obs+'_av1_'+key], 
                    #histos['h1_'+obs+'_av2_'+key], 
                    histos['h1_'+obs+'_m_'+key]
                    ],
                 legend=legend)
        #raw_input()
    return

##############################################################


cfg_base = {
    'gaus' : {
        'n_bins'       :  10,
        'x_low'        : -2.0,
        'x_high'       : +2.0,
        'n_events'     : 500,
        'n_toys'       : 1000,
        'pdf_x'        : "Gaus",
        'pdf_x_par'    : [0.0,1.0,-99.],
        'pdf_smear'    :"Gaus", 
        'pdf_smear_par': [0.05,0.01,-99.],
        'postfix'      : "gaus_5pc",
        'legend'       : "Gauss #pm5%, N=500",
        },
    'exp' : {
        'n_bins'       :  5,
        'x_low'        :  1.0,
        'x_high'       :  5.0,
        'n_events'     : 500,
        'n_toys'       : 1000,
        'pdf_x'        : "Exp",
        'pdf_x_par'    : [+1.0,-99,-99.],
        'pdf_smear'    :"Gaus", 
        'pdf_smear_par': [0.05,0.01,-99.],
        'postfix'      : "exp_5pc",
        'legend'       : "Exp #pm5%, N=50000",
        },    
}


cfg = {}
count = 0
for n_ev in [10, 
             25, 
             50, 
             100,
             500 
             ]:
    for n_b in [5,  
                #10
                ]:
        for pdf in ['gaus', 
                    #'exp'
                    ]:
            c = cfg_base[pdf].copy()
            c['n_toys'] = 200
            c['n_events'] = n_ev
            c['n_bins'] = n_b
            c['legend'] = (pdf+" #pm5%, N="+str(n_ev))
            c['postfix'] = 'TEST'
            cfg[pdf+'_'+str(count)] = c
            count += 1
        

for key, value in cfg.iteritems():

    #if key!='exp_1':
    #    continue

    create = True
    if create:
        print "Create tree with toy events..."
        create_tree(n_bins=value['n_bins'], x_low=value['x_low'], x_high=value['x_high'], 
                    n_events=value['n_events'], n_toys=value['n_toys'], 
                    pdf_x=value['pdf_x'], pdf_x_par=value['pdf_x_par'],
                    pdf_smear=value['pdf_smear'], pdf_smear_par=value['pdf_smear_par'],
                    postfix=value['postfix'])
        print "Run minimzer on toy events..."
        test_minimizer(fin_name='gen_toy_'+str(value['n_bins'])+'bin_'+str(value['n_events'])+'ev_'+value['postfix'], fout_name='result_'+str(value['n_bins'])+'bin_'+str(value['n_events'])+'ev_'+value['postfix'])

    print "Make plots..."
    make_plots(fin_name=str(value['n_bins'])+'bin_'+str(value['n_events'])+'ev_'+value['postfix'], n_bins=value['n_bins'], legend=value['legend'])
    #raw_input()
