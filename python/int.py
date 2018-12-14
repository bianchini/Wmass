import ROOT
import math
from sys import argv

class FixPoint:
    def __init__(self, x_low=0.92, x_high=1.08, nbins=201, p0=0.5, p1=1.0, p3=3.0):
        print 'Initialize'
        self.out = ROOT.TFile('fixpoint_'+argv[2]+'.root', 'RECREATE')
        self.histo = ROOT.TH1D('histo', 'Toys', nbins, x_low, x_high)
        self.p0 = p0
        self.p1 = p1
        self.p3 = p3
        self.histo.Sumw2()

    def pdf_gamma(self, gamma):
        pdf = math.pow((gamma-self.p0),2.0)*math.exp(-self.p1*(gamma-self.p0)) if (gamma<self.p3 and gamma>=1.0) else 0.0 
        #pdf = 1.0 if (gamma<self.p3 and gamma>=1.0) else 0.0 
        return pdf

    def pdf_y(self, y, width=0.01):
        return 1./math.pi/width/2./((y-1.)*(y-1.) + width*width/4)

    def gammabeta2(self,gamma):
        return gamma*gamma-1
    def gammabeta(self,gamma):
        return math.sqrt(self.gbeta2(gamma))
    def beta(self,gamma):
        return self.gbeta(gamma)/gamma

    def A0(self, gamma):
        #return 2./3.
        return math.log(gamma)*0.1
        #return math.tanh(gamma-1)*0.2
        #return float(argv[4])

    def A4(self, gamma):
        return math.log(gamma)
        #return float(argv[5])

    def fill_histo(self, ntoys=100, width=0.):
        step = (self.p3-1.0)/ntoys
        for ig in range(ntoys):                
            if ig%1000==0:
                print 'Toy ', ig ,'/', ntoys
            gamma = 1.0 + step*(ig+0.5)
            B0 = 3./8.*(1+0.5*self.A0(gamma))/math.sqrt(gamma*gamma-1)
            B1 = 3./8.*(self.A4(gamma))/(gamma*gamma-1)
            B2 = 3./8.*(1-3./2.*self.A0(gamma))/(gamma*gamma-1)/math.sqrt(gamma*gamma-1)        
            if width==0:
                for ix in range(self.histo.GetNbinsX()):
                    x = self.histo.GetBinCenter(ix+1)            
                    if gamma >= 0.5*(x+1./x):
                        integrand = self.pdf_gamma(gamma=gamma) *(B0 + B1*(x-gamma) + B2*(x-gamma)*(x-gamma))
                        self.histo.Fill(x, integrand*step)
            else:
                n_points_y = 12
                step_y = 6*width/n_points_y
                for ix in range(self.histo.GetNbinsX()):
                    x = self.histo.GetBinCenter(ix+1) 
                    for iy in range(n_points_y):
                        y = -3*width + iy*6*width/n_points_y + 1.0 
                        z = x/y
                        if gamma >= 0.5*(z+1./z):
                            integrand = self.pdf_gamma(gamma=gamma) *(B0 + B1*(z-gamma) + B2*(z-gamma)*(z-gamma))*self.pdf_y(y, width=width)/y
                            self.histo.Fill(x, integrand*step*step_y)
                        #self.histo.Fill(x, integrand*step*step_y)

        print 'End of loop'
        self.out.cd()
        self.histo.Write('', ROOT.TObject.kOverwrite)
        self.out.Close()

fixpoint = FixPoint(x_low=0.5, x_high=1.5, nbins=101, p0=float(argv[3]), p1=1.0, p3=4)
fixpoint.fill_histo(ntoys=int(argv[1]), width=0.0)

argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(False)
argv.remove( '-b-' )

f = ROOT.TFile.Open('fixpoint_'+argv[2]+'.root', 'READ')
c = ROOT.TCanvas("c", "canvas", 1200, 400)
c.Divide(3,1)
c.cd(1)
ROOT.gPad.SetLeftMargin(0.2)
leg = ROOT.TLegend(0.30,0.75,0.35,0.88, "","brNDC")
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)
leg.SetFillColor(10)
h = f.Get('histo')
h.SetTitle('f(x)')
h.SetStats(0)
h.SetXTitle('x')
h.SetYTitle('A.U.')
h.Draw('HISTE')
if len(argv)>4:
    leg.SetHeader('#splitline{g^{(0)}'+('>0' if float(argv[3])<1 else '=0')+'}{A_{0}='+'{:0.1f}'.format(float(argv[4]))+', A_{4}='+'{:0.1f}'.format(float(argv[5]))+'}')
else:
    leg.SetHeader('#splitline{g^{(0)}'+('>0' if float(argv[3])<1 else '=0')+'}{A_{0}=0.1log(#gamma), A_{4}=log(#gamma)}')
leg.Draw()

hD1= h.Clone('hD1')
hD1.Reset()
hD1.SetStats(0)
hD1.SetTitle('df/dx')
hD1.SetXTitle('x')
hD1.SetYTitle('A.U.')
hD1.SetLineColor(ROOT.kRed)
hD2 = h.Clone('hD2')
hD2.Reset()
hD2.SetStats(0)
hD2.SetTitle('d^{2}f/dx^{2}')
hD2.SetXTitle('x')
hD2.SetYTitle('A.U.')
hD2.SetLineColor(ROOT.kBlue)
for ib in range(2, h.GetNbinsX()):
    val0 = h.GetBinContent(ib)
    valM = h.GetBinContent(ib-1)
    valP = h.GetBinContent(ib+1)    
    d = h.GetBinWidth(ib)
    f_prime = (valP-valM)/d/2.0
    f_second = (valP+valM-2.0*val0)/d/d 
    hD1.SetBinContent(ib, f_prime  )
    hD2.SetBinContent(ib, f_second  )

c.cd(2)
ROOT.gPad.SetMargin(0.15,0.15,0.1,0.1)
hD1.Draw('HISTE')
c.cd(3)
ROOT.gPad.SetRightMargin(0.2)
hD2.Draw('HISTE')

c.SaveAs('toy_spectrum_'+argv[2]+'.png')
c.SaveAs('toy_spectrum_'+argv[2]+'.pdf')

f.Close()
#raw_input()
