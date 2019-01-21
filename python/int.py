import ROOT
import math
from sys import argv

class FixPoint:
    def __init__(self, x_low=0.92, x_high=1.08, nbins=201, p0=0.5, p1=1.0, p2=0.0, pR=3.0, pdf_type='power', create_file=True):
        print 'Initialize'
        if create_file:
            self.out = ROOT.TFile('fixpoint_'+argv[2]+'.root', 'RECREATE')
        self.histo = ROOT.TH1D('histo', 'Toys', nbins, x_low, x_high)
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2
        self.pR = pR
        self.histo.Sumw2()
        self.pdf_type = pdf_type

    def pdf_gamma(self, gamma):
        if self.pdf_type=='exp':
            pdf = (gamma-self.p0)*math.exp(-2*self.p1*(gamma-self.p0)) if (gamma<self.pR and gamma>=1.0) else 0.0 
        elif self.pdf_type=='power':
            pdf = math.pow(gamma-self.p0,self.p1) + self.p2 if (gamma<self.pR and gamma>=1.00) else 0.0 
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
        #return math.tanh(gamma-1)*0.2
        if len(argv)>7:
            return float(argv[7])
        else:
            #return math.log(gamma)*0.1
            return 1.*math.tanh(4*(gamma-1))

    def A4(self, gamma):
        #return math.log(gamma)
        if len(argv)>8:
            return float(argv[8])
        else:
            return 1*math.tanh(0.3*(gamma-1))


    def fill_histo(self, ntoys=100, width=0.):
        step = (self.pR-1.0)/ntoys
        for ig in range(ntoys):                
            if ig%1000==0:
                print 'Toy ', ig ,'/', ntoys
            gamma = 1.0 + step*(ig+0.5)
            B0 = 3./8.*(1+0.5*self.A0(gamma))/math.sqrt(gamma*gamma-1)
            B1 = 3./8.*(self.A4(gamma))/(gamma*gamma-1)
            B2 = 3./8.*(1-3./2.*self.A0(gamma))/(gamma*gamma-1)/math.sqrt(gamma*gamma-1)        
            B3 = 3./8.*1.0/(gamma*gamma-1)*(gamma*gamma-1)
            if width==0:
                for ix in range(self.histo.GetNbinsX()):
                    x = self.histo.GetBinCenter(ix+1)            
                    if gamma >= 0.5*(x+1./x):
                        integrand = self.pdf_gamma(gamma=gamma) *(B0 + B1*(x-gamma) + B2*(x-gamma)*(x-gamma)) #+ B3*(x-gamma)*(x-gamma)*(x-gamma)
                        self.histo.Fill(x, integrand*step)
            else:
                n_points_y = 20
                step_y = 6*width/n_points_y
                for ix in range(self.histo.GetNbinsX()):
                    x = self.histo.GetBinCenter(ix+1) 
                    for iy in range(n_points_y):
                        y = -3*width + iy*6*width/n_points_y + 1.0 
                        z = x/y                        
                        #integrand = self.pdf_y(y, width=width)
                        #self.histo.Fill(y, integrand*step*step_y)
                        #continue
                        if gamma >= 0.5*(z+1./z):
                            integrand = self.pdf_gamma(gamma=gamma) *(B0 + B1*(z-gamma) + B2*(z-gamma)*(z-gamma))*self.pdf_y(y, width=width)/y
                            self.histo.Fill(x, integrand*step*step_y)

        print 'End of loop'
        self.out.cd()
        self.histo.Write('', ROOT.TObject.kOverwrite)
        self.out.Close()

fixpoint = FixPoint(x_low=0.5, x_high=1.5, nbins=101, p0=float(argv[3]), p1=float(argv[4]), p2=float(argv[5]), pR=float(argv[6]), pdf_type='power', create_file=True)
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
ROOT.gPad.SetGrid()
leg = ROOT.TLegend(0.25,0.73,0.85,0.89, "","brNDC")
leg.SetFillStyle(1001)
leg.SetBorderSize(1)
leg.SetTextSize(0.05)
leg.SetFillColor(0)
h = f.Get('histo')
h.SetLineColor(ROOT.kBlack)
h.SetMaximum(h.GetMaximum()*1.25)
h.SetTitle('f^{(0)}(x)')
h.SetLineWidth(2)
h.SetStats(0)
h.SetXTitle('x')
h.SetTitleSize(0.05, 'X')
h.SetTitleSize(0.05, 'Y')
h.SetLabelSize(0.05, 'X')
h.SetLabelSize(0.05, 'Y')
h.SetNdivisions(5,'X')
h.SetNdivisions(10,'Y')
h.SetYTitle('A.U.')
h.Draw('HISTE')
if len(argv)>7:    
    if fixpoint.pdf_type=='exp':
        leg.SetHeader('#splitline{g(#gamma)=(#gamma-'+'{:0.1f}'.format(fixpoint.p0)+')e^{[-'+'{:0.1f}'.format(fixpoint.p1)+'(#gamma-'+'{:0.1f}'.format(fixpoint.p0)+')]}'+', #gamma#leq'+'{:0.1f}'.format(fixpoint.pR)+'}{A_{0}='+('{:0.1f}'.format(float(argv[7])) if float(argv[7])!=0.6667 else '2/3')+', A_{4}='+'{:0.1f}'.format(float(argv[8]))+'}')
    elif fixpoint.pdf_type=='power':
        leg.SetHeader('#splitline{g(#gamma)=(#gamma-'+'{:0.1f}'.format(fixpoint.p0)+')^{'+'{:0.1f}'.format(fixpoint.p1)+'}'+', '+'#gamma#leq'+'{:0.1f}'.format(fixpoint.pR)+'}{A_{0}='+('{:0.1f}'.format(float(argv[7])) if float(argv[7])!=0.6667 else '2/3')+', A_{4}='+'{:0.1f}'.format(float(argv[8]))+'}')
else:
    if fixpoint.pdf_type=='exp':
        leg.SetHeader('#splitline{g(#gamma)=(#gamma-'+'{:0.1f}'.format(fixpoint.p0)+')e^{[-'+'{:0.1f}'.format(fixpoint.p1)+'(#gamma-'+'{:0.1f}'.format(fixpoint.p0)+')]}'+', #gamma#leq'+'{:0.1f}'.format(fixpoint.pR)+'}{A_{0}=th[4(#gamma-1)], A_{4}=th[0.3(#gamma-1)]}')
    elif fixpoint.pdf_type=='power':
        leg.SetHeader('#splitline{g(#gamma)=(#gamma-'+'{:0.1f}'.format(fixpoint.p0)+')^{'+'{:0.1f}'.format(fixpoint.p1)+'}, #gamma#leq'+'{:0.1f}'.format(fixpoint.pR)+'}{A_{0}=th[4(#gamma-1)], A_{4}=th[0.3(#gamma-1)]}')
        
leg.Draw()

hD1=ROOT.TH1D('hD1','', h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
hD1.Reset()
hD1.SetStats(0)
hD1.SetLineWidth(2)
hD1.SetTitle('f^{(1)}(x)')
hD1.SetTitleSize(0.05, 'X')
hD1.SetTitleSize(0.05, 'Y')
hD1.SetLabelSize(0.05, 'X')
hD1.SetLabelSize(0.05, 'Y')
hD1.SetNdivisions(5,'X')
hD1.SetNdivisions(10,'y')
hD1.SetXTitle('x')
#hD1.SetYTitle('A.U.')
hD1.SetYTitle('')
hD1.SetLineColor(ROOT.kRed)
hD2=ROOT.TH1D('hD2','', h.GetNbinsX(), h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
hD2.Reset()
hD2.SetStats(0)
hD2.SetLineWidth(2)
hD2.SetTitle('f^{(2)}(x)')
hD2.SetXTitle('x')
#hD2.SetYTitle('A.U.')
hD2.SetYTitle('')
hD2.SetTitleSize(0.05, 'X')
hD2.SetTitleSize(0.05, 'Y')
hD2.SetLabelSize(0.05, 'X')
hD2.SetLabelSize(0.05, 'Y')
hD2.SetNdivisions(5,'X')
hD2.SetNdivisions(10,'y')
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
hD1.SetBinContent(1, hD1.GetBinContent(2))
hD1.SetBinContent(h.GetNbinsX(), hD1.GetBinContent(h.GetNbinsX()-1))
hD2.SetBinContent(1, hD2.GetBinContent(2))
hD2.SetBinContent(h.GetNbinsX(), hD2.GetBinContent(h.GetNbinsX()-1))

c.cd(2)
ROOT.gPad.SetMargin(0.15,0.15,0.1,0.1)
ROOT.gPad.SetGrid()
max_y = hD1.GetMaximum()
min_y = hD1.GetMinimum()
hD1.SetMaximum((max_y+min_y)*0.5 + 1.2*(max_y-min_y)*0.5)
hD1.SetMinimum((max_y+min_y)*0.5 - 1.2*(max_y-min_y)*0.5)
hD1.Draw('HIST')
c.cd(3)
ROOT.gPad.SetRightMargin(0.2)
ROOT.gPad.SetGrid()
max_y = hD2.GetMaximum()
min_y = hD2.GetMinimum()
hD2.SetMaximum((max_y+min_y)*0.5 + 1.2*(max_y-min_y)*0.5)
hD2.SetMinimum((max_y+min_y)*0.5 - 1.2*(max_y-min_y)*0.5)
hD2.Draw('HIST')

c.SaveAs('toy_spectrum_'+argv[2]+'.png')
c.SaveAs('toy_spectrum_'+argv[2]+'.pdf')

f.Close()
raw_input()
