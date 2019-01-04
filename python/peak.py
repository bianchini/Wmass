import ROOT
import math
import sys
import os
from pprint import pprint
import numpy as np
import numpy.random as ran

from scipy import optimize
from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import rcParams
rcParams['figure.figsize'] = 8,7

from copy import copy

scale_lumi = 1

def polynomial(x, coeff, order, der=0):
    if type(x)!=np.ndarray:
        xx = np.array([x])
    else:
        xx = x
    val = np.zeros(xx.size)
    for p in range(der, order+1):
        val +=  math.factorial(p)/math.factorial(p-der)*coeff[p]*np.power(xx,p-der)
    return val

class Peak:

    def __init__(self, f_name='out.root', var=['e'], cut='all', M=80.419, charges='both', derivables=[1,3], tag='', verbose=True ):
        self.f = ROOT.TFile.Open(f_name)
        self.cut = cut
        self.var = var
        self.charges=charges
        self.derivables = derivables
        self.tag = tag
        self.M = M
        self.verbose=verbose

    def initialize_central(self, deg=2, addbins=2, Wtype='preFSR'):
        if self.charges=='toy':
            h_plus  = self.f.Get('histo').Clone('h_plus')
        elif self.charges!='both':
            h_plus  = self.f.Get('W'+Wtype+'_mu_'+self.var[0]+'_'+self.charges+'_'+self.cut+'_w0').Clone('h_plus')
            if len(self.var)>1:
                h_plus.Add(  self.f.Get('W'+Wtype+'_mu_'+self.var[1]+'_'+self.charges+'_'+self.cut+'_w0') )
        else:
            h_plus  = self.f.Get('W'+Wtype+'_mu_'+self.var[0]+'_plus_'+self.cut+'_w0').Clone('h_plus')
            h_minus = self.f.Get('W'+Wtype+'_mu_'+self.var[0]+'_minus_'+self.cut+'_w0').Clone('h_minus')
            h_plus.Add(h_minus)    
            if len(self.var)>1:
                h_plus.Add(  self.f.Get('W'+Wtype+'_mu_'+self.var[1]+'_plus_'+self.cut+'_w0') )
                h_plus.Add(  self.f.Get('W'+Wtype+'_mu_'+self.var[1]+'_minus_'+self.cut+'_w0') )

        (self.nbins, self.x_min, self.x_max, self.bin_size, self.ib_max) =  (h_plus.GetNbinsX(), h_plus.GetXaxis().GetXmin(), h_plus.GetXaxis().GetXmax(), h_plus.GetBinWidth(1)/self.M, h_plus.FindBin(self.M*0.5) )
        self.ib_low  = self.ib_max - deg/2 - deg%2 - addbins
        self.ib_high = self.ib_max + deg/2  + addbins
        self.x_min_fit = h_plus.GetBinLowEdge(self.ib_low)
        self.x_max_fit = h_plus.GetBinLowEdge(self.ib_high+1)
        self.Wtype = Wtype
        if self.verbose:
            print 'Fitting ', self.ib_high-self.ib_low+1, 'bins in the range '+self.var[0]+' in [', self.x_min_fit, ',',  self.x_max_fit , ']' 

    def get_peak(self, deg=2, addbins=2, der=0, weight=0, ntoys=-1, do_plot=True, do_mean=False, title='', do_random_smearing=-1):        

        fit_tag = 'der'+str(deg)+'_addbins'+str(addbins)+'_w'+str(weight)
        str_weight = str(weight)
        if do_random_smearing>0:
            str_weight = '0'

        if self.charges=='toy':
            h_plus  = self.f.Get('histo').Clone('h_plus')
        elif self.charges!='both':
            h_plus  = self.f.Get('W'+self.Wtype+'_mu_'+self.var[0]+'_'+self.charges+'_'+self.cut+'_w'+str_weight).Clone('h_plus')
            if len(self.var)>1:
                h_plus.Add(  self.f.Get('W'+self.Wtype+'_mu_'+self.var[1]+'_'+self.charges+'_'+self.cut+'_w'+str_weight))            
        else:
            h_plus = self.f.Get('W'+self.Wtype+'_mu_'+self.var[0]+'_plus_'+self.cut+'_w'+str_weight).Clone('h_plus')
            h_minus = self.f.Get('W'+self.Wtype+'_mu_'+self.var[0]+'_minus_'+self.cut+'_w'+str_weight).Clone('h_minus')
            h_plus.Add(h_minus)
            if len(self.var)>1:
                h_plus.Add(  self.f.Get('W'+self.Wtype+'_mu_'+self.var[1]+'_plus_'+self.cut+'_w'+str_weight) )
                h_plus.Add(  self.f.Get('W'+self.Wtype+'_mu_'+self.var[1]+'_minus_'+self.cut+'_w'+str_weight) )

        x     = np.array( [h_plus.GetBinCenter(ib)/self.M*2.0 - 1.0  for ib in range(self.ib_low, self.ib_high+1) ] )
        y     = np.array( [h_plus.GetBinContent(ib)/h_plus.GetBinContent(self.ib_max) for ib in range(self.ib_low, self.ib_high+1) ] )        
        w     = np.array( [h_plus.GetBinContent(ib) for ib in range(self.ib_low, self.ib_high+1) ] )
        w2    = np.array( [h_plus.GetBinError(ib)*h_plus.GetBinError(ib) for ib in range(self.ib_low, self.ib_high+1) ] )
        y_err = np.array( [h_plus.GetBinError(ib)/h_plus.GetBinContent(self.ib_max)/math.sqrt(scale_lumi) for ib in range(self.ib_low, self.ib_high+1) ] )

        if do_random_smearing>0:
            y += np.random.normal( 0., do_random_smearing*y_err )
        
        self.mean_val = (x*y).sum()/np.sum(y)
        self.mean_err = math.sqrt((x*x*y).sum()/np.sum(y) - math.pow(self.mean_val,2))/math.sqrt( (w.sum()*w.sum())/w2.sum() )

        Vinv = np.diag(1./y_err/y_err)
        A = np.zeros((x.size,deg+1))
        for i in range(x.size):
            for j in range(deg+1):
                A[i,j] = math.pow((x[i]),j)
        self.cov = np.linalg.inv( np.linalg.multi_dot([A.T, Vinv, A]) )
        self.coeff = np.linalg.multi_dot([self.cov, A.T, Vinv, y])
        self.chi2 =  np.linalg.multi_dot([(np.linalg.multi_dot([A,self.coeff])-y).T,Vinv,np.linalg.multi_dot([A,self.coeff])-y])
        self.ndf = x.size-(deg+1)
        if ntoys>0 and self.verbose:
            print 'Chi2=', '{:0.2f}'.format(self.chi2), 'ndf=', self.ndf, 'chi2/ndf=', '{:0.3f}'.format(self.chi2/self.ndf)

        extremal = [-1.0]
        if der in self.derivables:
            try:
                extremal = optimize.newton( polynomial, 0.0, fprime=lambda x, coeff, deg, der : polynomial(x=x, coeff=coeff, order=deg, der=(der+1)), 
                                            args=(self.coeff, deg, der), tol=1.48e-08, maxiter=50, 
                                            fprime2=lambda x, coeff, deg, der : polynomial(x=x, coeff=self.coeff, order=deg, der=(der+2)) )
            except RuntimeError:
                print 'No root for','Deg:', deg, ', derivative:', der, ', add bins:', addbins
                extremal = [-1.0]

        if do_mean:
            self.points_err = [self.mean_err]
            extremal = [self.mean_val]       
            
        points = np.array( extremal )
        setattr(self, 'points', points)

        if ntoys<0 or do_mean:
            return

        y_rnd      = np.zeros( (x.size, ntoys) )                
        points_rnd = np.zeros( (1,ntoys) )    
        newton_errors = 0
        for itoy in range(ntoys):
            coeff_rnd = np.random.multivariate_normal(self.coeff, self.cov) 
            y_rnd[:, itoy] = polynomial(x=x, coeff=coeff_rnd, order=deg, der=der)
            if der in self.derivables:
                try:
                    point_rnd = optimize.newton( polynomial, 0.0, 
                                                 fprime=lambda x, coeff_rnd, deg, der : polynomial(x=x, coeff=coeff_rnd, order=deg, der=(der+1)), 
                                                 args=(coeff_rnd, deg, der), tol=1.48e-08, maxiter=50, 
                                                 fprime2=lambda x, coeff_rnd, deg, der : polynomial(x=x, coeff=coeff_rnd, order=deg, der=(der+2)))
                    points_rnd[0, itoy] = point_rnd[0]
                except RuntimeError:
                    newton_errors += 1
                    continue

        #print 'Fraction of errors: ', newton_errors, '/', ntoys, '=',  newton_errors/ntoys

        # get errors
        if der in self.derivables:
            points_err = np.std( points_rnd, axis=1)
            self.points_err = points_err
            if self.verbose:
                print 'Deg:', deg, ', derivative:', der, ', add bins:', addbins, ' => ', '{:0.5f}'.format(points[0]), '+/-', '{:0.5f}'.format(points_err[0]), \
                    ' => ', '{:0.5f}'.format((points[0]+1)*self.M), '+/-', '{:0.5f}'.format((points_err[0])*self.M)
        else:
            points_err = [-1.0]

        if not do_plot:
            return

        # plotting
        plt.figure()
        fig, ax = plt.subplots()
        ax.fill_between(x, polynomial(x=x, coeff=self.coeff, order=deg, der=der)-np.std(y_rnd, axis=1), 
                        polynomial(x=x, coeff=self.coeff, order=deg, der=der)+np.std(y_rnd, axis=1), 
                        color='y', linestyle='-', label=r'$\pm1\sigma$')
        ax.plot(x, polynomial(x=x, coeff=self.coeff, order=deg, der=der), 'r--', label=r'Fit ($\mathrm{pol}_{'+str(deg)+'}$), $\chi^2$/ndof='+'{:0.2f}'.format(self.chi2/self.ndf), linewidth=3.0)
        
        if der==0:
            ax.errorbar(x=x, y=y, xerr=[self.bin_size]*x.size, yerr=y_err, fmt='o', color='black', label='')
            #plt.axis( [x[0]-0.01, x[-1]+0.01, 0.96, 1.04] )

        ax.set_xlim( x[0]-0.01, x[-1]+0.01 )
        plt.grid(True)
        charge_label = '$f^{('+str(der)+')}$'
        if self.charges=='both':
            charge_label += (', $W^\pm$')
        else:
            charge_label += (', $W^+$' if self.charges=='plus' else '$W^-$')
        legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
        legend.set_title(charge_label+', '+title, prop = {'size':'x-large'})
        #plt.title('Polynom of deg '+str(deg)+', +/- '+str(addbins)+'bins from 0.0, derivative '+str(der))
        #plt.title(charge_label+', $d^{'+str(der+1)+'}\sigma/dx^{'+str(der+1)+'}$', fontsize='x-large')
        plt.xlabel('$E/E_{W}-1$', fontsize=15)
        #plt.ylabel('$d^{'+str(der+1)+'}\sigma/dx^{'+str(der+1)+'}$', fontsize=20)
        plt.ylabel('A.U.',  fontsize=15)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()
        plt.savefig('fit_pol'+str(deg)+'_add'+str(addbins)+'_der'+str(der)+'_'+self.charges+self.tag+'.png')
        plt.savefig('fit_pol'+str(deg)+'_add'+str(addbins)+'_der'+str(der)+'_'+self.charges+self.tag+'.eps')
        plt.close('all') 

    def run_fits(self, degs=[4], addbins=[5], do_systematics=False, weights_scale=[0], weights_pdf=[9], Wtype='preFSR', title=''):
        for deg in degs:
            for iaddbin,addbin in enumerate(addbins):
                self.initialize_central(deg=deg, addbins=addbin, Wtype=Wtype)
                for der in range(0, deg):
                    self.get_peak(deg=deg, addbins=addbin, der=der, weight=0, ntoys=1000, title=title)                    
                    nominal = self.points
                    if not do_systematics:
                        continue                    
                    extremals_scale = np.zeros(len(weights_scale))
                    if der not in self.derivables:
                        continue
                    for iw,w in enumerate(weights_scale):
                        self.get_peak(deg=deg, addbins=addbin, der=der, weight=w, ntoys=1000,  do_plot=False, do_random_smearing=-1)
                        extremals_scale[iw] = (self.points[0]-nominal[0])/self.points_err[0]
                    #scale_err = (np.max(extremals_scale)-np.min(extremals_scale))/2*self.M*1000
                    scale_err = np.mean(extremals_scale) , '+/-', np.std(extremals_scale)
                    print 'RMS (scale)', scale_err #, 'MeV'
                    extremals_pdf = np.zeros(len(weights_pdf))
                    for iw,w in enumerate(weights_pdf):
                        self.get_peak(deg=deg, addbins=addbin, der=der, weight=w, ntoys=1000,  do_plot=False, do_random_smearing=-1)
                        extremals_pdf[iw] = (self.points[0]-nominal[0])/self.points_err[0]
                    #pdf_err = np.std(extremals_pdf)*self.M*1000
                    pdf_err = np.mean(extremals_pdf) , '+/-', np.std(extremals_pdf)
                    print 'RMS (pdf)  ', pdf_err #, 'MeV'


    def run_fits_vs_mass(self, degs=[4], addbins=[5], do_systematics=False, weights_mass=[], Wtype='preFSR', title=''):
        for deg in degs:
            for iaddbin,addbin in enumerate(addbins):
                colors = ['magenta', 'red', 'yellow', 'blue']
                fmts=['o', '-', '+', '--']
                plt.figure()
                fig, ax = plt.subplots()
                self.initialize_central(deg=deg, addbins=addbin, Wtype=Wtype)
                for ider,der in enumerate(range(0, deg)):
                    do_mean=(True if der==0 else False)
                    self.get_peak(deg=deg, addbins=addbin, der=der, weight=len(weights_mass)/2, ntoys=1000, do_plot=False, do_mean=do_mean)
                    if not do_systematics:
                        continue                    
                    extremals_mass = np.zeros(len(weights_mass))
                    extremals_err_mass = np.zeros(len(weights_mass))
                    if der not in self.derivables:
                        continue
                    for iw,w in enumerate(weights_mass):
                        self.get_peak(deg=deg, addbins=addbin, der=der, weight=iw, ntoys=1000, do_plot=False, do_mean=do_mean)
                        extremals_mass[iw] = self.points[0]
                        extremals_err_mass[iw] = self.points_err[0]
                    lsq = np.linalg.lstsq( np.vstack([weights_mass, np.ones(len(weights_mass))]).T, (extremals_mass+1)*self.M, rcond=None )[0]
                    ax.errorbar(x=np.array(weights_mass)+ider*0.002, y=(extremals_mass+1)*self.M, xerr=0, yerr=extremals_err_mass*self.M, fmt=fmts[ider], color=colors[ider], linewidth=1,
                                label=(('$\hat{x}_{'+str(der)+'}$=' if der>0 else '$x_{\mu}$=')+'{:0.3f}'.format(lsq[0])+'$\cdot M_{W}$ '+('+' if lsq[1]>0 else '')+'{:0.1f}'.format(lsq[1])+' GeV')
                                )
                plt.grid(True)
                charge_label = '$W^+$' if self.charges=='plus' else '$W^-$'
                if self.charges=='both':
                    charge_label = '$W^\pm$'                    
                ax.plot(weights_mass, weights_mass, 'k--', label=r'Ideal', linewidth=3.0)
                legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
                legend.set_title(charge_label+','+title, prop = {'size':'x-large'})
                plt.axis([80.100, 80.750, 78.500, 80.800])
                plt.show()
                #plt.title(charge_label)
                plt.xlabel('$M_{W}$ (GeV)', fontsize=15)
                plt.ylabel('$\hat{x}_{i}$ (GeV)', fontsize=15)
                plt.xticks(fontsize=15)
                plt.yticks(fontsize=15)
                plt.savefig('calibration_'+self.var[0]+'_'+self.charges+self.tag+'.png')
                plt.savefig('calibration_'+self.var[0]+'_'+self.charges+self.tag+'.eps')
                plt.close('all') 


    def run_mean_vs_range(self, degs=[4], addbins=[5], do_systematics=False, weights_scale=[0], weights_pdf=[9], Wtype='preFSR'):
        ranges       = np.zeros(addbins.size)
        values_scale = np.zeros(addbins.size)
        values_pdf   = np.zeros(addbins.size)
        values_stat  = np.zeros(addbins.size)
        for deg in degs:
            for iaddbin,addbin in enumerate(addbins):
                self.initialize_central(deg=deg, addbins=addbin, Wtype=Wtype)
                ranges[iaddbin] = (self.x_max_fit-self.x_min_fit)*0.5
                self.get_peak(deg=deg, addbins=addbin, der=-1, weight=0, ntoys=-1)
                values_stat[iaddbin] = self.mean_err*self.M*1000
                extremals_scale = np.zeros(len(weights_scale))
                for iw,w in enumerate(weights_scale):
                    self.get_peak(deg=deg, addbins=addbin, der=-1, weight=w, ntoys=-1)
                    extremals_scale[iw] = self.mean
                scale_err = (np.max(extremals_scale)-np.min(extremals_scale))/2*self.M*1000
                print 'RMS (scale)', scale_err, 'MeV'
                values_scale[iaddbin] = scale_err
                extremals_pdf = np.zeros(len(weights_pdf))
                for iw,w in enumerate(weights_pdf):
                    self.get_peak(deg=deg, addbins=addbin, der=-1, weight=w, ntoys=-1)
                    extremals_pdf[iw] = self.mean
                pdf_err = np.std(extremals_pdf)*self.M*1000
                print 'RMS (pdf)  ', pdf_err, 'MeV'
                values_pdf[iaddbin] = pdf_err
            plt.figure()
            fig, ax = plt.subplots()
            ax.plot(ranges, values_stat, 'k-', label=r'Stat', linewidth=3.0)
            ax.plot(ranges, values_scale, 'r--', label=r'Scale', linewidth=3.0)
            ax.plot(ranges, values_pdf, 'b--', label=r'PDF', linewidth=3.0)
            plt.grid(True)
            legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
            plt.axis( [0, ranges[-1], 0, 100] )
            plt.show()
            plt.title(self.var[0])
            plt.xlabel('GeV from $M_W/2$', fontsize=20)
            plt.ylabel('$\Delta\mu$ (MeV)', fontsize=20)
            plt.savefig('systs_'+self.var[0]+'_'+self.charges+'.png')
            plt.close('all') 

    def close(self):
        self.f.Close()
        
########################

weights_scale = [0,1,2,3,4,6,8]
weights_pdf = range(9,109)
weights_mass = [80.119, 80.169, 80.219, 80.269, 80.319, 80.369, 80.419, 80.469, 80.519, 80.569, 80.619, 80.669, 80.719]
addbins  = np.array([38])

for charges in ['plus', 
                'minus', 
                'both']:
    #continue
    peak = Peak( f_name='out_pt_weights_acc.root',  var=['e'], cut='all', M=80.419, charges=charges, derivables=[1,3], tag='_acc', verbose=False)
    peak.run_fits(degs=[4], addbins=addbins, do_systematics=True, weights_scale=weights_scale, weights_pdf=weights_pdf, Wtype='preFSR', title='pre-FSR lepton')
    peak.close()

for charges in ['plus', 'minus', 'both']:
    continue
    peak = Peak( f_name='out_pt_weights_mass_acc.root',  var=['e'], cut='all', M=80.419, charges=charges, derivables=[0,1,3], tag='_acc')
    peak.run_fits_vs_mass(degs=[4], addbins=addbins, do_systematics=True, weights_mass=weights_mass, Wtype='preFSR', title='pre-FSR lepton in acceptance')
    peak.close()

    peak = Peak( f_name='out_pt_weights_mass_full.root',  var=['e'], cut='all', M=80.419, charges=charges, derivables=[0,1,3], tag='_full')
    peak.run_fits_vs_mass(degs=[4], addbins=addbins, do_systematics=True, weights_mass=weights_mass, Wtype='preFSR', title='pre-FSR lepton')
    peak.close()

    peak = Peak( f_name='out_pt_weights_mass_acc_bare.root',  var=['e'], cut='all', M=80.419, charges=charges, derivables=[0,1,3], tag='_acc_bare')
    peak.run_fits_vs_mass(degs=[4], addbins=addbins, do_systematics=True, weights_mass=weights_mass, Wtype='bare', title='bare lepton in acceptance')
    peak.close()

    peak = Peak( f_name='out_pt_weights_mass_full_bare.root',  var=['e'], cut='all', M=80.419, charges=charges, derivables=[0,1,3], tag='_full_bare')
    peak.run_fits_vs_mass(degs=[4], addbins=addbins, do_systematics=True, weights_mass=weights_mass, Wtype='bare', title='bare lepton')
    peak.close()
