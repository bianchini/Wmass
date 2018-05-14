import ROOT
import math
import sys
import os
import pickle
import time
from pprint import pprint
import numpy as np
import numpy.random as ran

from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from array import array;
import copy


class TemplateFitter:

    def __init__(self, 
                 DY='CC_FxFx', 
                 charge='Wplus', 
                 var='WpreFSR', 
                 job_name='TEST', 
                 mc_mass=80.419, 
                 verbose=True, 
                 fixed_parameters=[],
                 ):

        self.in_dir  = os.environ['CMSSW_BASE']+'/src/Wmass/test/plots/'
        self.out_dir = os.environ['CMSSW_BASE']+'/src/Wmass/test/'
        self.job_name = job_name
        self.verbose = verbose
        self.fixed_parameters = fixed_parameters
        self.mc_mass = mc_mass
        self.map_params = {}

        self.res_coeff = np.load(open(self.in_dir+'fit_results_'+DY+'_'+charge+'_'+var+'_all_A0-4_forced.pkl', 'r'))        
        self.cov_coeff = np.load(open(self.in_dir+'covariance_'+DY+'_'+charge+'_'+var+'_stat_plus_syst_all_A0-4_forced.npy', 'r'))        

        templates = np.load(self.in_dir+'template_'+charge+'_'+var+'_val.npz')
        templates_files = {'template': 0, 'masses': 1, 'bins_qt' : 2, 
                           'bins_y'  : 3, 'coefficients' : 4, 'bins_eta': 5, 
                           'bins_pt' : 6, 'mc_acceptances' : 7 }

        # for template, save the size of (pt,eta) plane 
        # for bins, save len(bins)-1, since bins are edges
        # for coeff, save len(coeff)-2, since coeff[-1] = MC and coeff[-2] = UL
        for key,p in templates_files.items():             
            template_file = templates['arr_'+str(p)]
            print ('Adding file arr_'+str(p)+' with shape...'), template_file.shape, (' to self with key... '+key)
            size = -1
            if key=='template':
                size = template_file[0][0][0][0].size
                setattr(self, key,template_file)
            elif key in ['bins_qt', 'bins_y']:
                #size = template_file.size - 2
                #setattr(self, key,template_file[:-1])
                size = 1
                setattr(self, key,template_file[:2])
            elif key in ['bins_pt', 'bins_eta']:
                size = template_file.size - 1
                setattr(self, key,template_file)
            elif key=='coefficients':
                size = template_file.size - 2
                setattr(self, key,template_file[:-2])
            else:
                size = template_file.size
                setattr(self, key,template_file)
            setattr(self, key+'_size', size)

        self.mc_mass_index = np.where(self.masses==mc_mass)[0][0]
        self.mc = self.template.sum(axis=(1,2))[self.mc_mass_index,-1]
        self.save_template_snapshot(data=self.mc, title='MC', tag='mc')

        self.data = copy.deepcopy(self.mc)
        self.ndof = self.data.size
        self.save_template_snapshot(data=self.mc, title='Data', tag='data')

        self.fit_results = {}
        self.out_file_res = open(self.out_dir+'result_'+DY+'_'+charge+'_'+job_name+'.pkl','wb')
        #pickle.dump(self.fit_results, self.out_file_res)
        #self.out_file_res.close()

        self.dim_A = 0
        self.dim_beta = 0
        self.n_param = 0
        self.book_parameters()
        self.build_aux_matrices()

    def define_parameter( self, par_name='', start=0., step=0., par_range=() ):
        self.gMinuit.DefineParameter( self.n_param, par_name, start, step, par_range[0], par_range[1] )        
        self.map_params[par_name] = self.n_param
        match = False
        for fix in self.fixed_parameters:
            if fix in par_name:
                match = True
        if (par_name in self.fixed_parameters) or \
                ('OF' in par_name) or \
                match:
            self.gMinuit.FixParameter(self.n_param)
        else:
            self.ndof -= 1
        self.n_param += 1


    def book_parameters(self):

        self.gMinuit = ROOT.TMinuit(500)  
        self.gMinuit.SetPrintLevel(1)
        if not self.verbose:
            self.gMinuit.SetPrintLevel(-1)
        self.gMinuit.SetFCN( self.fcn ) 

        self.arglist = array( 'd', 10*[0.] )
        self.ierflg = ROOT.Long(0)
        
        # 0.5: chi2, 1.0: -2nll
        self.arglist[0] = 1.0
        self.gMinuit.mnexcm( "SET ERR", self.arglist, 1, self.ierflg )

        self.define_parameter( par_name='mass', start=self.mc_mass, step=0.020, par_range=(self.mc_mass-0.500,self.mc_mass+0.500) )
        
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for iqt in range(self.bins_qt_size):
                qt_bin = self.get_qt_bin(iqt)
                mc_norm = self.template[self.mc_mass_index][iqt][iy][-1].sum()/self.mc_acceptances[self.mc_mass_index][iqt][iy]
                self.define_parameter( par_name=y_bin+'_'+qt_bin+'_norm', 
                                       start=mc_norm, 
                                       step=math.sqrt(mc_norm), 
                                       par_range=(max(mc_norm-math.sqrt(mc_norm)*10., 0.), mc_norm+math.sqrt(mc_norm)*10.) )
                for coeff in self.coefficients:
                    self.dim_A += 1
                    if iqt>0:
                        continue
                    (valid_orders, order) = self.get_orders(coeff, y_bin)
                    for o in valid_orders:
                        self.define_parameter( par_name=y_bin+'_'+coeff+'_'+'pol'+str(order)+'_p'+str(o), 
                                               start=self.res_coeff[coeff+'_'+y_bin+'_fit'][o], 
                                               step=0.01, 
                                               par_range=(-2., +2.) )
                        self.dim_beta += 1


    def get_orders(self, coeff='', y_bin=''):               
        valid_orders = []
        order = len(self.res_coeff[coeff+'_'+y_bin+'_fit'])-1
        for io,o in enumerate(self.res_coeff[coeff+'_'+y_bin+'_fit']):
            if o!=0.:
                valid_orders.append(io)
        return (valid_orders, order)                

    def get_y_bin(self, iy=-1):
        return 'y{:03.2f}'.format(self.bins_y[iy])+'_'+'y{:03.2f}'.format(self.bins_y[iy+1])

    def get_qt_bin(self, iqt=-1):
        return 'qt{:03.1f}'.format(self.bins_qt[iqt])+'_'+'qt{:03.1f}'.format(self.bins_qt[iqt+1])

    def get_y_qt_bin(self, iy=-1, iqt=-1):
        return self.get_y_bin(iy)+'_'+self.get_qt_bin(iqt)

    def build_aux_matrices(self):

        self.n = self.data.flatten()
        if self.n[self.n<10].size > 0:
            print "Bins with less than 10 entries!"
        self.Vinv = np.diag(1./self.n)
 
        self.K = np.zeros( ( self.dim_A, self.dim_beta) )

        # first loop
        idx_A = 0
        for iy1 in range(self.bins_y_size):
            for iqt in range(self.bins_qt_size):
                for coeff1 in self.coefficients:

                    # second loop
                    idx_beta = 0
                    for iy2 in range(self.bins_y_size):                    
                        y_bin2 = self.get_y_bin( iy2 )
                        for coeff2 in self.coefficients:
                            (valid_orders, order) = self.get_orders(coeff2, y_bin2)
                            for o in valid_orders:
                                if iy1==iy2 and coeff1==coeff2:
                                    self.K[idx_A, idx_beta] += math.pow( 0.5*(self.bins_qt[iqt]+self.bins_qt[iqt+1]), o)
                                idx_beta += 1
                idx_A += 1

        #np.set_printoptions(threshold=np.inf)
        #print self.K[0]
        #print self.K[1]

    def fcn(self, npar, gin, f, par, iflag ):
        #nll = math.pow((par[0]-self.mc_mass)/0.100, 2.0)
        nll = self.chi2(par=par)
        print 'nll: ', nll/self.ndof
        f[0] = nll
        return

    def chi2(self, par):
        Am = np.zeros((self.data.size, self.dim_A))
        bm = np.zeros(self.n.shape)

        idx_A = 0
        for iy in range(self.bins_y_size):
            for iqt in range(self.bins_qt_size):
                inorm = par[self.map_params[self.get_y_qt_bin(iy,iqt)+'_norm']]
                tUL = copy.deepcopy(self.template[0][iqt][iy][-2])
                bm += tUL.flatten()*inorm
                print 'UL*',inorm
                for icoeff,coeff in enumerate(self.coefficients):
                    tjA = copy.deepcopy(self.template[0][iqt][iy][icoeff])
                    tjA -= tUL
                    print coeff+'*', inorm
                    Am[:, idx_A] += tjA.flatten()*inorm
                    idx_A += 1

        np.set_printoptions(threshold=np.inf)
        print Am
        print bm

        nprime = self.n - bm
        aux1 =  np.linalg.multi_dot([self.K.T, Am.T, self.Vinv, Am, self.K])        
        aux1_inv = np.linalg.inv(aux1)
        aux2 = np.linalg.multi_dot( [ self.K.T, Am.T, self.Vinv, nprime ] )
        beta = np.linalg.multi_dot( [aux1_inv, aux2 ] )        
        res1 = (nprime - np.linalg.multi_dot( [Am, self.K, beta]) )
        res2 = (beta_prior - beta)
        chi2min = np.linalg.multi_dot( [res1.T, self.Vinv, res1 ] )
        return chi2min
    

    def run(self, n_points=500000, run_minos=False):

        self.arglist[0] = n_points
        self.arglist[1] = 1.0
        print("Convergence at EDM %s" % (self.arglist[1]*0.001))
        
        self.gMinuit.mnexcm( "HES", self.arglist, 1, self.ierflg )
        self.arglist[0] = 1
        self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
        status = self.gMinuit.Migrad()
        if status>0:
            self.arglist[0] = 2
            self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
            status = self.gMinuit.Migrad()

        if run_minos:
            massL, massH =  ROOT.Double(0.), ROOT.Double(0.) 
            self.gMinuit.mnmnot(1, 1, massH, massL)

        amin, edm, errdef    = ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.)
        nvpar, nparx, icstat = ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985)
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        if self.verbose:
            self.gMinuit.mnprin( 3, amin )

        self.cov = ROOT.TMatrixDSym(self.n_param)
        self.gMinuit.mnemat(self.cov.GetMatrixArray(), self.n_param)
        self.plot_cov_matrix()

    def save_template_snapshot(self, data=np.array([]), title='', tag=''):
        xx,yy = np.meshgrid(self.bins_pt, self.bins_eta)        
        plt.pcolormesh(yy, xx, data)
        plt.colorbar(format='%.0e')
        plt.axis([self.bins_eta.min(), self.bins_eta.max(), self.bins_pt.min(), self.bins_pt.max()])
        plt.title(title)
        plt.show()
        plt.savefig(self.out_dir+'snapshot_'+tag+'_'+self.job_name+'.png')
        plt.close('all')

    def plot_cov_matrix(self):
        c = ROOT.TCanvas("canvas", "canvas", 600, 600) 
        n_free = self.gMinuit.GetNumFreePars()
        h2 = ROOT.TH2D('cov', '', n_free, 0, n_free, n_free, 0, n_free)
        h2.SetStats(0) 
        for i in range(n_free):
            for j in range(n_free):
                rho_ij = self.cov(i,j)/math.sqrt(self.cov(i,i)*self.cov(j,j)) if self.cov(i,i)>0.0 and self.cov(j,j)>0.0 else 0.0
                h2.SetBinContent(i+1, j+1, rho_ij )
        h2.Draw("COLZ")
        c.SaveAs(self.out_dir+'covariance_'+self.job_name+'.png')
        c.SaveAs(self.out_dir+'covariance_'+self.job_name+'.C')
        c.IsA().Destructor( c )






