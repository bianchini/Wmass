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

from scipy import stats

from array import array;
import copy


class TemplateFitter:

    def __init__(self, 
                 DY='CC_FxFx', 
                 charge='Wplus', 
                 var='WpreFSR', 
                 job_name='TEST', 
                 mc_mass=80.419, 
                 num_events=-1,
                 verbose=True, 
                 fixed_parameters=[],
                 use_prior=False,
                 reduce_qt=-1,
                 reduce_y=-1,
                 do_parametric=True,
                 use_prefit=False,
                 add_nonclosure=True
                 ):

        self.in_dir  = os.environ['CMSSW_BASE']+'/src/Wmass/test/plots/'
        self.out_dir = os.environ['CMSSW_BASE']+'/src/Wmass/test/'
        self.job_name = job_name
        self.num_events = num_events 
        
        self.verbose = verbose
        self.fixed_parameters = fixed_parameters
        self.use_prior = use_prior
        self.mc_mass = mc_mass
        self.do_parametric = do_parametric
        self.map_params = {}
        self.release_for_hesse = False

        self.res_coeff = np.load(open(self.in_dir+'fit_results_'+DY+'_'+charge+'_'+var+'_all_A0-4_forced.pkl', 'r'))        

        templates = np.load(self.in_dir+'template_'+charge+'_'+var+'_val_masses.npz')
        templates_files = {'template': 0, 'masses': 1, 'bins_qt' : 2, 
                           'bins_y'  : 3, 'coefficients' : 4, 'bins_eta': 5, 
                           'bins_pt' : 6, 'mc_acceptances' : 7 }

        # for template, save the size of (pt,eta) plane 
        # for bins, save len(bins)-1, since bins are edges
        # for coeff, save len(coeff)-2, since coeff[-1] = MC and coeff[-2] = UL
        for key,p in templates_files.items():             
            template_file = templates['arr_'+str(p)]
            if verbose:
                print ('Adding file arr_'+str(p)+' with shape...'), template_file.shape, (' to self with key... '+key)
            size = -1
            if key=='template':
                size = template_file[0][0][0][0].size
                setattr(self, key,template_file)
            elif key=='bins_qt':
                size = template_file.size - 1 + reduce_qt
                setattr(self, key, template_file[:reduce_qt] if reduce_qt<0 else template_file)
            elif key=='bins_y':
                size = template_file.size - 1 + reduce_y
                setattr(self, key, template_file[:reduce_y] if reduce_y<0 else template_file)
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

        print 'Fitting these bins:'
        print '\ty :', self.bins_y
        print '\tqt:', self.bins_qt

        # remove UL from all templates
        for icoeff,coeff in enumerate(self.coefficients):
            self.template[:,:,:,icoeff] -= self.template[:,:,:,-2]

        self.mc_mass_index = np.where(self.masses==mc_mass)[0][0]
        print 'MC mass at index...', self.mc_mass_index
        self.mc = copy.deepcopy(self.template.sum(axis=(1,2))[self.mc_mass_index,-1])
        self.ndof = self.mc.size        
        if add_nonclosure:
            print 'Removing non closure by adding Sum-of-templates-MC'
            self.run_closure_tests()
            self.mc += self.mc_nonclosure.sum(axis=(0,1))

        if self.num_events>0:            
            total_mc = self.mc.sum()
            print 'Scale down MC by...', self.num_events/total_mc
            self.mc *= (self.num_events/total_mc)
            self.template[:,:,:,-1] *= (self.num_events/total_mc)

        print 'Total MC in acceptance: ', self.mc.sum()

        self.save_template_snapshot(data=self.mc, title='MC', tag='mc')

        # build template with all OF bins
        self.overflow_template = np.zeros( self.mc.shape )
        for burn_qt in range(reduce_qt, 0):
            for y in range(self.bins_y_size - reduce_y):
                self.overflow_template += self.template[self.mc_mass_index, burn_qt, y, -1]
        for burn_y in range(reduce_y, 0):
            for qt in range(self.bins_qt_size):
                self.overflow_template += self.template[self.mc_mass_index, qt, burn_y, -1]                
        self.save_template_snapshot(data=self.overflow_template, title='Overflow', tag='of')


        self.dim_A = 0
        self.dim_beta = 0
        self.n_param = 0

        if use_prefit and (os.path.isfile(self.in_dir+'prefit_beta.npy') and os.path.isfile(self.in_dir+'prefit_Vbeta.npy')):
            print 'Using prefit covariance matrix/values for beta'            
            self.Vbeta_prefit = np.load(self.in_dir+'prefit_Vbeta.npy') 
            self.beta_prefit  = np.load(self.in_dir+'prefit_beta.npy') 

        self.book_parameters()

        self.out_file = ROOT.TFile(self.out_dir+'result_'+DY+'_'+charge+'_'+job_name+'.root','RECREATE')
        self.out_tree = ROOT.TTree('tree', 'tree')
        self.fit_results = {}
        for key,p in self.map_params.items():                        
            if '_true' in key:
                continue
            self.fit_results[key] = array( 'f', 5*[ 0.0 ] ) 
            self.out_tree.Branch(key, self.fit_results[key], key+'[5]/F')
        self.fit_results['minuit'] = array( 'f', 6*[ 0.0 ] ) 
        self.out_tree.Branch('minuit', self.fit_results['minuit'], 'minuit[6]/F')

        self.build_aux_matrices()

        if self.use_prior:
            self.reshape_covariance( np.load(open(self.in_dir+'covariance_'+DY+'_'+charge+'_'+var+'_stat_plus_syst_all_A0-4_forced.npy', 'r')) )
        

    def run_closure_tests(self):

        from tree_utils import get_coeff_vals
        self.mc_nonclosure = np.zeros( (self.bins_qt_size, self.bins_y_size, self.bins_eta_size, self.bins_pt_size) )

        for iy in range(self.bins_y_size):
            for iqt in range(self.bins_qt_size):
                inorm = self.template[self.mc_mass_index][iqt][iy][-1].sum()/self.mc_acceptances[self.mc_mass_index][iqt][iy]
                tMC = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][-1])
                coeff_vals = get_coeff_vals(res=self.res_coeff, 
                                            coeff_eval=('fit' if self.do_parametric else 'val'), 
                                            bin_y=self.get_y_bin(iy), qt=self.mid_point(iqt), coeff=self.coefficients,
                                            np_bins_template_qt_new=self.bins_qt)
                tUL = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][-2])
                tjA = np.zeros(tUL.shape)
                normUL = 1.0
                for icoeff_val,coeff_val in enumerate(coeff_vals[:len(self.coefficients)]):
                    tjA += copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][icoeff_val])*coeff_val
                
                residual = (normUL*tUL+tjA)*inorm-tMC
                self.mc_nonclosure[iqt,iy] += residual
                #self.save_template_snapshot(data=residual, title='MC - sum of templates', tag=self.get_y_qt_bin(iy,iqt)+'_closure')
                #self.save_template_snapshot(data=(normUL*tUL+tjA)*inorm, title='Sum of templates', tag=self.get_y_qt_bin(iy,iqt)+'_template')
                

    def reshape_covariance(self, cov_coeff):

        self.cov_coeff = np.zeros((self.dim_beta,self.dim_beta))
        print 'Reshape covariance: make matrix with shape ', self.cov_coeff.shape

        idx_A1 = 0
        for coeffA1 in self.coefficients:
            for iyA1 in range(self.bins_y_size):            
                y_binA1 = self.get_y_bin( iyA1 )        
                (valid_ordersA1, orderA1) = self.get_orders(coeffA1, y_binA1)
                for oA1 in range(orderA1+1):

                    idx_A2 = 0
                    for coeffA2 in self.coefficients:
                        for iyA2 in range(self.bins_y_size):                    
                            y_binA2 = self.get_y_bin( iyA2 )        
                            (valid_ordersA2, orderA2) = self.get_orders(coeffA2, y_binA2)
                            for oA2 in range(orderA2+1):

                                idx_B1 = 0
                                for iyB1 in range(self.bins_y_size):                    
                                    y_binB1 = self.get_y_bin( iyB1 )
                                    for coeffB1 in self.coefficients:
                                        (valid_ordersB1, orderB1) = self.get_orders(coeffB1, y_binB1)
                                        for oB1 in valid_ordersB1:
                                            idx_B2 = 0
                                            for iyB2 in range(self.bins_y_size):                    
                                                y_binB2 = self.get_y_bin( iyB2 )
                                                for coeffB2 in self.coefficients:
                                                    (valid_ordersB2, orderB2) = self.get_orders(coeffB2, y_binB2)
                                                    for oB2 in valid_ordersB2:
                                                        if iyA1==iyB1 and iyA2==iyB2 and coeffA1==coeffB1 and coeffA2==coeffB2 and oA1==oB1 and oA2==oB2:
                                                            print '\tFill entry (', idx_B1, idx_B2, ')'
                                                            self.cov_coeff[idx_B1,idx_B2] = cov_coeff[idx_A1,idx_A2]
                                                        idx_B2 += 1
                                            idx_B1 += 1

                                idx_A2 += 1
                    idx_A1 += 1

        np.set_printoptions(threshold=np.inf)
        print self.cov_coeff        
        self.Vinv_prior =  np.linalg.inv(self.cov_coeff)
        return
    
    # tell Minuit to fix a parameter
    def fix_parameter(self, par_name):

        # condition for reducing by one the n.d.o.f.:
        # <=> a parameter may be fixed in the minimisation, but still count as --dof
        subtract_dof = (self.do_parametric and 'pol' in par_name) or \
            (not self.do_parametric and ('y' in par_name and 'qt' in par_name and 'A' in par_name) )

        # check for a regexpr matching
        for fix in self.fixed_parameters:
            if fix in par_name:
                if subtract_dof:
                    self.ndof -= 1
                return True

        # check if the parameter is in the list of fixed parameters:
        if (par_name in self.fixed_parameters):
            if subtract_dof:
                self.ndof -= 1
            return True

        # reduce ndof
        self.ndof -= 1
        return False
        

    def define_parameter( self, par_name='', start=0., step=0., par_range=(), true=0. ):
        self.gMinuit.DefineParameter( self.n_param, par_name, start, step, par_range[0], par_range[1] )        
        self.map_params[par_name] = self.n_param
        if self.fix_parameter(par_name):
            self.gMinuit.FixParameter(self.n_param)
        self.map_params[par_name+'_true'] = true
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

        map_betas = range(500)

        self.define_parameter( par_name='mass', start=self.mc_mass, step=0.020, par_range=(self.mc_mass-0.500,self.mc_mass+0.500), true=self.mc_mass )
        
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for iqt in range(self.bins_qt_size):
                qt_bin = self.get_qt_bin(iqt)
                mc_norm = self.template[self.mc_mass_index][iqt][iy][-1].sum()/self.mc_acceptances[self.mc_mass_index][iqt][iy]
                self.define_parameter( par_name=y_bin+'_'+qt_bin+'_norm', 
                                       start=mc_norm, 
                                       step=math.sqrt(mc_norm)*10., 
                                       par_range=(mc_norm/1.5, mc_norm*1.5),
                                       true=mc_norm)

                for coeff in self.coefficients:
                    self.dim_A += 1
                    if not self.do_parametric:
                        map_betas[self.dim_beta] = self.n_param
                        step_beta  = math.sqrt(self.Vbeta_prefit[self.dim_beta,self.dim_beta]) if hasattr(self, 'Vbeta_prefit') else 0.01
                        start_beta = self.beta_prefit[self.dim_beta] if hasattr(self, 'beta_prefit') else self.res_coeff[coeff+'_'+y_bin+'_val'][iqt] 
                        self.define_parameter( par_name=y_bin+'_'+qt_bin+'_'+coeff, 
                                               start=start_beta,
                                               step=step_beta,
                                               par_range=(-2.,+2.),
                                               true=self.res_coeff[coeff+'_'+y_bin+'_val'][iqt])
                        self.dim_beta += 1
                        
                    elif self.do_parametric and iqt==0:
                        (valid_orders, order) = self.get_orders(coeff, y_bin)
                        for o in valid_orders:
                            map_betas[self.dim_beta] = self.n_param
                            step_beta  = math.sqrt(self.Vbeta_prefit[self.dim_beta,self.dim_beta]) if hasattr(self, 'Vbeta_prefit') else 0.1*math.pow(10,-(o+2))
                            start_beta = self.beta_prefit[self.dim_beta] if hasattr(self, 'beta_prefit') else self.res_coeff[coeff+'_'+y_bin+'_fit'][o]
                            self.define_parameter( par_name=y_bin+'_'+coeff+'_'+'pol'+str(order)+'_p'+str(o), 
                                                   start=start_beta,
                                                   step=step_beta,
                                                   par_range=(start_beta-5*step_beta, start_beta+5*step_beta),
                                                   true=self.res_coeff[coeff+'_'+y_bin+'_fit'][o])
                            self.dim_beta += 1                            

        self.map_betas = map_betas[:self.dim_beta]

        # the coefficients
        self.beta = np.zeros(self.dim_beta) 
        # approximate cov. matrix of beta
        self.Vbeta = np.identity(self.dim_beta) 


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

    def mid_point(self, iqt=-1):
         return 0.5*(self.bins_qt[iqt]+self.bins_qt[iqt+1])

    def build_aux_matrices(self):

        if self.do_parametric:
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
                                        self.K[idx_A, idx_beta] += math.pow( self.mid_point(iqt), o)
                                    self.map_betas[idx_beta] = self.map_params[y_bin2+'_'+coeff2+'_'+'pol'+str(order)+'_p'+str(o)]
                                    idx_beta += 1
                        idx_A += 1
                
        else:
            self.K = np.identity(self.dim_A)

        return 

    def fcn(self, npar, gin, f, par, iflag ):
        nll = self.chi2(par=par)
        print 'Chi2: ', '{:0.6f}'.format(nll), '/', self.ndof, ' dof = ', '{:0.4f}'.format(nll/self.ndof)
        f[0] = nll
        return

    def chi2(self, par):

        # construction of Am and Bm matrices
        (Am, bm) = (np.zeros((self.data.size, self.dim_A)), np.zeros(self.nsub.shape))
        bm = np.zeros(self.nsub.shape)
        idx_A = 0
        for iy in range(self.bins_y_size):
            for iqt in range(self.bins_qt_size):
                inorm = par[self.map_params[self.get_y_qt_bin(iy,iqt)+'_norm']]
                tUL = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][-2])
                bm += tUL.flatten()*inorm
                for icoeff,coeff in enumerate(self.coefficients):
                    tjA = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][icoeff])
                    Am[:, idx_A] += tjA.flatten()*inorm
                    idx_A += 1
        
        nsub = (self.nsub - bm)

        if self.release_for_hesse:
            beta = np.zeros( self.dim_beta )
            for ib in range(self.dim_beta):
                beta[ib] = par[ self.map_betas[ib] ]
            aux = np.linalg.multi_dot([Am, self.K, beta])
            chi2min = np.linalg.multi_dot( [(nsub-aux).T, self.Vinv, (nsub-aux)] )
            if self.use_prior:
                res2 = (beta_prior-beta)
                chi2min += np.linalg.multi_dot( [res2.T, self.Vinv_prior, res2] )
            return chi2min

        # construction of chi2
        aux1 =  np.linalg.multi_dot([self.K.T, Am.T, self.Vinv, Am, self.K])        
        if self.use_prior:
            aux1 += self.Vinv_prior
        aux1_inv = np.linalg.inv(aux1)
        aux2 = np.linalg.multi_dot( [ self.K.T, Am.T, self.Vinv, nsub ] )
        if self.use_prior:
            aux2 +=  np.linalg.multi_dot( [ self.Vinv_prior, beta_prior ])
        beta = np.linalg.multi_dot( [aux1_inv, aux2 ] )        
        res1 = (nsub - np.linalg.multi_dot( [Am, self.K, beta]) )
        chi2min = np.linalg.multi_dot( [res1.T, self.Vinv, res1 ] )
        if self.use_prior:
            res2 = (beta_prior-beta)
            chi2min += np.linalg.multi_dot( [res2.T, self.Vinv_prior, res2] )

        # save the intermediate values of beta and Vbeta
        self.beta  = beta
        self.Vbeta = aux1_inv

        return chi2min
    
    
    def load_data(self, dataset='asymov',  postfix=''):
        
        self.data = np.zeros( self.mc.shape )
        if dataset=='asymov':
            self.data += copy.deepcopy(self.mc)
        elif dataset=='random':
            self.data += np.random.poisson( self.mc )
        else:
            print 'Option not implemented'
        self.save_template_snapshot(data=self.data, title='Data', tag='data'+postfix)

        n = self.data.flatten()
        if n[n<10].size > 0:
            print n[n<10].size, ' bins with less than 10 entries!'
        self.Vinv = np.diag(1./n)
        self.nsub = n - self.overflow_template.flatten()
        self.save_template_snapshot(data=self.nsub.reshape(self.mc.shape), title='Data-overflow', tag='data_sub'+postfix)


    def run(self, n_points=500000, run_minos=False, run_post_hesse=False):

        clock = time.time()

        # Run HESSE --> MIGRAD --> MINOS
        self.arglist[0] = n_points
        self.arglist[1] = 1.0
        print("Convergence at EDM %s" % (self.arglist[1]*0.001))        

        self.release_for_hesse = False

        print 'Running HESSE...'
        self.gMinuit.mnexcm( "HES", self.arglist, 1, self.ierflg )
        self.arglist[0] = 1
        self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
        print 'Running MIGRAD with status=1...'
        status = self.gMinuit.Migrad()
        if status>0:
            self.arglist[0] = 2
            self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
            print 'Running MIGRAD with status=2...'
            status = self.gMinuit.Migrad()

        np.save(self.in_dir+'prefit_Vbeta.npy', self.Vbeta)
        np.save(self.in_dir+'prefit_beta.npy', self.beta)

        if run_minos:
            print 'Running MINOS on parameter...mass'
            (massL, massH) =  (ROOT.Double(0.), ROOT.Double(0.) )
            self.gMinuit.mnmnot(1, 1, massH, massL)

        # stop the clock
        clock -= time.time()

        (amin, edm, errdef)    = (ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.))
        (nvpar, nparx, icstat) = (ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985))
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        self.gMinuit.mnprin( 3, amin )            
        pvalue = 1-stats.chi2.cdf(float(amin), self.ndof)

        print 'Bins: ', self.data.size, ', number of d.o.f.: ', self.ndof, \
            ' chi2: '+'{:4.1f}'.format(float(amin)), ' (p-value: '+'{:4.3f}'.format(pvalue)+')'
        print('Fit done in '+'{:4.1f}'.format(-clock)+' seconds')

        self.fit_results['minuit'][0] = float(status)
        self.fit_results['minuit'][1] = float(edm)
        self.fit_results['minuit'][2] = float(amin)
        self.fit_results['minuit'][3] = float(self.ndof)
        self.fit_results['minuit'][4] = float(pvalue)
        self.fit_results['minuit'][5] = -float(clock)

        # Set coefficients to their best fit value from chi2
        for ib,b in enumerate(self.beta):
            self.arglist[0] = self.map_betas[ib]+1
            self.arglist[1] = b
            self.gMinuit.mnexcm( "SET PAR", self.arglist, 2, self.ierflg )
            if run_post_hesse:
                continue
                self.gMinuit.Release( self.map_betas[ib] )                
                #err_b_approx = math.sqrt(self.Vbeta[ib,ib])
                #self.arglist[1] = b-5*err_b_approx
                #self.arglist[2] = b+5*err_b_approx
                #self.gMinuit.mnexcm( "SET LIM", self.arglist, 3, self.ierflg )

        if run_post_hesse:
            self.release_for_hesse = True
            self.arglist[0] = n_points
            self.gMinuit.mnexcm( "HES", self.arglist, 1, self.ierflg )

        self.cov = ROOT.TMatrixDSym(self.n_param)
        #self.gMinuit.mnemat(self.cov.GetMatrixArray(), self.n_param)
        #self.plot_cov_matrix()



    def update_results(self):

        (val,err,errL,errH)  = (ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.) )

        print 'Normalisations:'
        for key,p in self.map_params.items():                        
            if not ('_norm' in key or key=='mass'):
                continue
            if '_true' in key:
                continue
            self.gMinuit.GetParameter(p, val, err)
            true = self.map_params[key+'_true']
            pull = (val-true)/err if err>0. else 0.0
            self.fit_results[key][0] = float(val)
            self.fit_results[key][1] = -float(err)
            self.fit_results[key][2] = float(err)
            self.fit_results[key][3] = float(true)
            self.fit_results[key][4] = float(pull)
            print key+' = ', '{:0.2E}'.format(val), '+/-', '{:0.2E}'.format(err), \
                ' true = ', '{:0.2E}'.format(true), \
                ' pull = ', '{:0.2f}'.format(pull)

        print 'Coefficients:'
        for key,p in self.map_params.items():                        
            if '_norm' in key or key=='mass':
                continue
            if '_true' in key:
                continue
            self.gMinuit.GetParameter(p, val, err)                
            true = self.map_params[key+'_true']
            if not self.release_for_hesse:
                idxb = -1
                for ib,b in enumerate(self.map_betas):
                    if b==p:              
                        idxb = ib
                err = ROOT.Double(math.sqrt(self.Vbeta[idxb,idxb]))
            pull = (val-true)/err if err>0. else 0.0
            self.fit_results[key][0] = float(val)
            self.fit_results[key][1] = -float(err)
            self.fit_results[key][2] = float(err)
            self.fit_results[key][3] = float(true)
            self.fit_results[key][4] = float(pull)
            print key+' = ', '{:0.2E}'.format(val), '+/-', '{:0.2E}'.format(err), \
                ' true = ', '{:0.2E}'.format(true), \
                ' pull = ', '{:0.2f}'.format(pull)
        self.out_tree.Fill()

        # reset parameters
        for key,p in self.map_params.items():                        
            if '_true' in key:
                continue
            self.arglist[0] = p+1
            self.arglist[1] = self.map_params[key+'_true']  
            self.gMinuit.mnexcm( "SET PAR", self.arglist, 2, self.ierflg )

    # save the TTree and close
    def close(self):
        self.out_tree.Scan('minuit', '', '', 1)
        self.out_file.cd()
        self.out_tree.Write("tree", ROOT.TObject.kOverwrite)
        self.out_file.Close()

    def save_template_snapshot(self, data=np.array([]), title='', tag=''):
        xx,yy = np.meshgrid(self.bins_pt, self.bins_eta)        
        plt.pcolormesh(yy, xx, data)
        if data[data<=0.].size > 0:
            plt.colorbar()
        else:
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






