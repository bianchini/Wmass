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
                 input_tag_fit='all_A0-4_forced_v4',  
                 input_tag_templ='',
                 alternative_mc='',
                 job_name='TEST', 
                 mc_mass=80.419, 
                 num_events=-1,
                 verbose=True, 
                 fixed_parameters=[],
                 prior_options=['sum'],
                 reduce_qt=-1,
                 reduce_y=-1,
                 reduce_pt=0,
                 fit_mode='parametric',
                 do_analytic_fit=False,
                 interpolation='linear',
                 use_prefit=False,
                 add_nonclosure=True,
                 use_gradient=False,
                 save_plots=['mc', 'of'],
                 random_start=False,
                 print_evals=True,
                 run_on_crab=False
                 ):

        self.in_dir  = os.environ['CMSSW_BASE']+'/src/Wmass/data/'
        self.out_dir = os.environ['CMSSW_BASE']+'/src/Wmass/test/'
        if run_on_crab:
            self.out_dir = './'
        self.job_name = job_name
        self.num_events = num_events 
        
        self.verbose = verbose
        self.setup_norm_ranges = 'hardcoded'+input_tag_templ
        self.random_start = random_start
        self.use_gradient = use_gradient

        self.print_evals = print_evals

        self.fixed_parameters = fixed_parameters
        self.interpolate_mass = ('mass' not in fixed_parameters)
        self.interpolation = interpolation

        self.prior_options = prior_options
        self.use_prior = len(prior_options)>0

        self.mc_mass = mc_mass
        self.fit_mode = fit_mode
        self.do_analytic_fit = do_analytic_fit
        self.map_params = {}
        self.map_params_minuit = {}
        self.release_for_hesse = False

        self.reduce_pt = reduce_pt
        self.reduce_qt = reduce_qt
        self.reduce_y  = reduce_y
        self.charge = charge

        self.update = True

        templates_files = {'template': 0, 'masses': 1, 'bins_qt' : 2, 
                           'bins_y'  : 3, 'coefficients' : 4, 'bins_eta': 5, 
                           'bins_pt' : 6, 'mc_acceptances' : 7 }

        self.res_coeff = np.load(open(self.in_dir+'fit_results_'+DY+'_'+charge+'_'+var+'_'+input_tag_fit+'.pkl', 'r'))        
        templates = np.load(self.in_dir+'template_'+charge+'_'+var+'_val'+input_tag_templ+'.npz')

        #########################################
        #for key,p in self.res_coeff.items():
        #    if ('A2' in key and 'fit' in key):
        #        self.res_coeff[key] += [0.,0.,0.]
        #########################################

        if alternative_mc!='':
            print 'Loading alternative MC...'
            templates_alternative = np.load(self.in_dir+'template_'+charge+'_'+var+'_val_'+alternative_mc+'.npz')
            self.res_coeff = np.load(open(self.in_dir+'fit_results_'+DY[:2]+'_'+alternative_mc+'_'+charge+'_'+var+'_'+input_tag_fit+'_alternative.pkl', 'r'))        

        # for template, save the size of (pt,eta) plane 
        # for bins, save len(bins)-1, since bins are edges
        # for coeff, save len(coeff)-2, since coeff[-1] = MC and coeff[-2] = UL
        for key,p in templates_files.items():             
            template_file = templates['arr_'+str(p)]
            if verbose:
                print ('Adding file arr_'+str(p)+' with shape...'), template_file.shape, (' to self with key... '+key)
            size = -1
            if key=='template':
                size = template_file[0,0,0,0].size
                setattr(self, key,template_file[:,:,:,:,:,0:reduce_pt] if reduce_pt<0 else template_file)
                if reduce_pt<0:
                    self.update_acceptances = np.ones( template_file.shape[0:3] )
                    self.update_acceptances *= template_file[:,:,:,-1,:,0:reduce_pt].sum(axis=(3,4))/template_file[:,:,:,-1,:,:].sum(axis=(3,4))
                if alternative_mc!='':
                    setattr(self, key+'_alt', templates_alternative['arr_'+str(p)][:,:,:,:,:,0:reduce_pt] \
                                if reduce_pt<0 else templates_alternative['arr_'+str(p)])
            elif key=='bins_qt':
                size = template_file.size - 1 + reduce_qt
                setattr(self, key, template_file[:reduce_qt] if reduce_qt<0 else template_file)
            elif key=='bins_y':
                size = template_file.size - 1 + reduce_y
                setattr(self, key, template_file[:reduce_y] if reduce_y<0 else template_file)
            elif key=='bins_pt':
                size = template_file.size - 1 + reduce_pt
                setattr(self, key,template_file[0:reduce_pt] if reduce_pt<0 else template_file)
            elif key=='bins_eta':
                size = template_file.size - 1
                setattr(self, key,template_file)
            elif key=='coefficients':
                size = template_file.size - 2
                setattr(self, key,template_file[:-2])
            else:
                size = template_file.size
                setattr(self, key,template_file)
                if alternative_mc!='':
                    setattr(self, key+'_alt', templates_alternative['arr_'+str(p)])

            setattr(self, key+'_size', size)

        print 'Fit range:'
        print '\ty  :', self.bins_y
        print '\tqt :', self.bins_qt
        print 'Observables:'
        print '\tpt :', self.bins_pt
        print '\teta:', self.bins_eta

        self.mc_mass_index = np.where(self.masses==mc_mass)[0][0]
        print 'MC mass at index...', self.mc_mass_index

        if hasattr(self, 'update_acceptances'):
            print 'Updating acceptances to account for pt burnt bins'
            self.mc_acceptances *= self.update_acceptances

        if self.interpolation=='quadratic':
            M = np.zeros((self.masses_size,self.masses_size))
            for im1,m1 in enumerate(self.masses):
                for im2,m2 in enumerate(self.masses):
                    M[im1,im2]= math.pow(m1,im2)
            self.Minv = np.linalg.inv(M)

        # remove UL from all templates
        for icoeff,coeff in enumerate(self.coefficients):
            self.template[:,:,:,icoeff] -= self.template[:,:,:,-2]

        self.mc = copy.deepcopy(self.template.sum(axis=(1,2))[self.mc_mass_index,-1])
        self.ndof = self.mc.size        
        if add_nonclosure:
            print 'Removing non closure by adding Sum-of-templates-MC'
            self.run_closure_tests(save_plots=[])
            self.mc += self.mc_nonclosure.sum(axis=(0,1))

        if self.num_events>0:            
            total_mc = self.mc.sum()
            print 'Scale MC by...', self.num_events/total_mc
            self.mc *= (self.num_events/total_mc)
            self.template[:,:,:,-1] *= (self.num_events/total_mc)

        print 'Total MC in acceptance: ', self.mc.sum()

        if 'mc' in save_plots:
            self.save_template_snapshot(data=self.mc, title='MC', tag='mc')

        # build template with all OF bins
        self.overflow_template = np.zeros( self.mc.shape )
        for burn_qt in range(reduce_qt, 0):
            for y in range(self.bins_y_size - reduce_y):
                if alternative_mc=='':
                    self.overflow_template += getattr(self, 'template')[self.mc_mass_index, burn_qt, y, -1]
                else:
                    self.overflow_template += getattr(self, 'template_alt')[0, burn_qt, y, -1]
        for burn_y in range(reduce_y, 0):
            for qt in range(self.bins_qt_size):
                if alternative_mc=='':
                    self.overflow_template += getattr(self, 'template')[self.mc_mass_index, qt, burn_y, -1]                
                else:
                    self.overflow_template += getattr(self, 'template_alt')[0, qt, burn_y, -1]                
        if 'of' in save_plots:
            self.save_template_snapshot(data=self.overflow_template, title='Overflow', tag='of')
        print 'Overflow template has integral:', self.overflow_template.sum(), ' ('+'{:4.4f}'.format(self.overflow_template.sum()/self.mc.sum())+' of total)' 

        # overflow bins in qt
        self.overflow_template_qt = np.zeros( self.mc.shape )
        self.overflow_template_qt += getattr(self, 'template')[self.mc_mass_index, reduce_qt:, :reduce_y, -1].sum(axis=(0,1))
        if 'of' in save_plots:
            self.save_template_snapshot(data=self.overflow_template_qt, title='Overflow $q_T$', tag='of_qt')

        # overflow bins in y
        self.overflow_template_y  = np.zeros( self.mc.shape )
        self.overflow_template_y += getattr(self, 'template')[self.mc_mass_index, :, reduce_y:, -1].sum(axis=(0,1))
        if 'of' in save_plots:
            self.save_template_snapshot(data=self.overflow_template_y, title='Overflow $|y|$', tag='of_y')

        print 'Cross-check OF-templates:'
        print '\tAll:', self.overflow_template.sum()
        print '\tqt:', self.overflow_template_qt.sum()
        print '\ty:', self.overflow_template_y.sum()
        print '\tDiff:', (self.overflow_template - self.overflow_template_y - self.overflow_template_qt).sum()

        self.dim_A     = 0
        self.dim_alpha = 0
        self.dim_beta  = 0
        self.n_param   = 0
        self.n_param_minuit = 0

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
            if (self.fit_mode=='parametric' or self.fit_mode=='parametric2D') and '_norm' in key:
                for coeff in self.coefficients:
                    name = key[:-5]+'_'+coeff
                    self.fit_results[name] = array( 'f', 5*[ 0.0 ] )
                    self.out_tree.Branch(name, self.fit_results[name], name+'[5]/F')

        self.fit_results['minuit'] = array( 'f', 6*[ 0.0 ] ) 
        self.out_tree.Branch('minuit', self.fit_results['minuit'], 'minuit[6]/F')
        self.fit_results['weight'] = array( 'f', [ 0.0 ] ) 
        self.out_tree.Branch('weight', self.fit_results['weight'], 'weight/F')

        self.build_aux_matrices()

        if self.use_prior:
            print 'Use prior...'+self.prior_options['prior']
            self.reshape_covariance( cov=np.load(open(self.in_dir+'covariance_'+DY+'_'+charge+'_'+var+'_'+self.prior_options['prior']+'_'+input_tag_fit+'.npy', 'r')), 
                                     dictionary=pickle.load(open(self.in_dir+'covariance_dict_'+DY+'_'+charge+'_'+var+'_'+input_tag_fit+'.pkl', 'r')), 
                                     save_plots=save_plots )
            self.ndof += self.dim_beta
        

    def reshape_covariance(self, cov=np.array([]), dictionary={}, save_plots=[]):
        V_prior = np.zeros((self.dim_beta,self.dim_beta))
        self.beta_prior =  np.zeros(self.dim_beta)

        for ib1,b1 in enumerate(self.map_betas):
            par_name1 = ''
            for key1,p1 in self.map_params.items():                        
                if '_true' in key1:
                    continue
                if p1==b1:
                    par_name1 = key1
            self.beta_prior[ib1] = self.map_params[par_name1+'_true']
            match1 = any([sel in par_name1 for sel in self.prior_options['select']]) or self.prior_options['select']==[]
            idx1 = dictionary[par_name1]
            for ib2,b2 in enumerate(self.map_betas):
                par_name2 = ''
                for key2,p2 in self.map_params.items():                        
                    if '_true' in key2:
                        continue
                    if p2==b2:
                        par_name2 = key2                    
                match2 = any([sel in par_name2 for sel in self.prior_options['select']]) or self.prior_options['select']==[]
                idx2 = dictionary[par_name2]
                V_prior[ib1,ib2] = cov[idx1,idx2]
                
                if 'ybins' in self.prior_options['decorrelate']:
                    if par_name1[0:11]!=par_name2[0:11]:
                        if self.verbose:
                            print 'Removing correlations between '+par_name1+' and '+par_name2
                        V_prior[ib1,ib2] *= 0.
                if 'coeff' in self.prior_options['decorrelate']:
                    if par_name1[12:14]!=par_name2[12:14]:
                        if self.verbose:
                            print 'Removing correlations between '+par_name1+' and '+par_name2
                        V_prior[ib1,ib2] *= 0.                    
                if 'all' in self.prior_options['decorrelate']:
                    if par_name1!=par_name2:
                        if self.verbose:
                            print 'Removing correlations between '+par_name1+' and '+par_name2
                        V_prior[ib1,ib2] *= 0.                                            

                if not match1:
                    inflate = math.pow(self.prior_options['inflate'], 2.0 if not match2 else 1.0)
                    if self.verbose:
                        print 'Inflate V_prior['+par_name1+', '+par_name2+'] by... ', inflate
                    V_prior[ib1,ib2] *= inflate

        self.V_prior = V_prior
        self.Vinv_prior = np.linalg.inv(V_prior)
        self.beta_prior_eval = copy.deepcopy(self.beta_prior)

        if 'prior' in save_plots:
            self.plot_cov_matrix(n_free=V_prior.shape[0], cov=V_prior, name='prior_'+self.prior_options['prior'])


    def run_closure_tests(self, save_plots=['mc-sum', 'sum']):

        from tree_utils import get_coeff_vals
        self.mc_nonclosure = np.zeros( (self.bins_qt_size, self.bins_y_size, self.bins_eta_size, self.bins_pt_size) )

        for iy in range(self.bins_y_size):
            for iqt in range(self.bins_qt_size):
                inorm = self.template[self.mc_mass_index][iqt][iy][-1].sum()/self.mc_acceptances[self.mc_mass_index][iqt][iy]
                tMC = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][-1])
                coeff_vals = get_coeff_vals(res=self.res_coeff, 
                                            #coeff_eval=('fit' if self.fit_mode=='parametric' else 'val'), 
                                            coeff_eval='val', 
                                            bin_y=self.get_y_bin(iy), qt=self.mid_point_qt(iqt), coeff=self.coefficients,
                                            np_bins_template_qt_new=self.bins_qt)
                tUL = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][-2])
                tjA = np.zeros(tUL.shape)
                normUL = 1.0
                for icoeff_val,coeff_val in enumerate(coeff_vals[:len(self.coefficients)]):
                    tjA += copy.deepcopy(self.template[self.mc_mass_index][iqt][iy][icoeff_val])*coeff_val
                
                residual = (normUL*tUL+tjA)*inorm-tMC
                #print self.get_y_bin(iy), self.get_qt_bin(iqt), residual.sum()/math.sqrt(tMC.sum())
                self.mc_nonclosure[iqt,iy] += residual
                if 'mc-sum' in save_plots:
                    self.save_template_snapshot(data=residual, title='MC - sum of templates', tag=self.get_y_qt_bin(iy,iqt)+'_closure')
                if 'sum' in save_plots:
                    self.save_template_snapshot(data=(normUL*tUL+tjA)*inorm, title='Sum of templates', tag=self.get_y_qt_bin(iy,iqt)+'_template')
                
    # Create a random dataset with smeared values of the coefficents
    def randomize_coefficients(self, randomize_in_generation=False):
        from tree_utils import get_coeff_vals

        # Random instance of betas sampled from prior covariance matrix
        beta_rnd = np.random.multivariate_normal(self.beta_prior, self.V_prior)
        if self.verbose:
            print 'Old betas:', self.beta_prior
            print 'New betas:', beta_rnd
        # Use 'smeared' value of beta_prior in chi2 function
        self.beta_prior_eval = copy.deepcopy(beta_rnd)
        if not randomize_in_generation:
            return

        print 'Randomizing coefficients in toy generation...'
        self.mc_random_coeff = np.zeros( (self.bins_qt_size, self.bins_y_size, self.bins_eta_size, self.bins_pt_size) )
        for iy1 in range(self.bins_y_size):
            y_bin1 = self.get_y_bin(iy1)
            for iqt in range(self.bins_qt_size):
                tUL = copy.deepcopy(self.template[self.mc_mass_index][iqt][iy1][-2])
                tjA = np.zeros(tUL.shape)
                normUL = 1.0
                inorm = self.template[self.mc_mass_index][iqt][iy1][-1].sum()/self.mc_acceptances[self.mc_mass_index][iqt][iy1]
                for icoeff1,coeff1 in enumerate(self.coefficients): 
                    (valid_orders1, order1) = self.get_orders(coeff1, y_bin1)                    
                    coeff_val = 0.0
                    for io1,o1 in enumerate(valid_orders1):
                        counter_beta2 = 0
                        for iy2 in range(self.bins_y_size):
                            y_bin2 = self.get_y_bin(iy2)
                            for icoeff2,coeff2 in enumerate(self.coefficients):
                                (valid_orders2, order2) = self.get_orders(coeff2, y_bin2)
                                for io2,o2 in enumerate(valid_orders2):                                    
                                    if iy1==iy2 and icoeff1==icoeff2 and io1==io2:
                                       counter_beta1 = counter_beta2 
                                    counter_beta2 += 1
                        coeff_val += math.pow( self.mid_point_qt(iqt), o1)*beta_rnd[counter_beta1]
                        #coeff_val += math.pow( self.mid_point_qt(iqt), o1)*self.beta_prior_backup[counter_beta1]
                    tjA += copy.deepcopy(self.template[self.mc_mass_index][iqt][iy1][icoeff1])*coeff_val

                residual = (normUL*tUL+tjA)*inorm
                self.mc_random_coeff[iqt,iy1] += residual

    # tell Minuit to fix a parameter
    def fix_parameter(self, par_name):
        # condition for reducing by one the n.d.o.f.:
        # <=> a parameter may be fixed in the minimisation, but still count as --dof
        subtract_dof = (self.fit_mode=='parametric' and 'pol' in par_name) or \
            (self.fit_mode=='binned' and ('y' in par_name and 'qt' in par_name and 'A' in par_name) )
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
        

    def define_parameter( self, par_name='', start=0., step=0., par_range=(), true=0., par_type='' ):
        self.gMinuit.DefineParameter( self.n_param, par_name, start, step, par_range[0], par_range[1] )        
        self.map_params[par_name] = self.n_param
        
        if par_type=='beta':
            self.map_betas[self.dim_beta] = self.n_param
            self.dim_beta += 1
        elif par_type=='alpha':
             self.map_alphas[self.dim_alpha] = self.n_param
             self.dim_alpha += 1

        if self.fix_parameter(par_name):
            self.gMinuit.FixParameter(self.n_param)
        else:
            self.map_params_minuit[par_name] = self.n_param_minuit
            self.n_param_minuit += 1

        self.map_params[par_name+'_true'] = true
        self.n_param += 1

    # return (scale, scale_range_low, scale_range_high)
    def norm_ranges(self, par_name='', bin=()):
        out = (0.1, 1., 1.)

        if self.verbose:
            print self.setup_norm_ranges
        if self.setup_norm_ranges=='hardcoded':
            rel_err = max(self.mid_point_qt(bin[1])*0.01, 0.1)* \
                (1.0 + math.pow(self.mid_point_y(bin[0])/2.5, 2.0) )
            rel_err /= math.sqrt(self.num_events/1.5e+06 if self.num_events>0 else 1.0)
            n_sigmas = 5
            out = (rel_err, (1-rel_err*n_sigmas), (1+rel_err*n_sigmas))
                
        elif self.setup_norm_ranges=='hardcoded_finer_y':
            rel_err = max(self.mid_point_qt(bin[1])*0.02, 0.2)* \
                (1.0 + math.pow(self.mid_point_y(bin[0])/2.5, 2.0) )            
            rel_err /= math.sqrt(self.num_events/1.5e+06 if self.num_events>0 else 1.0)
            n_sigmas = 5
            out = (rel_err, (1-rel_err*n_sigmas), (1+rel_err*n_sigmas))            

        elif 'finer_y_qt32' in self.setup_norm_ranges:
            rel_err = max(self.mid_point_qt(bin[1])*0.03, 0.2)* \
                (1.0 + math.pow(self.mid_point_y(bin[0])/2.5, 2.0) )
            rel_err /= math.sqrt(self.num_events/1.5e+06 if self.num_events>0 else 1.0)
            n_sigmas = 5
            out = (rel_err, (1-rel_err*n_sigmas), (1+rel_err*n_sigmas))            

        if self.verbose:
            print par_name, '=> set step at ', '{:0.2f}'.format(out[0])+'*NORM and range at ['+ \
                '{:0.2f}'.format(out[1])+', '+'{:0.2f}'.format(out[2])+']*NORM'
        return out


    def book_parameters(self):

        self.gMinuit = ROOT.TMinuit(500)  
        self.gMinuit.SetPrintLevel(1)
        if not self.verbose:
            self.gMinuit.SetPrintLevel(-1)
        self.gMinuit.SetFCN( self.fcn ) 
    
        self.arglist = array( 'd', 10*[0.] )
        self.ierflg = ROOT.Long(0)
        
        # 0.5: chi2, 1.0: -nll
        self.arglist[0] = 1.0
        self.gMinuit.mnexcm( "SET ERR", self.arglist, 1, self.ierflg )

        self.map_alphas = range(1000)
        self.map_betas  = range(1000)

        # mass
        self.define_parameter( par_name='mass', start=self.mc_mass, step=0.100, 
                               par_range=(self.mc_mass-1.0,self.mc_mass+1.0), 
                               true=self.mc_mass, 
                               par_type='alpha' )

        # qt/y bins and coefficients
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for iqt in range(self.bins_qt_size):
                qt_bin = self.get_qt_bin(iqt)
                mc_norm = self.template[self.mc_mass_index][iqt][iy][-1].sum()/ \
                    self.mc_acceptances[self.mc_mass_index][iqt][iy]

                if hasattr(self, 'template_alt'):
                    mc_norm = self.template_alt[0][iqt][iy][-1].sum()/ \
                        self.mc_acceptances_alt[0][iqt][iy]
                    
                norm_range = self.norm_ranges(par_name=y_bin+'_'+qt_bin+'_norm', bin=(iy,iqt))
                self.define_parameter( par_name=y_bin+'_'+qt_bin+'_norm', 
                                       start=( mc_norm + np.random.normal(mc_norm, math.sqrt(mc_norm)) if self.random_start else mc_norm),
                                       step=mc_norm*norm_range[0], 
                                       par_range=(mc_norm*norm_range[1], mc_norm*norm_range[2]),
                                       true=mc_norm,
                                       par_type=('alpha' if not self.do_analytic_fit else 'beta'))

                if self.do_analytic_fit:
                    self.dim_A += 1

                for coeff in self.coefficients:
                    self.dim_A += 1

                    if self.fit_mode=='binned':
                        step_beta  = math.sqrt(self.Vbeta_prefit[self.dim_beta,self.dim_beta]) if hasattr(self, 'Vbeta_prefit') else 0.01
                        start_beta = self.beta_prefit[self.dim_beta] if hasattr(self, 'beta_prefit') else self.res_coeff[coeff+'_'+y_bin+'_val'][iqt] 
                        self.define_parameter( par_name=y_bin+'_'+qt_bin+'_'+coeff, 
                                               start=start_beta,
                                               step=step_beta,
                                               par_range=(-2.,+2.),
                                               true=self.res_coeff[coeff+'_'+y_bin+'_val'][iqt],
                                               par_type='beta')
                        
                    elif self.fit_mode=='parametric':
                        if iqt>0:
                            continue
                        (valid_orders, order) = self.get_orders(coeff, y_bin)
                        for o in valid_orders:
                            step_beta  = math.sqrt(self.Vbeta_prefit[self.dim_beta,self.dim_beta]) if hasattr(self, 'Vbeta_prefit') else 0.1*math.pow(10,-(o+2))
                            start_beta = self.beta_prefit[self.dim_beta] if hasattr(self, 'beta_prefit') else self.res_coeff[coeff+'_'+y_bin+'_fit'][o]
                            self.define_parameter( par_name=y_bin+'_'+coeff+'_'+'pol'+str(order)+'_p'+str(o), 
                                                   start=start_beta,
                                                   step=step_beta,
                                                   par_range=(start_beta-5*step_beta, start_beta+5*step_beta),
                                                   true=self.res_coeff[coeff+'_'+y_bin+'_fit'][o],
                                                   par_type='beta')

                    elif self.fit_mode=='parametric2D':
                        if iqt>0 or iy>0:
                            continue
                        (valid_orders_y, order_y)   = self.get_orders(coeff, '')
                        (valid_orders_qt, order_qt) = self.get_orders(coeff, y_bin)
                        for oy in valid_orders_y:
                            for oqt in valid_orders_qt:
                                step_beta  = 0.1
                                start_beta = 0.0
                                self.define_parameter( par_name=coeff+'_'+'pol'+str(order_y)+'_pol'+str(order_qt)+'_p'+str(oy)+'_p'+str(oqt), 
                                                       start=start_beta,
                                                       step=step_beta,
                                                       par_range=(start_beta-5*step_beta, start_beta+5*step_beta),
                                                       true=0.,
                                                       par_type='beta')

        self.map_alphas = self.map_alphas[:self.dim_alpha]
        self.map_betas  = self.map_betas[:self.dim_beta]

        # the coefficients
        self.beta = np.zeros(self.dim_beta) 
        # conditional cov. matrix of beta
        self.Vbeta = np.identity(self.dim_beta) 


    def get_orders(self, coeff='', y_bin=''):               
        valid_orders = []

        if y_bin=='':
            if coeff in ['A0','A2','A3']:
                return ([0,2,4], 4)
            else:
                return ([1,3,5], 5)

        order = len(self.res_coeff[coeff+'_'+y_bin+'_fit'])-1
        for io,o in enumerate(self.res_coeff[coeff+'_'+y_bin+'_fit']):
            if o!=0.:
                valid_orders.append(io)
        
        ###############################
        #if coeff in ['A2']:
        #    return ([0,1,2,3,4], 5)
            #return ([], 0)
        #if coeff in ['A1', 'A3']:
        #    return ([1,2,3], 3)
            #return ([], 0)
        #if coeff in ['A4']:
        #    return ([0], 0)                            
        #    return ([0,1], 1)   
        ###############################

        return (valid_orders, order)                

    def get_y_bin(self, iy=-1):
        return 'y{:03.2f}'.format(self.bins_y[iy])+'_'+'y{:03.2f}'.format(self.bins_y[iy+1])

    def get_qt_bin(self, iqt=-1):
        return 'qt{:03.1f}'.format(self.bins_qt[iqt])+'_'+'qt{:03.1f}'.format(self.bins_qt[iqt+1])

    def get_y_qt_bin(self, iy=-1, iqt=-1):
        return self.get_y_bin(iy)+'_'+self.get_qt_bin(iqt)

    def mid_point_qt(self, iqt=-1):
         return 0.5*(self.bins_qt[iqt]+self.bins_qt[iqt+1])

    def mid_point_y(self, iy=-1):
         return 0.5*(self.bins_y[iy]+self.bins_y[iy+1])

    def width_qt(self, iqt=-1):
         return -(self.bins_qt[iqt]-self.bins_qt[iqt+1])

    def width_y(self, iy=-1):
         return -(self.bins_y[iy]-self.bins_y[iy+1])


    def build_aux_matrices(self):

        if self.fit_mode=='parametric' and not self.do_analytic_fit:
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
                                        self.K[idx_A, idx_beta] += math.pow( self.mid_point_qt(iqt), o)
                                    idx_beta += 1
                        idx_A += 1

   
        elif self.fit_mode=='parametric2D':
            self.K = np.zeros( ( self.dim_A, self.dim_beta) )
            # first loop
            idx_A = 0
            for iy in range(self.bins_y_size):
                for iqt in range(self.bins_qt_size):
                    for coeff1 in self.coefficients:                        
                        # second loop
                        idx_beta = 0
                        for coeff2 in self.coefficients:
                            (valid_orders_y, order_y)   = self.get_orders(coeff2, '')
                            (valid_orders_qt, order_qt) = self.get_orders(coeff2, self.get_y_bin( 0 ))
                            for oy in valid_orders_y:
                                for oqt in valid_orders_qt:
                                    if coeff1==coeff2:
                                        self.K[idx_A, idx_beta] += math.pow( self.mid_point_y(iy), oy)*math.pow( self.mid_point_qt(iqt), oqt)
                                    idx_beta += 1
                        idx_A += 1            

        elif self.fit_mode=='binned':
            self.K = np.identity(self.dim_A)
        else:
            print self.fit_mode+' not supported'

        return 


    def interpolate(self, mass, iqt, iy, icoeff):
        if not self.interpolate_mass:
            return (self.template[self.mc_mass_index, iqt, iy, icoeff])
        if self.interpolation=='linear':
            (im_low,im_high)  = (np.where(self.masses<=mass)[0][-1],  np.where(self.masses>mass)[0][0])
            r = (mass-self.masses[im_low])/(self.masses[im_high]-self.masses[im_low])
            return (1-r)*self.template[im_low, iqt, iy, icoeff] + r*self.template[im_high, iqt, iy, icoeff]
        elif self.interpolation=='quadratic':
            theta = np.einsum('ij,jkl->ikl', self.Minv, self.template[:, iqt, iy, icoeff, :, :])
            return np.einsum('ikl,i', theta, np.array([1.0, mass, mass*mass]))


    def fcn(self, npar, gin, f, par, iflag ):
        nll = self.chi2(par=par)
        f[0] = nll
        if self.use_gradient:
            gin[0] = self.grad(par=par) 
        return

    def grad(self, par):
        der = [0.]*self.n_param
        return der
 
    def chi2(self, par):

        ##############################
        #return math.pow((par[0]-self.mc_mass)/0.1,2.0)
        ##############################

        # construction of Am and Bm matrices
        (Am, bm) = (np.zeros((self.data.size, self.dim_A)), np.zeros(self.nsub.shape))
        bm = np.zeros(self.nsub.shape)
 
        if self.do_analytic_fit:
            idx_A = 0
            for iy in range(self.bins_y_size):
                for iqt in range(self.bins_qt_size):
                    tUL = self.interpolate(mass=par[0], iqt=iqt, iy=iy, icoeff=-2)
                    Am[:, idx_A] += tUL.flatten()        
                    idx_A += 1
                    for icoeff,coeff in enumerate(self.coefficients):
                        tjA = self.interpolate(mass=par[0], iqt=iqt, iy=iy, icoeff=icoeff) 
                        Am[:, idx_A] += tjA.flatten()
                        idx_A += 1
            nsub = self.nsub

        else:
            idx_A = 0
            for iy in range(self.bins_y_size):
                for iqt in range(self.bins_qt_size):
                    inorm = par[self.map_params[self.get_y_qt_bin(iy,iqt)+'_norm']]
                    tUL = self.interpolate(mass=par[0], iqt=iqt, iy=iy, icoeff=-2) 
                    bm += tUL.flatten()*inorm
                    for icoeff,coeff in enumerate(self.coefficients):
                        tjA = self.interpolate(mass=par[0], iqt=iqt, iy=iy, icoeff=icoeff) 
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
                res2 = (self.beta_prior_eval-beta)
                chi2_prior = np.linalg.multi_dot( [res2.T, self.Vinv_prior, res2] )
                chi2min += chi2_prior
            if self.print_evals:
                print 'Chi2: ', '{:0.5f}'.format(chi2min), \
                    (' = '+'{:0.5f}'.format(chi2min-chi2_prior)+' (stat.) + '+'{:0.5f}'.format(chi2_prior)+' (prior)' if self.use_prior else ''), \
                    '/', self.ndof, ' dof = ', '{:0.5f}'.format(chi2min/self.ndof)

            #################################
            # Debug in case of negative norms
            #if self.verbose:
            #    for key,p in self.map_params.items():
            #        if 'true' in key or '_norm' not in key:
            #            continue
            #        if par[p]<0.:
            #            print '\t', key, par[p]
            #################################
            self.mu = -(nsub-aux)+self.n
            return chi2min

        # construction of chi2
        aux1 =  np.linalg.multi_dot([self.K.T, Am.T, self.Vinv, Am, self.K])        
        if self.use_prior:
            aux1 += self.Vinv_prior
        aux1_inv = np.linalg.inv(aux1)
        aux2 = np.linalg.multi_dot( [ self.K.T, Am.T, self.Vinv, nsub ] )
        if self.use_prior:
            aux2 +=  np.linalg.multi_dot( [ self.Vinv_prior, self.beta_prior_eval ])
        beta = np.linalg.multi_dot( [aux1_inv, aux2 ] )        
        res1 = (nsub - np.linalg.multi_dot( [Am, self.K, beta]) )
        chi2min = np.linalg.multi_dot( [res1.T, self.Vinv, res1 ] )

        if self.use_prior:
            res2 = (self.beta_prior_eval-beta)
            chi2_prior = np.linalg.multi_dot( [res2.T, self.Vinv_prior, res2] )
            chi2min += chi2_prior

        if self.print_evals:
            print 'Chi2: ', '{:0.5f}'.format(chi2min), \
                (' = '+'{:0.5f}'.format(chi2min-chi2_prior)+' (stat.) + '+'{:0.5f}'.format(chi2_prior)+' (prior)' if self.use_prior else ''), \
                '/', self.ndof, ' dof = ', '{:0.5f}'.format(chi2min/self.ndof)
            
        # save the intermediate values of beta and Vbeta
        self.beta  = beta
        self.Vbeta = aux1_inv
        self.mu = -res1+self.n

        #################################
        # debug in case of negative norms
        #if self.verbose:
        #    for key,p in self.map_params.items():
        #        if 'true' in key or '_norm' not in key:
        #            continue
        #        if par[p]<0.:
        #            print '\t', key, par[p]
        #################################

        return chi2min
    
    
    def load_data(self, dataset='asimov',  save_plots=['data', 'data-overflow'], postfix='', scale_id=0):
        
        self.dataset_type=dataset
        self.data = np.zeros( self.mc.shape )

        if dataset=='asimov':
            self.data += copy.deepcopy( self.mc )
            print 'Loading asimov dataset as DATA with', self.data.sum(), 'entries'
 
        elif ('asimov-scaled' in dataset or 'asimov-shifted' in dataset):
            self.fit_results['weight'][0] = scale_id 
            acc = dataset.split('-')[2]            
            sub_acc = '_'+dataset.split('-')[3] if acc=='out' and len(dataset.split('-'))>3 and dataset.split('-')[3] in ['y','qt'] else ''

            if 'scaled' in dataset:
                input_tag = 'mc_weights'
            elif 'shifted' in dataset:
                input_tag = 'pt_scales'

            template_scaled_file  = np.load(self.in_dir+'template_'+self.charge+'_'+input_tag+'.npz')
            template_scaled_full = template_scaled_file['arr_0']
            iweight = np.where(template_scaled_file['arr_4']==scale_id)[0][0]
            iweight_nominal = np.where(template_scaled_file['arr_4']==0)[0][0]
            print 'Nominal template at place', iweight_nominal
            
            imass = np.where( template_scaled_file['arr_1']==self.mc_mass)[0][0]
            if acc!='full':
                (all_qt, burnt_qt, all_y, burnt_y)   = (template_scaled_file['arr_2'].size-1 ,                                                             
                                                        template_scaled_file['arr_2'].size-np.where( template_scaled_file['arr_2']==self.bins_qt[-1])[0][0]-1,
                                                        template_scaled_file['arr_3'].size-1 , 
                                                        template_scaled_file['arr_3'].size-np.where( template_scaled_file['arr_3']==self.bins_y[-1])[0][0]-1 
                                                        )
                print '(imass, all_qt, burnt_qt, all_y, burnt_y) =', (imass, all_qt, burnt_qt, all_y, burnt_y) 

            template_scaled = np.zeros(template_scaled_full.shape[3:]) 
            if acc=='full':
                template_scaled += template_scaled_full[imass,:,:,:,:,:].sum(axis=(0,1))
            elif acc=='in':
                template_scaled += template_scaled_full[imass,:-burnt_qt,:-burnt_y:,:,:].sum(axis=(0,1))
            elif acc=='out':
                if sub_acc in ['', '_qt']:
                    template_scaled += template_scaled_full[imass,-burnt_qt:,:-burnt_y,:,:,:].sum(axis=(0,1))
                if sub_acc in ['', '_y']:
                    template_scaled += template_scaled_full[imass,:,-burnt_y:,:,:,:].sum(axis=(0,1))

            mc_nominal = template_scaled[iweight_nominal,       :,0:self.reduce_pt] if self.reduce_pt<0 else template_scaled[iweight_nominal,       :,:]
            mc_scaled  = template_scaled[iweight,:,0:self.reduce_pt] if self.reduce_pt<0 else template_scaled[iweight,:,:]
            mc_nominal += 1e+01
            mc_scaled  += 1e+01
            scaling = mc_scaled/mc_nominal
            if acc=='full':
                mc_distorted = self.mc*scaling
                tmp_data_norm = 0.0
            else:
                mc_distorted = (self.mc-getattr(self, 'overflow_template'+sub_acc))*scaling if acc=='in' else getattr(self, 'overflow_template'+sub_acc)*scaling
                self.data += getattr(self, 'overflow_template'+sub_acc) if acc=='in' else (self.mc-getattr(self, 'overflow_template'+sub_acc))
                tmp_data_norm = self.data.sum()
            self.data += mc_distorted/mc_distorted.sum()*(self.num_events-tmp_data_norm)
            if 'pulls' in save_plots:
                self.save_template_snapshot(data=(self.data-self.mc)/np.sqrt(self.mc), title='Prefit pulls', tag='prefit_pulls_'+acc+sub_acc+'_weight'+str(scale_id))            
                self.save_template_snapshot(data=scaling, title='Scaling', tag='scaling_'+acc+sub_acc+'_weight'+str(scale_id))            
            print 'Loading asimov dataset scaled by weight '+str(scale_id)+' as DATA with', self.data.sum(), 'entries'

        elif dataset=='random':
            if self.use_prior:
                self.randomize_coefficients( randomize_in_generation=False )
            self.data += np.random.poisson( self.mc )
            print 'Loading poisson dataset as DATA with', self.data.sum(), 'entries'

        elif dataset=='random-coefficients':
            self.randomize_coefficients( randomize_in_generation=True )
            self.data += np.random.poisson(self.overflow_template + self.mc_random_coeff.sum(axis=(0,1)))
            print 'Loading dataset with coeff ~ V_prior(coeff) as DATA with', self.data.sum(), 'entries'

        elif dataset=='normal':
            data_tmp = np.zeros( self.mc.shape ) 
            while (data_tmp[data_tmp<=0.].size > 0):
                data_tmp = np.random.normal( self.mc, np.sqrt(self.mc) )
            self.data += data_tmp
            print 'Loading normal dataset as DATA with', self.data.sum(), 'entries'

        elif dataset=='alternative':
            self.data += copy.deepcopy(self.template_alt.sum(axis=(1,2))[0,-1])
            print 'Loading alternative dataset as DATA with', self.data.sum(), 'entries'

        else:
            print 'Option not implemented'

        if 'data' in save_plots:
            self.save_template_snapshot(data=self.data, title='Data', tag='data'+postfix)

        n = self.data.flatten()
        if n[n<10].size > 0:
            print n[n<10].size, ' bins with less than 10 entries!'
        self.Vinv = np.diag(1./n)
        self.nsub = n - self.overflow_template.flatten()
        self.n = copy.deepcopy(n)
        if 'data-overflow' in save_plots:
            self.save_template_snapshot(data=self.nsub.reshape(self.mc.shape), title='Data-overflow', tag='data_sub'+postfix)
        self.mu = np.zeros(self.nsub.shape)

        return


    def run(self, n_points=500000, strategy=2, tolerance=0.1, run_minos=False, run_post_hesse=False):

        clock = time.time()

        if self.do_analytic_fit:
            chi2min =  self.chi2(par=[self.mc_mass])
            status = 1
            # stop the clock
            clock -= time.time()
            pvalue = 1-stats.chi2.cdf(chi2min, self.ndof)            
            print 'Bins: ', self.data.size, ', number of d.o.f.: ', self.ndof, \
                ' chi2: '+'{:4.1f}'.format(chi2min), ' (p-value: '+'{:4.3f}'.format(pvalue)+')'
            print('Fit done in '+'{:4.1f}'.format(-clock)+' seconds')
            idx_beta = 0
            for iy in range(self.bins_y_size):
                y_bin = self.get_y_bin(iy)
                for iqt in range(self.bins_qt_size):                    
                    qt_bin = self.get_qt_bin(iqt)
                    val = self.beta[idx_beta]
                    true = self.map_params[y_bin+'_'+qt_bin+'_norm_true']
                    err = math.sqrt(self.Vbeta[idx_beta,idx_beta])
                    pull = (val-true)/err
                    idx_beta += 1
                    print y_bin+'_'+qt_bin+'_norm'+' = ', '{:0.2E}'.format(val), '+/-', '{:0.2E}'.format(err), \
                        ' true = ', '{:0.2E}'.format(true), \
                        ' pull = ', '{:0.2f}'.format(pull)
                    for icoeff,coeff in enumerate(self.coefficients):
                        idx_beta += 1        
            self.plot_cov_matrix(n_free=self.Vbeta.shape[0], cov=self.Vbeta, name='analytic')
            self.mu_postfit = copy.deepcopy(self.mu)
            pulls  = (self.data-self.mu_postfit.reshape(self.data.shape))/np.sqrt(self.data)
            print 'pulls sum = '+'{:3.3f}'.format(pulls.sum())+', '+'pulls2 sum = '+'{:3.3f}'.format(np.power(pulls,2).sum())
            self.save_template_snapshot(data=pulls, title='Postfit pulls', tag='postfit_pulls')
            self.save_template_snapshot(data=np.power(pulls,2), title='Postfit pulls squared', tag='postfit_pulls2')
            return status

        # this flag freezes beta's in chi2
        self.release_for_hesse = False

        self.arglist[0] = 0
        if self.use_gradient:
            self.gMinuit.mnexcm( "SET GRA", self.arglist, 1, self.ierflg)
        else:
            self.gMinuit.mnexcm( "SET NOG", self.arglist, 0, self.ierflg)

        # Run HESSE --> MIGRAD --> MINOS
        print 'Running HESSE...'
        self.arglist[0] = n_points
        self.gMinuit.mnexcm( "HES", self.arglist, 1, self.ierflg )

        print 'Running MIGRAD with strategy', strategy, '...'
        self.arglist[0] = strategy
        self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
        #status = self.gMinuit.Migrad()
        self.arglist[0] = n_points
        self.arglist[1] = tolerance
        print("Convergence at EDM %s" % (self.arglist[1]*0.001))         
        self.gMinuit.mnexcm( "MIG", self.arglist, 2, self.ierflg )
        status = copy.deepcopy(int(self.ierflg))

        if status>0:
            return 1
            self.arglist[0] = 2
            self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
            print 'Running MIGRAD with strategy', int(self.arglist[0]), '...'
            status = self.gMinuit.Migrad()

        self.Vbeta_min = copy.deepcopy(self.Vbeta)
        self.beta_min = copy.deepcopy(self.beta)
        self.mu_postfit = copy.deepcopy(self.mu)

        np.save(self.in_dir+'prefit_Vbeta.npy', self.Vbeta_min)
        np.save(self.in_dir+'prefit_beta.npy', self.beta_min)

        if run_minos:
            print 'Running MINOS on parameter...mass'
            (self.massL, self.massH) =  (ROOT.Double(0.), ROOT.Double(0.) )
            self.gMinuit.mnmnot(1, 1, self.massH, self.massL)

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
        for ib,b in enumerate(self.beta_min):
            self.arglist[0] = self.map_betas[ib]+1
            self.arglist[1] = b
            self.gMinuit.mnexcm( "SET PAR", self.arglist, 2, self.ierflg )
            if run_post_hesse:
                continue
                self.gMinuit.Release( self.map_betas[ib] )                

        if run_post_hesse:
            self.release_for_hesse = True
            self.arglist[0] = n_points
            self.gMinuit.mnexcm( "HES", self.arglist, 1, self.ierflg )

        self.TMatrix = ROOT.TMatrixDSym(self.n_param_minuit)
        self.gMinuit.mnemat(self.TMatrix.GetMatrixArray(), self.n_param_minuit)

        return 0

    def scan(self, par_name='', step=1, par_low=-1, par_high=+1):
        self.print_evals = False
        par = [0.]*self.n_param 
        (val,err)  = (ROOT.Double(0.), ROOT.Double(0.))
        nsteps = int((par_high-par_low)/step)
        print par_name
        for s in range(nsteps):
            param_val = par_low + step*s
            for key,p in self.map_params.items(): 
                if '_true' in key:
                    continue
                self.gMinuit.GetParameter(p, val, err)            
                par[p] = float(val) if key!=par_name else param_val
            chi2 = self.chi2(par=par)
            print param_val, ' ===> ',  '{:0.2f}'.format(chi2) 


    def update_results(self, print_results=True, save_plots=['coeff', 'polynom', 'norm', 'cov', 'postfit'], propagate_covariance=False):

        if 'cov' in save_plots:
            self.plot_cov_matrix(n_free=self.gMinuit.GetNumFreePars(), cov=self.TMatrix, name='fit')

        (val,err,errL,errH)  = (ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.) )

        if propagate_covariance:
            self.propagate_covariance()

        self.gMinuit.GetParameter(self.map_params['mass'], val, err)
        true = self.map_params['mass'+'_true']
        pull = (val-true)/err if err>0. else 0.0
        self.fit_results['mass'][0] = float(val)
        self.fit_results['mass'][1] = -float(err) if not hasattr(self, 'massL') else self.massL-val
        self.fit_results['mass'][2] = float(err)  if not hasattr(self, 'massH') else self.massH-val
        self.fit_results['mass'][3] = float(true)
        self.fit_results['mass'][4] = float(pull)
        if print_results:
            print 'Mass:'
            print 'mass'+' = ', '{:0.3f}'.format(val), '+/-', '{:0.3f}'.format(err), \
                ' true = ', '{:0.3f}'.format(true), \
                ' pull = ', '{:0.2f}'.format(pull), (' ' if not hasattr(self, 'massL') else ' MINOS: ['+'{:0.3f}'.format(self.massL-val)+', +'+'{:0.3f}'.format(self.massH-val)+']')

        if print_results:
            print 'Normalisations:'
        for key,p in self.map_params.items():                        
            if not '_norm' in key :
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
            if print_results:
                print key+' = ', '{:0.2E}'.format(val), '+/-', '{:0.2E}'.format(err), \
                    ' true = ', '{:0.2E}'.format(true), \
                    ' pull = ', '{:0.2f}'.format(pull)

        if print_results:
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
                err = ROOT.Double(math.sqrt(self.Vbeta_min[idxb,idxb]))
            pull = (val-true)/err if err>0. else 0.0
            self.fit_results[key][0] = float(val)
            self.fit_results[key][1] = -float(err)
            self.fit_results[key][2] = float(err)
            self.fit_results[key][3] = float(true)
            self.fit_results[key][4] = float(pull)
            if print_results:
                print key+' = ', '{:0.2E}'.format(val), '+/-', '{:0.2E}'.format(err), \
                    ' true = ', '{:0.2E}'.format(true), \
                    ' pull = ', '{:0.2f}'.format(pull)


        print 'Propagate errors to Ai'
        rnd_As = np.zeros( (self.dim_A, 500) )
        for itoy in range(500):
            rnd_betas = np.random.multivariate_normal(self.beta_min, self.Vbeta_min) 
            rnd_A = np.linalg.multi_dot([self.K, rnd_betas])
            rnd_As[:, itoy] = rnd_A
        (cov_A, mean_A) = (np.cov(rnd_As), rnd_As.mean(axis=1))
        idx_A = 0
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for iqt in range(self.bins_qt_size):
                qt_bin = self.get_qt_bin(iqt)                               
                for icoeff,coeff in enumerate(self.coefficients):
                    self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][0] = mean_A[idx_A] 
                    self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][1] = -math.sqrt(cov_A[idx_A,idx_A])
                    self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][2] = +math.sqrt(cov_A[idx_A,idx_A])
                    self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][3] = self.res_coeff[coeff+'_'+y_bin+'_val'][iqt]
                    self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][4] = \
                        (mean_A[idx_A]-self.res_coeff[coeff+'_'+y_bin+'_val'][iqt])/math.sqrt(cov_A[idx_A,idx_A])
                    if print_results:
                        print y_bin+'_'+qt_bin+'_'+coeff+' = ', '{:0.2E}'.format(self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][0]), '+/-', \
                            '{:0.2E}'.format(self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][2]), \
                            ' true = ', '{:0.2E}'.format(self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][3]), \
                            ' pull = ', '{:0.2f}'.format(self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][4] )
                    idx_A += 1


        if 'polynom' in save_plots:
            if self.fit_mode=='parametric':
                self.plot_results_polynom_y(var='resolution')        

        if 'coeff' in save_plots:
            if self.fit_mode=='parametric2D' or self.fit_mode=='parametric':
                self.plot_results_coeff_y_qt(var='resolution')
            if self.fit_mode=='parametric':
                self.plot_results_coeff_qt()
                #self.plot_results_coeff_qt_eigen()
            if self.fit_mode=='parametric2D':
                self.plot_results_coeff_qt_2D()

        if 'norm' in save_plots:
            self.plot_results_norm_y_qt(var='resolution')        

        if 'qt-spectra' in save_plots:
            self.plot_spectra(A='qt',B='y')
        if 'y-spectra' in save_plots:
            self.plot_spectra(A='y', B='qt')

        if 'postfit' in save_plots:
            pulls  = (self.data-self.mu_postfit.reshape(self.data.shape))/np.sqrt(self.data)
            pulls2 = np.power(self.data-self.mu_postfit.reshape(self.data.shape), 2)/self.data 
            print 'Pulls sum = '+'{:3.3f}'.format(pulls.sum())
            print 'Pulls2 sum = '+'{:3.3f}'.format(pulls2.sum())
            self.save_template_snapshot(data=pulls, title='Postfit pulls', tag='postfit_pulls')
            self.save_template_snapshot(data=pulls2, title='Postfit pulls squared', tag='postfit_pulls2')

        # fill the tree
        self.out_tree.Fill()

        # reset parameters
        print 'Reset parameters to their initial values'
        for key,p in self.map_params.items():                        
            if '_true' in key:
                continue
            self.arglist[0] = p+1
            self.arglist[1] = self.map_params[key+'_true']  
            self.gMinuit.mnexcm( "SET PAR", self.arglist, 2, self.ierflg )

        return


    def propagate_covariance(self, ntoys=200):

        print "Running post-fit toys on alpha's to propagate error on beta's..."

        (val,err)  = ( ROOT.Double(0.), ROOT.Double(0.) )

        self.alpha = np.zeros(self.dim_alpha)
        for ia,a in enumerate(self.map_alphas):
            par_name = ''
            for key,p in self.map_params.items():
                if '_true' in key:
                    continue
                if p==a:
                    par_name = key
            self.gMinuit.GetParameter(self.map_params[par_name], val, err)
            self.alpha[ia] = val        

        self.Valpha = np.zeros((self.dim_alpha,self.dim_alpha))        
        for ia1,a1 in enumerate(self.map_alphas):
            par_name1 = ''
            for key1,p1 in self.map_params.items():
                if '_true' in key1:
                    continue
                if p1==a1:
                    par_name1 = key1            
            aa1 = self.map_params_minuit[par_name1] if self.map_params_minuit.has_key(par_name1) else -1
            for ia2,a2 in enumerate(self.map_alphas):
                par_name2 = ''
                for key2,p2 in self.map_params.items():
                    if '_true' in key2:
                        continue
                    if p2==a2:
                        par_name2 = key2
                aa2 = self.map_params_minuit[par_name2] if self.map_params_minuit.has_key(par_name2) else -1
                self.Valpha[ia1,ia2] = self.TMatrix(aa1,aa2) if aa1>=0 and aa2>=0 else 0.
        
        rnd_betas = np.zeros( (self.dim_beta, ntoys) )

        self.print_evals = False
        for itoy in range(ntoys):
            if itoy%50==0:
                print '\tToy ',itoy, '/', ntoys
            rnd_alpha = np.random.multivariate_normal(self.alpha,self.Valpha)
            if self.interpolation=='linear' and \
                    rnd_alpha[0]<self.mc_mass-0.499 or rnd_alpha[0]>self.mc_mass+0.499:
                continue
            rnd_par = [0.]*self.n_param
            for ia,a in enumerate(self.map_alphas):
                par_name = ''
                for key,p in self.map_params.items():
                    if '_true' in key:
                        continue
                    if p==a:
                        rnd_par[p] = rnd_alpha[ia]
            self.chi2(par=rnd_par)
            rnd_betas[:, itoy] = self.beta

        self.print_evals = False
        cov_betas = np.cov(rnd_betas)
        self.Vbeta_min += cov_betas


    # save the TTree and close
    def close(self):
        self.out_tree.Scan('minuit', '', '', 1)
        self.out_file.cd()
        self.out_tree.Write("tree", ROOT.TObject.kOverwrite)
        self.out_file.Close()

    # plotting functions
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

    def plot_results_norm_y_qt(self, var='resolution'):
        c = ROOT.TCanvas("canvas", "canvas", 600, 600) 
        c.SetRightMargin(0.15)
        h2 = ROOT.TH2F('norm_'+var, var+';|y|;q_{T} (GeV)',  
                       self.bins_y_size, array('f',self.bins_y), 
                       self.bins_qt_size, array('f',self.bins_qt) )
        h2.SetStats(0)         
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for iqt in range(self.bins_qt_size):
                qt_bin = self.get_qt_bin(iqt)
                fit_res = self.fit_results[y_bin+'_'+qt_bin+'_norm']
                val = 0.0
                if var=='resolution':
                    val = fit_res[2]/fit_res[3] if fit_res[3]>0. else 0.
                elif var=='pull':
                    val =  fit_res[4]
                h2.SetBinContent(iy+1,iqt+1, val)
        h2.Draw('COLZ')
        c.SaveAs(self.out_dir+'norm_'+var+'_'+self.job_name+'.png')
        c.SaveAs(self.out_dir+'norm_'+var+'_'+self.job_name+'.C')
        c.IsA().Destructor( c )                

    def plot_results_polynom_y(self, var='resolution'):
        max_order = 5
        histos = {}
        for icoeff,coeff in enumerate(self.coefficients):
            histos[coeff] = ROOT.TH2F('polynom_'+coeff+'_'+var, coeff+' '+var+';|y|;c_{n} (GeV^{-n})',  
                                      self.bins_y_size, array('f',self.bins_y), 
                                      max_order+1, array('f', range(max_order+2) ) )
            histos[coeff].SetStats(0)         

        c = ROOT.TCanvas("canvas", "canvas", 900, 600) 
        c.Divide(3,2)
        for icoeff,coeff in enumerate(self.coefficients):
            h2 = histos[coeff] 
            for iy in range(self.bins_y_size):
                y_bin = self.get_y_bin(iy)
                (valid_orders, order) = self.get_orders(coeff, y_bin)
                for o in range(order+1):
                    if o not in valid_orders:
                        continue                                                   
                    fit_res = self.fit_results[y_bin+'_'+coeff+'_'+'pol'+str(order)+'_p'+str(o)]
                    val = 0.
                    if var=='resolution':
                        val = fit_res[2]
                    elif var=='pull':
                        val =  fit_res[4]
                    h2.SetBinContent(iy+1,o+1, val)
            c.cd(icoeff+1)
            c.SetRightMargin(0.15)
            h2.Draw("COLZ")            
        
        c.cd()
        c.SaveAs(self.out_dir+'polynom_'+var+'_'+self.job_name+'.png')            
        c.SaveAs(self.out_dir+'polynom_'+var+'_'+self.job_name+'.C')            

    def plot_results_coeff_y_qt(self, var='resolution'):
        histos = {}
        for icoeff,coeff in enumerate(self.coefficients):
            histos[coeff] = ROOT.TH2F('coeff_'+coeff+var, coeff+' '+var+';|y|;q_{T}',  
                                      self.bins_y_size,  array('f',self.bins_y), 
                                      self.bins_qt_size, array('f',self.bins_qt), )
            histos[coeff].SetStats(0)         

        c = ROOT.TCanvas("canvas", "canvas", 1200, 800) 
        c.Divide(3,2)
        ROOT.gStyle.SetPaintTextFormat('4.3f')

        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for iqt in range(self.bins_qt_size):
                qt_bin = self.get_qt_bin(iqt)                               
                for icoeff,coeff in enumerate(self.coefficients):
                    val = 0.
                    if var=='resolution':
                        val = self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][2]
                    elif var=='value':
                        val = self.fit_results[y_bin+'_'+qt_bin+'_'+coeff][0]
                    histos[coeff].SetBinContent(iy+1,iqt+1, val  )

        for icoeff,coeff in enumerate(self.coefficients):
            c.cd(icoeff+1)
            if var=='resolution':
                ROOT.gPad.SetLogz()
                histos[coeff].SetMinimum(1e-03)
                histos[coeff].SetMaximum(5.)
            c.SetRightMargin(0.15)
            histos[coeff].Draw('COLZ' if var=='value' else 'COLZ')            
        
        c.cd()
        c.SaveAs(self.out_dir+'coeff_'+var+'_'+self.job_name+'.png')
        c.SaveAs(self.out_dir+'coeff_'+var+'_'+self.job_name+'.C')


    def plot_results_coeff_qt(self):
        from fit_utils import polynomial

        idx_beta = 0
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for icoeff,coeff in enumerate(self.coefficients):

                plt.figure()
                fig, ax = plt.subplots()

                x = np.array( [0.] + [ (self.bins_qt[iqt]+self.bins_qt[iqt+1])*0.5 for iqt in range(self.bins_qt_size)] + [self.bins_qt[-1]] )
                
                (valid_orders, order) = self.get_orders(coeff, y_bin)
                p = np.zeros(order+1)
                cov = np.zeros((order+1,order+1))
                idx_o1 = 0
                for o1 in range(order+1):
                    if o1 not in valid_orders:
                        continue         
                    idx_o2 = 0
                    for o2 in range(order+1):
                        if o2 not in valid_orders:
                            continue                             
                        cov[o1,o2] = self.Vbeta_min[idx_beta+idx_o1,idx_beta+idx_o2]
                        idx_o2 += 1
                    p[o1] = self.fit_results[y_bin+'_'+coeff+'_'+'pol'+str(order)+'_p'+str(o1)][0]                    
                    idx_o1 += 1
                idx_beta += len(valid_orders)

                ntoys = 200
                y_rnd = np.zeros( (x.size, ntoys) )                
                for itoy in range(ntoys):
                    p_rnd = np.random.multivariate_normal(p, cov)
                    y_rnd[:, itoy] = polynomial(x=x, coeff=p_rnd, order=order)

                ax.fill_between(x,  polynomial(x=x, coeff=p, order=order)-np.std(y_rnd, axis=1), 
                                polynomial(x=x, coeff=p, order=order)+np.std(y_rnd, axis=1), 
                                color='y', linestyle='-', label=r'$\pm1\sigma$')
                
                ax.plot(x, polynomial(x=x, coeff=p, order=order), 'r--', 
                        label=r'Fit ($\mathrm{pol}_{'+str(order)+'}$)', linewidth=3.0)

                ax.errorbar([ (self.bins_qt[iqt]+self.bins_qt[iqt+1])*0.5 for iqt in range(self.bins_qt_size)], 
                            self.res_coeff[coeff+'_'+y_bin+'_val'][0:self.bins_qt_size], 
                            xerr=[(self.bins_qt[iqt+1]-self.bins_qt[iqt])*0.5 for iqt in range(self.bins_qt_size)], 
                            yerr=self.res_coeff[coeff+'_'+y_bin+'_val_err'][0:self.bins_qt_size],
                            fmt='o', color='black', label='$'+coeff[0]+'_{'+coeff[1]+'}$ true')            

                plt.axis( [0.0, self.bins_qt[-1], -1, 1] )
                plt.grid(True)
                legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
                plt.title('Dataset: '+self.dataset_type+', $|y| \in ['+y_bin[1:5]+','+y_bin[7:11]+']$', fontsize=20)
                plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
                plt.ylabel('$'+coeff[0]+'_{'+coeff[1]+'}$', fontsize=20)
                plt.show()
                plt.savefig(self.out_dir+'/coefficient_'+coeff+'_'+y_bin+'_'+self.job_name+'_fit.png')
                plt.close('all')            

    def plot_results_coeff_qt_eigen(self):
        from fit_utils import polynomial
        
        # binning
        x = np.array( [0.] + [ (self.bins_qt[iqt]+self.bins_qt[iqt+1])*0.5 \
                                   for iqt in range(self.bins_qt_size)] + [self.bins_qt[-1]] )
        idx_beta = 0
        for iy in range(self.bins_y_size):
            y_bin = self.get_y_bin(iy)
            for icoeff,coeff in enumerate(self.coefficients):

                plt.figure()
                fig, ax = plt.subplots()

                (valid_orders, order) = self.get_orders(coeff, y_bin)
                
                # coefficients
                p = np.zeros(order+1)
                p_valid = np.zeros(len(valid_orders))
                cov = np.zeros((order+1,order+1))

                idx_o1 = 0
                for o1 in range(order+1):
                    if o1 not in valid_orders:
                        continue         
                    idx_o2 = 0
                    for o2 in range(order+1):
                        if o2 not in valid_orders:
                            continue                             
                        cov[o1,o2] = self.Vbeta_min[idx_beta+idx_o1,idx_beta+idx_o2]
                        idx_o2 += 1
                    p[o1] = self.fit_results[y_bin+'_'+coeff+'_'+'pol'+str(order)+'_p'+str(o1)][0]                    
                    p_valid[idx_o1] =  p[o1]
                    idx_o1 += 1
                
                cov_valid = self.Vbeta_min[idx_beta:(idx_beta+len(valid_orders)), 
                                           idx_beta:(idx_beta+len(valid_orders))]

                idx_beta += len(valid_orders)

                # eigenvalues/vectors of sub-cov matrix
                eigs =  np.linalg.eig(cov_valid)
                cov_valid_prime = eigs[0]
                U = eigs[1]
                # transformed of best-fit value of params
                p_valid_prime = np.dot(U.T, p_valid)

                ntoys = 200
                y_rnd = np.zeros( (x.size, ntoys) )                
                colors  = ['m', 'b', 'g', 'y', 'c']
                hatches = ['/', '/', '|', '|', '']
                for ieig in range(p_valid.size):
                    for itoy in range(ntoys):
                        # random value of the ieig'th eigenvector
                        p_valid_prime_rnd = copy.deepcopy(p_valid_prime)
                        p_valid_prime_rnd[ieig] = np.random.normal(p_valid_prime[ieig], math.sqrt(cov_valid_prime[ieig]) ) 
                        # transform back to p
                        p_valid_rnd = np.dot(U, p_valid_prime_rnd)
                        p_rnd = np.zeros(order+1)
                        idx_o = 0
                        for o in range(order+1):
                            if o not in valid_orders:
                                continue
                            p_rnd[o] = p_valid_rnd[idx_o]
                            idx_o += 1
                        y_rnd[:, itoy] = polynomial(x=x, coeff=p_rnd, order=order)
                    ax.fill_between(x,  
                                    y_rnd.mean(axis=1)-np.std(y_rnd, axis=1),
                                    y_rnd.mean(axis=1)+np.std(y_rnd, axis=1),
                                    facecolor='none' if ieig<p_valid.size-1 else colors[ieig] , 
                                    hatch=hatches[ieig] if ieig<p_valid.size-1 else None, 
                                    linestyle='-', 
                                    edgecolor=colors[ieig] if ieig<p_valid.size-1 else 'none',
                                    label=r'Eig. '+str(ieig+1))
                
                ax.plot(x, polynomial(x=x, coeff=p, order=order), 'r--', 
                        label=r'Fit ($\mathrm{pol}_{'+str(order)+'}$)', linewidth=3.0)
                ax.errorbar([ (self.bins_qt[iqt]+self.bins_qt[iqt+1])*0.5 for iqt in range(self.bins_qt_size)], 
                            self.res_coeff[coeff+'_'+y_bin+'_val'][0:self.bins_qt_size], 
                            xerr=[(self.bins_qt[iqt+1]-self.bins_qt[iqt])*0.5 for iqt in range(self.bins_qt_size)], 
                            yerr=self.res_coeff[coeff+'_'+y_bin+'_val_err'][0:self.bins_qt_size],
                            fmt='o', color='black', label='$'+coeff[0]+'_{'+coeff[1]+'}$ true')            

                plt.axis( [0.0, self.bins_qt[-1], -1, 1] )
                plt.grid(True)
                legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
                plt.title('Dataset: '+self.dataset_type+', $|y| \in ['+y_bin[1:5]+','+y_bin[7:11]+']$', fontsize=20)
                plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
                plt.ylabel('$'+coeff[0]+'_{'+coeff[1]+'}$', fontsize=20)
                plt.show()
                plt.savefig(self.out_dir+'/coefficient_'+coeff+'_'+y_bin+'_'+self.job_name+'_eigenvector_fit.png')
                plt.close('all')            


    def plot_results_coeff_qt_2D(self):
        from fit_utils import polynomial

        x = np.array( [0.] + [ (self.bins_qt[iqt]+self.bins_qt[iqt+1])*0.5 for iqt in range(self.bins_qt_size)] + [self.bins_qt[-1]] )

        for iy in range(self.bins_y_size):
            (y_bin, y) = ( self.get_y_bin(iy), self.mid_point_y(iy) )
             
            idx_beta = 0
            for icoeff,coeff in enumerate(self.coefficients):

                plt.figure()
                fig, ax = plt.subplots()

                (valid_orders_y, order_y)   = self.get_orders(coeff, '')
                (valid_orders_qt, order_qt) = self.get_orders(coeff, self.get_y_bin( 0 ))

                p = np.zeros( (order_y+1)*(order_qt+1) )
                cov = np.zeros( ((order_y+1)*(order_qt+1) , (order_y+1)*(order_qt+1)) )

                idx_o1 = 0
                for oy1 in range(order_y+1):
                    if oy1 not in valid_orders_y:
                        continue                             
                    for oqt1 in range(order_qt+1):
                        if oqt1 not in valid_orders_qt:
                            continue                         
                        idx_o2 = 0
                        for oy2 in range(order_y+1):
                            if oy2 not in valid_orders_y:
                                continue                    
                            for oqt2 in range(order_qt+1):
                                if oqt2 not in valid_orders_qt:
                                    continue   
                                cov[oy1*(order_qt+1) + oqt1 , oy2*(order_qt+1) + oqt2] = self.Vbeta_min[idx_beta+idx_o1,idx_beta+idx_o2]
                                idx_o2 += 1
                        p[oy1*(order_qt+1) + oqt1] = self.fit_results[coeff+'_'+'pol'+str(order_y)+'_pol'+str(order_qt)+'_p'+str(oy1)+'_p'+str(oqt1)][0]
                        idx_o1 += 1

                p_y = np.zeros( (order_qt+1) )
                for oy in range(order_y+1):
                    for oqt in range(order_qt+1):
                        p_y[oqt] += p[oy*(order_qt+1) + oqt]*math.pow(y, oy)

                ntoys = 200
                y_rnd = np.zeros( (x.size, ntoys) )                
                for itoy in range(ntoys):
                    p_rnd = np.random.multivariate_normal(p, cov)
                    p_y_rnd = np.zeros( (order_qt+1) )
                    for oy in range(order_y+1):
                        for oqt in range(order_qt+1):
                            p_y_rnd[oqt] += p_rnd[oy*(order_qt+1) + oqt]*math.pow(y, oy)
                    y_rnd[:, itoy] = polynomial(x=x, coeff=p_y_rnd, order=order_qt)

                ax.fill_between(x,  polynomial(x=x, coeff=p_y, order=order_qt)-np.std(y_rnd, axis=1), 
                                polynomial(x=x, coeff=p_y, order=order_qt)+np.std(y_rnd, axis=1), 
                                color='y', linestyle='-', label=r'$\pm1\sigma$')
                
                ax.plot(x, polynomial(x=x, coeff=p_y, order=order_qt), 'r--', 
                        label=r'Fit ($\mathrm{pol}_{'+str(order_qt)+'}$)', linewidth=3.0)

                ax.errorbar([ (self.bins_qt[iqt]+self.bins_qt[iqt+1])*0.5 for iqt in range(self.bins_qt_size)], 
                            self.res_coeff[coeff+'_'+y_bin+'_val'][0:self.bins_qt_size], 
                            xerr=[(self.bins_qt[iqt+1]-self.bins_qt[iqt])*0.5 for iqt in range(self.bins_qt_size)], 
                            yerr=self.res_coeff[coeff+'_'+y_bin+'_val_err'][0:self.bins_qt_size],
                            fmt='o', color='black', label='$'+coeff[0]+'_{'+coeff[1]+'}$ true')            

                plt.axis( [0.0, self.bins_qt[-1], -1, 1] )
                plt.grid(True)
                legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
                plt.title('Dataset: '+self.dataset_type+', $|y| \in ['+y_bin[1:5]+','+y_bin[7:11]+']$', fontsize=20)
                plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
                plt.ylabel('$'+coeff[0]+'_{'+coeff[1]+'}$', fontsize=20)
                plt.show()
                plt.savefig(self.out_dir+'/coefficient_'+coeff+'_'+y_bin+'_'+self.job_name+'_fit.png')
                plt.close('all')            

                idx_beta += (len(valid_orders_y)*len(valid_orders_qt))

    def plot_spectra(self, A='qt', B='y'):

        # plot A | B         
        for iB in range( getattr(self, 'bins_'+B+'_size')):
            B_bin = getattr(self,'get_'+B+'_bin')(iB)
            plt.figure()
            fig, ax = plt.subplots()
            vals = np.zeros( (getattr(self, 'bins_'+A+'_size')) )
            errs = np.zeros( (getattr(self, 'bins_'+A+'_size')) )
            trues = np.zeros( (getattr(self, 'bins_'+A+'_size')) )
            for iA in range(getattr(self, 'bins_'+A+'_size')):                    
                A_bin = getattr(self,'get_'+A+'_bin')(iA)
                bin_name = B_bin+'_'+A_bin+'_norm' if A=='qt' else  A_bin+'_'+B_bin+'_norm'
                vals[iA] = self.fit_results[bin_name][0]/getattr(self, 'width_'+A)(iA)
                errs[iA] = self.fit_results[bin_name][2]/getattr(self, 'width_'+A)(iA)
                trues[iA] = self.fit_results[bin_name][3]/getattr(self, 'width_'+A)(iA)
            ax.errorbar([ getattr(self, 'mid_point_'+A)(iA) for iA in range(getattr(self, 'bins_'+A+'_size'))], 
                        trues, 
                        xerr=[getattr(self, 'width_'+A)(iA)*0.5  for iA in range(getattr(self, 'bins_'+A+'_size'))], 
                        yerr=[0.]*getattr(self, 'bins_'+A+'_size'),
                        fmt='none', ecolor='r', elinewidth=3.0, label='True')            
            ax.errorbar([ getattr(self, 'mid_point_'+A)(iA) for iA in range(getattr(self, 'bins_'+A+'_size'))], 
                        vals, 
                        xerr=[getattr(self, 'width_'+A)(iA)*0.5  for iA in range(getattr(self, 'bins_'+A+'_size'))], 
                        yerr=errs,
                        fmt='o', color='black', label='Fit $\pm1\sigma$')            
            plt.axis( [0.0, getattr(self, 'bins_'+A)[-1], 0.0, vals.max()*1.1] )
            plt.grid(True)
            legend = ax.legend(loc='best', shadow=False, fontsize='x-large')
            if A=='qt':
                plt.title('Dataset: '+self.dataset_type+', $|y| \in ['+B_bin.split('_')[0][1:]+','+B_bin.split('_')[1][1:]+']$', fontsize=20)
                plt.xlabel('$q_{T}$ (GeV)', fontsize=20)
                plt.ylabel('$dN/dq_{T}$', fontsize=20)
            elif A=='y':
                print B_bin
                plt.title('Dataset: '+self.dataset_type+', $q_{T} \in ['+B_bin.split('_')[0][2:]+','+B_bin.split('_')[1][2:]+']$', fontsize=20)
                plt.xlabel('$|y|$', fontsize=20)
                plt.ylabel('$dN/d|y|$', fontsize=20)                
            plt.show()
            plt.savefig(self.out_dir+'/spectrum_'+A+'_'+B_bin+'_'+self.job_name+'.png')

                
    def plot_cov_matrix(self, n_free=0, cov=None, name=''):
        c = ROOT.TCanvas("canvas", "canvas", 600, 600) 
        h2 = ROOT.TH2F('cov', '', n_free, 0, n_free, n_free, 0, n_free)
        h2.SetStats(0) 
        for i in range(n_free):
            for j in range(n_free):
                if isinstance(cov, np.ndarray):
                    rho_ij = cov[i,j]/math.sqrt(cov[i,i]*cov[j,j]) \
                        if cov[i,i]>0.0 and cov[j,j]>0.0 else 0.0
                else:
                    rho_ij = cov(i,j)/math.sqrt(cov(i,i)*cov(j,j)) \
                        if cov(i,i)>0.0 and cov(j,j)>0.0 else 0.0                    
                h2.SetBinContent(i+1, j+1, rho_ij )
        h2.SetMinimum(-1.0)
        h2.SetMaximum(+1.0)
        h2.Draw("COLZ")
        c.SaveAs(self.out_dir+'covariance_'+name+'_'+self.job_name+'.png')
        c.SaveAs(self.out_dir+'covariance_'+name+'_'+self.job_name+'.C')
        c.IsA().Destructor( c )






