import ROOT
import math
import sys
import pickle
import time
from pprint import pprint
import numpy as np
import numpy.random as ran

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from array import array;
import copy

from template_parameters import accept_point, pdf_test, bin_ndarray

class Unfolder:

    def __init__(self, input_dir='../data/', params={},  rebin=(1,1), mass=80.0000, num_events=1000000, fix=[], interp_deg=1, n_points=500000, job_name='TEST', verbose=True, prior_coeff=0.3, prior_xsec=0.3, strategy=0, decorrelate=True, decorrelate_full=True,  do_semianalytic=True, do_taylor_expansion=False, n_taylor=2, add_constant_A4=True, run_minos=False, gen_toy=[0.0, 0.0, 0.0, 0.0, 0.0] ):

        self.run_minos=run_minos
        self.do_taylor_expansion = do_taylor_expansion
        self.n_taylor = n_taylor
        self.add_constant_A4 = add_constant_A4
        self.do_semianalytic = do_semianalytic
        self.decorrelate = decorrelate
        self.decorrelate_full = decorrelate_full
        self.input_dir = input_dir
        self.fix = fix
        self.mass = mass
        self.interp_deg = interp_deg
        self.num_events = num_events
        self.n_points = n_points
        self.verbose = verbose
        self.prior_coeff = prior_coeff
        self.prior_xsec = prior_xsec
        self.job_name = job_name
        self.strategy = strategy
        self.gen_toy = gen_toy
        self.truth = {}

        self.loads = [ [0.0, 0.0, 0.0, 0.0, 0.0],
                       [2.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0, 0.0, 0.0],
                       [0.0, 0.0, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0, 1.0, 0.0],
                       [0.0, 0.0, 0.0, 0.0, 2.0]
                       ]

        self.load_files(params=params, rebin=rebin)
        self.load_data()
        self.book_parameters(params=params)

        self.result = {}

        self.result['fix'] = self.fix
        self.result['n_points'] = n_points
        self.result['ntoys'] = 0
        self.result['num_events'] = num_events
        self.result['prior_coeff'] = prior_coeff
        self.result['prior_xsec'] = prior_xsec 
        self.f = open('result_'+job_name+'.pkl','wb')

    def load_files(self, params={}, rebin=(1,1)):

        self.input_y_bins  = np.linspace(params['params_lep']['y_range'][0], 
                                         params['params_lep']['y_range'][1],  
                                         params['params_lep']['y_bins']/rebin[0]+1)
        self.input_pt_bins = np.linspace(params['params_lep']['pt_range'][0], 
                                         params['params_lep']['pt_range'][1],  
                                         params['params_lep']['pt_bins']/rebin[1]+1)

        self.input_shapes_pt = params['params_template']['pt']
        self.input_shapes_y  = params['params_template']['y']
        self.input_shapes_A0  = params['params_template']['A0']
        self.input_shapes_A1  = params['params_template']['A1']
        self.input_shapes_A2  = params['params_template']['A2']
        self.input_shapes_A3  = params['params_template']['A3']
        self.input_shapes_A4  = params['params_template']['A4']
        self.input_shapes_mass  = params['params_template']['mass']

        self.num_pt_bins = len(self.input_shapes_pt)-1
        self.num_y_bins = len(self.input_shapes_y)-1

        if self.verbose:
            print "Loading files..."

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                for iA0,A0 in enumerate(self.input_shapes_A0):
                    for iA1,A1 in enumerate(self.input_shapes_A1):
                        for iA2,A2 in enumerate(self.input_shapes_A2):
                            for iA3,A3 in enumerate(self.input_shapes_A3):
                                for iA4,A4 in enumerate(self.input_shapes_A4):
                                    if not accept_point(coeff=[A0,A1,A2,A3,A4]):
                                        continue
                                    for im,m in enumerate(self.input_shapes_mass):
                                        input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(m)
                                        for ic,c in enumerate([A0,A1,A2,A3,A4]):
                                            input_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                                        grid_raw = np.load(self.input_dir+'/'+input_name+'.npy')
                                        grid = bin_ndarray( grid_raw, (grid_raw.shape[0]/rebin[0], grid_raw.shape[1]/rebin[1]), 'sum' )
                                        norm = grid.sum() 
                                        setattr(self, input_name, grid)
                                        setattr(self, input_name+'_norm', norm)
                                        self.shape = grid.shape                                         

    # generate a rnd sample
    def load_data(self):
        self.toy_data()
        if self.decorrelate and not self.decorrelate_full:
            self.prefit_covariance()
        elif self.decorrelate_full:
            self.prefit_covariance_full()

    # generate rnd samples
    def toy_data(self, ntoy=0):

        self.ntoy = ntoy
        self.data = np.zeros(self.shape)
        self.data_asymov = np.zeros(self.shape)

        normalisation = 0.0
        #total_weight = 0.0
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                #total_weight += pdf_test(pt=pt_bin[0], y=y_bin[0])
                input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate(self.gen_toy):
                    input_name += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                for pt in np.linspace(pt_bin[0],pt_bin[1]-1,4):
                    normalisation += getattr(self, input_name+'_norm')*pdf_test(pt, y=y_bin[0])
        
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                weight = 0.0
                for pt in np.linspace(pt_bin[0],pt_bin[1]-1,4):
                    weight += pdf_test(pt, y=y_bin[0]) #/total_weight
                name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
                input_name = 'mixed_dataset_'+name+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate(self.gen_toy):
                    input_name += '_A'+str(ic)+('{:03.2f}'.format(c))

                n_true = self.num_events*weight/normalisation
                data_asymov = getattr(self, input_name)*n_true 
                data_rnd = np.random.poisson( data_asymov ) #* (max(np.random.normal(1.0, self.prior_xsec), 0.0))

                self.truth[name] = n_true
                self.truth[name+'_gen'] = data_rnd.sum()/getattr(self, input_name+'_norm')
                for ic,c in enumerate(self.gen_toy):
                    self.truth[name+'_A'+str(ic)] = c
                self.data += data_rnd
                self.data_asymov += data_asymov


        self.truth['mass'] = self.mass
        self.truth['num_events_toy'] = self.data.sum()
        if ntoy==0:
            xx, yy = np.meshgrid(self.input_pt_bins, self.input_y_bins)        
            plt.pcolormesh(yy, xx, self.data)
            plt.colorbar()
            plt.show()
            plt.savefig('data_toy_'+str(ntoy)+'_'+self.job_name+'.png')
            plt.close()
        #print self.data
        
    def book_parameters(self, params={}):

        self.ndof = self.data.size
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

        self.map_params={}

        # add mass
        self.n_param = 0
        self.gMinuit.DefineParameter( self.n_param, "mass", self.mass, 0.020, self.mass-0.490, self.mass+0.490 )        
        self.map_params['mass'] = self.n_param
        if 'mass' in self.fix:
            self.gMinuit.FixParameter(0)
        else:
            self.ndof -= 1
        self.n_param += 1

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                par_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
                self.map_params[par_name] = self.n_param

                # 1/3 -- 3 with 10 step
                if not (self.decorrelate or self.decorrelate_full):
                    self.gMinuit.DefineParameter( self.n_param, par_name, self.truth[par_name], 
                                                  100, 
                                                  self.truth[par_name]/2.0, 
                                                  self.truth[par_name]*2.0  )
                elif self.decorrelate:
                    index = ipt*(len(self.input_shapes_y)-1) + iy
                    scale_range = self.sigmaU[index]*(2.0 if self.do_semianalytic else 1.0)
                    self.gMinuit.DefineParameter( self.n_param, par_name, 
                                                  self.truth[par_name+'_prime'], 
                                                  scale_range*0.5, 
                                                  self.truth[par_name+'_prime']-scale_range*20, self.truth[par_name+'_prime']+scale_range*20 )
                elif self.decorrelate_full:
                    index = (ipt*(len(self.input_shapes_y)-1) + iy)*len(self.loads)
                    scale_range = self.sigmaU_full[index]
                    self.gMinuit.DefineParameter( self.n_param, par_name, 
                                                  self.truth[par_name+'_prime'], 
                                                  scale_range*1.0, 
                                                  self.truth[par_name+'_prime']-scale_range*10, self.truth[par_name+'_prime']+scale_range*10 )

                if 'pt_y' in self.fix: 
                    self.gMinuit.FixParameter(self.n_param)
                else:
                    self.ndof -= 1

                self.n_param += 1

                for icoeff,coeff in enumerate(['A0', 'A1', 'A2', 'A3', 'A4']):
                # add one parameter per pt/y/Ai bin                    
                    par_name_coeff = par_name+'_'+coeff 
                    self.map_params[par_name_coeff] = self.n_param            

                    if not self.decorrelate_full:
                        self.gMinuit.DefineParameter( self.n_param, par_name_coeff, self.truth[par_name_coeff], 0.01, self.truth[par_name_coeff]-0.5, self.truth[par_name_coeff]+0.5)
                    else:
                        index = (ipt*(len(self.input_shapes_y)-1) + iy)*len(self.loads) + icoeff + 1
                        scale_range = self.sigmaU_full[index]
                        self.gMinuit.DefineParameter( self.n_param, par_name_coeff, 
                                                      self.truth[par_name_coeff+'_prime'], 
                                                      scale_range*1.0, 
                                                      self.truth[par_name_coeff+'_prime']-scale_range*5, self.truth[par_name_coeff+'_prime']+scale_range*5 )                        
                    if coeff in self.fix:
                        self.gMinuit.FixParameter(self.n_param)
                    if (coeff not in self.fix) or self.do_semianalytic:
                        if self.do_taylor_expansion and ipt==0:
                            self.ndof -= (self.n_taylor)
                            if self.add_constant_A4 and coeff=='A4':
                                self.ndof -= 1

                        elif not self.do_taylor_expansion:
                            self.ndof -= 1

                    self.n_param += 1            

        # dimension of Ai parameter vector
        self.dim_A = (len(self.input_shapes_pt)-1)*(len(self.input_shapes_y)-1)*(len(self.loads)-1)

        self.dim_beta = self.n_taylor*(len(self.input_shapes_y)-1)*(len(self.loads)-1)
        if self.add_constant_A4:
            self.dim_beta = self.n_taylor*(len(self.input_shapes_y)-1)*(len(self.loads)-2) + (self.n_taylor+1)*(len(self.input_shapes_y)-1)

        if self.do_semianalytic:
            self.theta = np.zeros((self.dim_A))
            self.theta_err = np.diag( np.ones(self.dim_A) )
        if self.do_taylor_expansion:
            self.beta = np.zeros((self.dim_beta))
            self.beta_err = np.diag( np.ones(self.dim_beta) )
            


    def prefit_covariance(self):
        self.dim = (len(self.input_shapes_pt)-1)*(len(self.input_shapes_y)-1)
        covinv = np.zeros((self.dim,self.dim))
        x =  np.zeros(self.dim)
        # first index
        for ipt1 in range(len(self.input_shapes_pt)-1):
            pt_bin1=[ self.input_shapes_pt[ipt1], self.input_shapes_pt[ipt1+1] ]
            for iy1 in range(len(self.input_shapes_y)-1):
                y_bin1=[ self.input_shapes_y[iy1], self.input_shapes_y[iy1+1] ]                
                index1 = ipt1*(len(self.input_shapes_y)-1) + iy1
                name1 = 'pt{:02.1f}'.format(pt_bin1[0])+'-'+'{:02.1f}'.format(pt_bin1[1])+'_'+'y{:03.2f}'.format(y_bin1[0])+'-'+'{:03.2f}'.format(y_bin1[1])
                input_name1 = 'mixed_dataset_'+name1+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate(self.loads[0]):
                    input_name1 += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                tj = copy.deepcopy(getattr(self, input_name1))
                self.ignore_from_sum(tj)
                x[index1] = self.truth[name1]

                # second index
                for ipt2 in range(len(self.input_shapes_pt)-1):
                    pt_bin2=[ self.input_shapes_pt[ipt2], self.input_shapes_pt[ipt2+1] ]
                    for iy2 in range(len(self.input_shapes_y)-1):
                        y_bin2=[ self.input_shapes_y[iy2], self.input_shapes_y[iy2+1] ]                
                        index2 = ipt2*(len(self.input_shapes_y)-1) + iy2
                        name2 = 'pt{:02.1f}'.format(pt_bin2[0])+'-'+'{:02.1f}'.format(pt_bin2[1])+'_'+'y{:03.2f}'.format(y_bin2[0])+'-'+'{:03.2f}'.format(y_bin2[1])
                        input_name2 = 'mixed_dataset_'+name2+'_M'+'{:05.3f}'.format(self.mass)
                        for ic,c in enumerate(self.loads[0]):
                            input_name2 += '_A'+str(ic)+('{:03.2f}'.format(c))                                                
                        tk = copy.deepcopy(getattr(self, input_name2))
                        self.ignore_from_sum(tk)

                        asymov = copy.deepcopy(self.data_asymov)
                        # make sure there are no 0.0
                        if not asymov.all():
                            asymov += np.ones(self.shape)*sys.float_info.min

                        covinv[index1][index2] = (tj*tk/asymov).sum()

        covinv *= self.data_asymov.sum()
        cov = np.linalg.inv(covinv)
        np.set_printoptions(threshold=np.inf)
        #print cov
        res =  np.linalg.eig(cov)
        #print res
        self.U = res[1]
        self.sigmaU = np.sqrt(res[0]*self.data_asymov.sum())
        #print x.sum()
        x_prime = np.dot((self.U).T, x)

        for ipt1 in range(len(self.input_shapes_pt)-1):
            pt_bin1=[ self.input_shapes_pt[ipt1], self.input_shapes_pt[ipt1+1] ]
            for iy1 in range(len(self.input_shapes_y)-1):
                y_bin1=[ self.input_shapes_y[iy1], self.input_shapes_y[iy1+1] ]                
                index1 = ipt1*(len(self.input_shapes_y)-1) + iy1
                name1 = 'pt{:02.1f}'.format(pt_bin1[0])+'-'+'{:02.1f}'.format(pt_bin1[1])+'_'+'y{:03.2f}'.format(y_bin1[0])+'-'+'{:03.2f}'.format(y_bin1[1])
                self.truth[name1+"_prime"] = x_prime[index1]
                print name1+"_prime", x_prime[index1]
        #print x_prime
        #print np.dot(U, x_prime) - x
        #print np.dot(np.dot(U.T, cov), U)


    def ignore_from_sum(self, a):
        np.place(a, self.data_asymov<=0., 0.0)      

    def prefit_covariance_full(self):
        self.dim = (len(self.input_shapes_pt)-1)*(len(self.input_shapes_y)-1)*len(self.loads)
        covinv = np.zeros((self.dim,self.dim))
        x =  np.zeros(self.dim)

        # first index
        for ipt1 in range(len(self.input_shapes_pt)-1):
            pt_bin1=[ self.input_shapes_pt[ipt1], self.input_shapes_pt[ipt1+1] ]
            for iy1 in range(len(self.input_shapes_y)-1):
                y_bin1=[ self.input_shapes_y[iy1], self.input_shapes_y[iy1+1] ]                
                name1 = 'pt{:02.1f}'.format(pt_bin1[0])+'-'+'{:02.1f}'.format(pt_bin1[1])+'_'+'y{:03.2f}'.format(y_bin1[0])+'-'+'{:03.2f}'.format(y_bin1[1])
                input_name_all1   = 'mixed_dataset_'+name1+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate(self.loads[0]):
                    input_name_all1 += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                tj = copy.deepcopy(getattr(self, input_name_all1))
                self.ignore_from_sum(tj)

                for iload1,load1 in enumerate(self.loads):
                    input_name_coeff1 = 'mixed_dataset_'+name1+'_M'+'{:05.3f}'.format(self.mass)
                    for ic,c in enumerate( load1 ):
                        input_name_coeff1 += '_A'+str(ic)+('{:03.2f}'.format(c))
                    tjC = copy.deepcopy(getattr(self, input_name_coeff1))
                    tjC -= tj
                    tjC /= (load1[iload1-1] if iload1>0 else 1.0)
                    self.ignore_from_sum(tjC)

                    index1 =  (ipt1*(len(self.input_shapes_y)-1) + iy1)*len(self.loads) + iload1
                    par_name1 = name1+('' if iload1==0 else '_A'+str(iload1-1))
                    par_name_norm1 = name1
                    #print index1 , '-->', par_name1, '=', self.truth[par_name1]
                    x[index1] = self.truth[par_name1]

                    # second index
                    for ipt2 in range(len(self.input_shapes_pt)-1):
                        pt_bin2=[ self.input_shapes_pt[ipt2], self.input_shapes_pt[ipt2+1] ]
                        for iy2 in range(len(self.input_shapes_y)-1):
                            y_bin2=[ self.input_shapes_y[iy2], self.input_shapes_y[iy2+1] ]                
                            name2 = 'pt{:02.1f}'.format(pt_bin2[0])+'-'+'{:02.1f}'.format(pt_bin2[1])+'_'+'y{:03.2f}'.format(y_bin2[0])+'-'+'{:03.2f}'.format(y_bin2[1])
                            input_name_all2   = 'mixed_dataset_'+name2+'_M'+'{:05.3f}'.format(self.mass)
                            for ic,c in enumerate(self.loads[0]):
                                input_name_all2 += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                            tk = copy.deepcopy(getattr(self, input_name_all2))
                            self.ignore_from_sum(tk)

                            for iload2,load2 in enumerate(self.loads):
                                input_name_coeff2 = 'mixed_dataset_'+name2+'_M'+'{:05.3f}'.format(self.mass)
                                for ic,c in enumerate( load2 ):
                                    input_name_coeff2 += '_A'+str(ic)+('{:03.2f}'.format(c))
                                tkD = copy.deepcopy(getattr(self, input_name_coeff2))
                                tkD -= tk
                                tkD /= (load2[iload2-1] if iload2>0 else 1.0)
                                self.ignore_from_sum(tkD)

                                index2 =  (ipt2*(len(self.input_shapes_y)-1) + iy2)*len(self.loads) + iload2
                                par_name2 = name2+('' if iload2==0 else '_A'+str(iload2-1))
                                par_name_norm2 = name2
                                #print '\t', index2 , '-->', par_name2, '=', self.truth[par_name2]

                                asymov = copy.deepcopy(self.data_asymov)
                                # make sure there are no 0.0
                                if not asymov.all():
                                    asymov += np.ones(self.shape)*sys.float_info.min

                                covinv[index1][index2] = 0.0
                                if iload1==0 and iload2==0:
                                    covinv[index1][index2] = (tj*tk/asymov).sum()
                                elif iload1==0 and iload2>0:
                                    covinv[index1][index2] = (tj*tkD/asymov).sum()*self.truth[par_name_norm2]
                                elif iload1>0 and iload2==0:
                                    covinv[index1][index2] = (tjC*tk/asymov).sum()*self.truth[par_name_norm1]
                                else: 
                                    covinv[index1][index2] = (tjC*tkD/asymov).sum()*self.truth[par_name_norm1]*self.truth[par_name_norm2]
                                #print '\t\t (', index1,index2, ') ===>', covinv[index1][index2]

        np.set_printoptions(threshold=np.inf)
        #print covinv

        cov = np.linalg.inv(covinv)
        #print cov
        res =  np.linalg.eig(cov)
        #U = res[1]
        #print res

        self.U_full = res[1]
        self.sigmaU_full = np.sqrt(res[0])
        print self.sigmaU_full
        x_prime = np.dot((self.U_full).T, x)

        for ipt1 in range(len(self.input_shapes_pt)-1):
            pt_bin1=[ self.input_shapes_pt[ipt1], self.input_shapes_pt[ipt1+1] ]
            for iy1 in range(len(self.input_shapes_y)-1):
                y_bin1=[ self.input_shapes_y[iy1], self.input_shapes_y[iy1+1] ]                
                name1 = 'pt{:02.1f}'.format(pt_bin1[0])+'-'+'{:02.1f}'.format(pt_bin1[1])+'_'+'y{:03.2f}'.format(y_bin1[0])+'-'+'{:03.2f}'.format(y_bin1[1])
                for iload1,load1 in enumerate(self.loads):
                    index1 =  (ipt1*(len(self.input_shapes_y)-1) + iy1)*len(self.loads) + iload1
                    par_name1 = name1+('' if iload1==0 else '_A'+str(iload1-1))
                    self.truth[par_name1+"_prime"] = x_prime[index1]
                    print par_name1+"_prime", x_prime[index1]

        #x_prime = np.dot( U.T, x )
        #print np.dot(U, x_prime) - x
                    

    def chi2_analytic(self, par):

        # find closest teplates for interpolation
        mass = par[0]
        mass_neighbours = self.find_neighbours(mass=mass)
        if len(mass_neighbours)<2:
            print "No points for mass interpolation at m=%s" % mass
            return 1.0

        x = np.zeros(len(mass_neighbours))
        grid_points_A = []
        grid_points_b = []

        n = self.data.flatten()
        if n[n<10].size > 0:
            print "Bins with less than 10 entries!"

        Vinv = np.diag(1./n)

        # prior
        priors = np.ones(self.dim_A)*self.prior_coeff
        Vinv_prior = np.diag(1./priors)
        theta_prior = np.zeros(self.dim_A)

        for iim,im in enumerate(mass_neighbours):
            mass_point = self.input_shapes_mass[im]
            x[iim] = mass_point

            Am = np.zeros((self.data.size, self.dim_A))
            bm = np.zeros(n.shape)

            # first index
            for ipt in range(len(self.input_shapes_pt)-1):
                pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
                for iy in range(len(self.input_shapes_y)-1):
                    y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                    name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
                    input_name_all   = 'mixed_dataset_'+name+'_M'+'{:05.3f}'.format(mass_point)
                    for ic,c in enumerate( self.loads[0] ):
                        input_name_all += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                    tj = copy.deepcopy(getattr(self, input_name_all))
                    #self.ignore_from_sum(tj)

                    for iload,load in enumerate(self.loads):
                        index     =  (ipt*(len(self.input_shapes_y)-1) + iy)*(len(self.loads)-1) + (iload-1)
                        index_par =  (ipt*(len(self.input_shapes_y)-1) + iy)*(len(self.loads))
                        if iload==0:
                            bm += tj.flatten()*par[index_par+1]
                            continue
                        input_name_coeff = 'mixed_dataset_'+name+'_M'+'{:05.3f}'.format(mass_point)
                        for ic,c in enumerate( load ):
                            input_name_coeff += '_A'+str(ic)+('{:03.2f}'.format(c))
                        tjC = copy.deepcopy(getattr(self, input_name_coeff))
                        tjC -= tj
                        tjC /= load[iload-1]
                        #self.ignore_from_sum(tjC)

                        for i in range(n.size):
                            Am[i][index] = tjC.flatten()[i]*par[index_par+1] 

            #print "=>", mass_point, Am .sum()
            grid_points_A.append( Am )
            grid_points_b.append( bm )

        nprime =  n-self.interpolate_mass(mass=mass, x=x, grid_points=grid_points_b, deg=self.interp_deg, shape=(self.data.size))
        A = self.interpolate_mass(mass=mass, x=x, grid_points=grid_points_A, deg=self.interp_deg, shape=(self.data.size, self.dim_A))
        #nprime = n-b
        #np.set_printoptions(threshold=np.inf)
        #print A.shape
        #print "A:", A
        #print "V-1:", Vinv
        #print "b:", b
        #print "n':", nprime

        aux1 =  np.linalg.multi_dot([A.T, Vinv, A])
        if self.prior_coeff>0.0:
            aux1 += Vinv_prior
        aux1_inv = np.linalg.inv(aux1)
        aux2 = np.linalg.multi_dot( [ A.T, Vinv, nprime ] )
        aux3 = np.linalg.multi_dot( [ Vinv_prior, theta_prior ])
        if self.prior_coeff>0.0:
            aux2 += aux3
        theta = np.linalg.multi_dot( [aux1_inv, aux2 ] )
        #theta = theta_prior
        self.theta = theta

        # FIX!!!
        self.theta_err = aux1_inv

        res1 = (nprime - np.dot(A, theta))
        res2 = (theta_prior - theta)
        chi2min = np.linalg.multi_dot( [res1.T, Vinv, res1 ] )
        if self.prior_coeff>0.0:
            chi2min += np.linalg.multi_dot( [res2.T, Vinv_prior, res2] )
        return chi2min


    def chi2_analytic_parametric(self, par):

        mass = par[0]
        mass_neighbours = self.find_neighbours(mass=mass)
        if len(mass_neighbours)<2:
            print "No points for mass interpolation at m=%s" % mass
            return 1.0

        x = np.zeros(len(mass_neighbours))
        grid_points_A = []
        grid_points_b = []

        n = self.data.flatten()
        if n[n<10].size > 0:
            print "Bins with less than 10 entries!"

        Vinv = np.diag(1./n)

        K = np.zeros( ( self.dim_A, self.dim_beta) )

        # prior
        priors = np.ones( self.dim_beta )
        for i in range(self.dim_beta):
            if i%3==0:
                priors[i] *= 1e-03
            elif i%3==1:
                priors[i] *= 1e-04
            elif i%3==2:
                priors[i] *= 1e-05
        Vinv_prior = np.diag(1./priors)
        beta_prior = np.zeros(self.dim_beta)

        for iim,im in enumerate(mass_neighbours):
            mass_point = self.input_shapes_mass[im]
            x[iim] = mass_point

            Am = np.zeros((self.data.size, self.dim_A))
            bm = np.zeros(n.shape)

            # first index
            for ipt in range(len(self.input_shapes_pt)-1):
                pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
                for iy in range(len(self.input_shapes_y)-1):
                    y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                    name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
                    input_name_all   = 'mixed_dataset_'+name+'_M'+'{:05.3f}'.format(mass_point)
                    for ic,c in enumerate( self.loads[0] ):
                        input_name_all += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                    tj = copy.deepcopy(getattr(self, input_name_all))

                    for iload,load in enumerate(self.loads):
                        index     =  (ipt*(len(self.input_shapes_y)-1) + iy)*(len(self.loads)-1) + (iload-1)
                        index_par =  (ipt*(len(self.input_shapes_y)-1) + iy)*(len(self.loads))
                        if iload==0:
                            bm += tj.flatten()*par[index_par+1]
                            continue

                        if not self.add_constant_A4:
                            # all coefficients
                            for e in range(self.n_taylor):
                                index_e = (iy*(len(self.loads)-1) + (iload-1))*self.n_taylor + e
                                power = e+2
                                K[index, index_e] = math.pow( 0.5*(pt_bin[0]+pt_bin[1]), power)/math.factorial(power)
                        else:
                            # A4
                            if iload==5:
                                for e in range(self.n_taylor+1):
                                    index_e = iy*( (len(self.loads)-2)*self.n_taylor + 1*(self.n_taylor+1) ) + (len(self.loads)-2)*self.n_taylor + e
                                    power = e+1 if e>0 else 0
                                    K[index, index_e] = math.pow( 0.5*(pt_bin[0]+pt_bin[1]), power)/math.factorial(power)                          
                            # all others
                            else:
                                for e in range(self.n_taylor):
                                    index_e = iy*( (len(self.loads)-2)*self.n_taylor + 1*(self.n_taylor+1) ) + (iload-1)*self.n_taylor + e
                                    power = e+2
                                    K[index, index_e] = math.pow( 0.5*(pt_bin[0]+pt_bin[1]), power)/math.factorial(power)


                        input_name_coeff = 'mixed_dataset_'+name+'_M'+'{:05.3f}'.format(mass_point)
                        for ic,c in enumerate( load ):
                            input_name_coeff += '_A'+str(ic)+('{:03.2f}'.format(c))
                        tjC = copy.deepcopy(getattr(self, input_name_coeff))
                        tjC -= tj
                        tjC /= load[iload-1]

                        for i in range(n.size):
                            Am[i][index] = tjC.flatten()[i]*par[index_par+1] 

            grid_points_A.append( Am )
            grid_points_b.append( bm )

        nprime =  n-self.interpolate_mass(mass=mass, x=x, grid_points=grid_points_b, deg=self.interp_deg, shape=(self.data.size))
        A = self.interpolate_mass(mass=mass, x=x, grid_points=grid_points_A, deg=self.interp_deg, shape=(self.data.size, self.dim_A))

        aux1 =  np.linalg.multi_dot([K.T, A.T, Vinv, A, K])
        if self.prior_coeff>0.0:
            aux1 += Vinv_prior
        aux1_inv = np.linalg.inv(aux1)
        aux2 = np.linalg.multi_dot( [ K.T, A.T, Vinv, nprime ] )
        aux3 = np.linalg.multi_dot( [ Vinv_prior, beta_prior ])
        if self.prior_coeff>0.0:
            aux2 += aux3
        beta = np.linalg.multi_dot( [aux1_inv, aux2 ] )        

        self.beta = beta
        # FIX!!!
        self.beta_err = aux1_inv

        res1 = (nprime - np.linalg.multi_dot( [A, K, beta]) )
        res2 = (beta_prior -beta)
        chi2min = np.linalg.multi_dot( [res1.T, Vinv, res1 ] )
        if self.prior_coeff>0.0:
            chi2min += np.linalg.multi_dot( [res2.T, Vinv_prior, res2] )

        return chi2min



    def reset_parameters(self):
        self.arglist[0] = 1
        self.arglist[1] = self.mass 
        self.gMinuit.mnexcm("SET PAR", self.arglist, 2, self.ierflg )

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                par_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
                self.arglist[0] = self.map_params[par_name]+1
                if not (self.decorrelate or self.decorrelate_full):
                    self.arglist[1] = self.truth[par_name] 
                else:
                    self.arglist[1] = self.truth[par_name+'_prime'] 

                self.gMinuit.mnexcm("SET PAR", self.arglist, 2, self.ierflg )

                for icoeff,coeff in enumerate(['A0', 'A1', 'A2', 'A3', 'A4']):
                    par_name_coeff = par_name+'_'+coeff 
                    self.arglist[0] = self.map_params[par_name_coeff]+1
                    if not self.decorrelate_full:
                        self.arglist[1] = self.truth[par_name_coeff] 
                    else:
                        self.arglist[1] = self.truth[par_name_coeff+'_prime'] 

                    self.gMinuit.mnexcm("SET PAR", self.arglist, 2, self.ierflg )


    def find_neighbours(self, mass=80.000):
        mass_neighbours=[]

        for im in range(len(self.input_shapes_mass)-1):            
            mass_low = self.input_shapes_mass[im]
            mass_high = self.input_shapes_mass[im+1]
            if mass>=mass_low and mass<mass_high:
                if self.interp_deg==1:
                    mass_neighbours.extend( [im,im+1] )
                elif self.interp_deg==2:
                    if im>len(self.input_shapes_mass)-3:
                        mass_neighbours.extend( [im-1,im,im+1] )
                    else:
                        mass_neighbours.extend( [im,im+1,im+2] )
                elif self.interp_deg==3:
                    if im>len(self.input_shapes_mass)-3:
                        mass_neighbours.extend( [im-2,im-1,im,im+1] )
                    else:
                        mass_neighbours.extend( [im-1,im,im+1,im+2] )
                else:
                    return mass_neighbours
        return mass_neighbours

    def convert_to_norm(self, par):
        #print "convert_to_norm..."
        x_prime = []
        for p in range(self.n_param):
            if (p-1)%(len(self.loads[0])+1) == 0:
                x_prime.append(par[p])
                #print "Append par" , p, "with value", par[p]
        x = np.dot(self.U, x_prime) 
        counter = 0
        for p in range(self.n_param):
            if (p-1)%(len(self.loads[0])+1) == 0:
                par[p] = x[counter]
                counter += 1
                #print "Set new par" , p, "with value", par[p]

    def convert_to_norm_full(self, par):
        #print "convert_to_norm_full..."
        x_prime = []
        for p in range(self.n_param):
            if p>0:
                x_prime.append(par[p])
                #print "Append par" , p, "with value", par[p]
        x = np.dot(self.U_full, x_prime) 
        counter = 0
        for p in range(self.n_param):
            if p>0:
                par[p] = x[counter]
                counter += 1
                #print "Set new par" , p, "with value", par[p]


    # compute chi2
    def chi2(self, par):
        nll = 0.0
        # parameters
        if self.decorrelate:
            self.convert_to_norm(par)
        if self.do_taylor_expansion:
            nll += self.chi2_analytic_parametric(par)
        else:
            nll += self.chi2_analytic(par)
        print "chi2/ndof:",  nll/self.ndof 
        return nll


    # compute -2*log likelihood
    def nll(self, par):

        nll = 0.0

        #print "Evaluating fcn..."
        pdf = np.zeros(self.shape)

        # parameters
        if self.decorrelate and not self.decorrelate_full:
            self.convert_to_norm(par)
        elif self.decorrelate_full:
            self.convert_to_norm_full(par)

        # find closest teplates for interpolation
        mass = par[0]
        mass_neighbours = self.find_neighbours(mass=mass)
        if len(mass_neighbours)<2:
            print "No points for mass interpolation at m=%s" % mass
            return 1.0

        # keet track of par mapping: par[0] is mass
        count_param = 1
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                

                # the normalisation of a given pt_y bin
                norm_bin = par[count_param]

                # the pdf of a given pt_y bin
                pdf_tmp_pt_y = np.zeros(self.shape)

                load_name = 'mixed_dataset_pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])

                # mass interpolation and Ai averaging
                x = np.zeros(len(mass_neighbours))
                grid_points = []
                for iim,im in enumerate(mass_neighbours):
                    mass_point = self.input_shapes_mass[im]
                    x[iim] = mass_point

                    # the pdf at a given pt_y bin is obtained by linearly adding the pdf's at various Ai
                    pdf_tmp_coeff = np.zeros(self.shape)

                    for iload, load in enumerate(self.loads):
                        grid_name = load_name + '_M'+'{:05.3f}'.format(mass_point) 
                        for ic,c in enumerate( load ):
                            grid_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                        grid = getattr(self, grid_name)
                        # this is [0,0,0,0,0]
                        if iload==0:
                            norm = 1.0
                            for p in range(0, len(load)):
                                shift_gen = self.loads[p+1][p]
                                norm -= par[count_param + (p+1)]/shift_gen
                            pdf_tmp_coeff += norm*grid
                        # these are [2,0,0,0,0], [0,1,0,0,0], ...
                        else:
                            shift_gen = load[iload-1]
                            pdf_tmp_coeff += par[count_param+iload]/shift_gen*grid                            

                    # avoid negative values
                    grid_points.append( pdf_tmp_coeff )

                # interpolate templates of different mass                
                pdf_tmp_pt_y += self.interpolate_mass(mass=mass, x=x, grid_points=grid_points, deg=self.interp_deg, shape=pdf.shape)
                # normalise bin using norm_bin
                pdf_tmp_pt_y *= norm_bin
                # add tis bin to the total pdf
                pdf += pdf_tmp_pt_y
                # increase par counter by 3 units: (pt_y, A0, A1, ..., A4)
                count_param += (1+len(self.loads[0]))


        # prevent from any zeros in the pdf
        pdf += (np.ones(self.shape)*sys.float_info.min)
        pdf = np.absolute(pdf)

        nll += 2*(-self.data*np.log(pdf) + pdf).sum()

        # prior
        if self.prior_coeff>0.0 or self.prior_xsec>0.0:
            for key,p in self.map_params.items():            
                if self.prior_xsec>0.0 and ('_A' not in key and 'mass' not in key):                
                    true = self.truth[key]
                    if self.decorrelate:
                        true = self.truth[key+'_prime']
                    err = self.prior_xsec
                    nll += math.pow((par[p]-true)/(err*true), 2.0)
                elif self.prior_coeff>0.0 and ('_A' in key):
                    fit = par[p]
                    true = self.truth[key]
                    err = self.prior_coeff
                    nll += math.pow((fit-true)/err, 2.0)        

        print "nll:", nll
        return nll


    def fcn(self, npar, gin, f, par, iflag ):

        nll = 0.0

        if self.do_semianalytic:
            nll = self.chi2(par)
        else:
            nll = self.nll(par)

        #print("{:16f}".format(nll))
        #for key,p in self.map_params.items():          
        #    #if '_A' not in key and 'mass' not in key:                
        #    print '\t', key, par[p], self.truth[key]

        f[0] = nll
        return

    def interpolate_mass(self, mass=80.000, x=np.array([]), grid_points=[], deg=1, shape=()):

        #print "Interpolate %s in range %s" % (mass, x[:])        
        if deg==1:
            fit = (mass-x[0])/(x[1]-x[0])
            return (1-fit)*grid_points[0] + fit*grid_points[1]
        else:
            y = np.zeros( (len(grid_points), grid_points[0].size) )
            for m in range(len(x)):
                for i in range(grid_points[0].size):
                    y[m][i] = grid_points[m].flatten()[i]
            fit = np.polyfit(x, y, deg=deg)
            p = np.zeros(grid_points[0].size)
            for i in range(p.size):
                p[i] = fit[0][i]*mass*mass + fit[1][i]*mass + fit[2][i]
            return p.reshape(shape)        
        return grid_points[0]


    def run(self):

        self.result['ntoys'] += 1 

        # function calls
        self.arglist[0] = self.n_points
        # convergence
        self.arglist[1] = 10.0
        print("Convergence at EDM %s" % (self.arglist[1]*0.001))
        ##########################################
        self.update_result('time', -time.time())
        
        # make a global fit (DEFAULT)
        massL, massH =  ROOT.Double(0.), ROOT.Double(0.) 
        status = -1
        if self.strategy == 0:
            self.gMinuit.mnexcm( "HES", self.arglist, 1, self.ierflg )
            self.arglist[0] = 1
            self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
            status = self.gMinuit.Migrad()
            #self.arglist[0] = self.n_points
            #self.arglist[1] = 10.0
            #self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )                        

            if status>0:
                self.arglist[0] = 2
                self.gMinuit.mnexcm( "SET STR", self.arglist, 1, self.ierflg )
                status = self.gMinuit.Migrad()

            if self.run_minos:
                #self.arglist[0] = self.n_points
                #self.arglist[1] = 1
                #self.gMinuit.mnexcm( "MINO", self.arglist, 2, self.ierflg )
                self.gMinuit.mnmnot(1, 1, massH, massL)

        # first fit for the cross sections, then fit for the coefficients
        elif self.strategy == 1:
            for key,p in self.map_params.items():            
                if 'A' in key:
                    self.gMinuit.FixParameter(p)
            status = self.gMinuit.Migrad()
            for key,p in self.map_params.items():            
                if 'A' in key:
                    self.gMinuit.Release(p)
            status = self.gMinuit.Migrad()
        else:
            "Strategy not implemented"

        self.update_result('status', status)        
        self.result['time'][-1] += time.time()
        ##########################################

        # Print results
        amin, edm, errdef = ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.)
        nvpar, nparx, icstat = ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985)
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        if self.verbose:
            self.gMinuit.mnprin( 3, amin )

        self.cov = ROOT.TMatrixDSym(self.n_param)
        self.gMinuit.mnemat(self.cov.GetMatrixArray(), self.n_param)

        self.update_result('edm', float(edm))
        self.update_result('amin', float(amin))
        self.update_result('ndof', self.ndof)

        val = ROOT.Double(0.)
        err = ROOT.Double(0.)
        errL = ROOT.Double(0.)
        errH = ROOT.Double(0.)

        for key,p in self.map_params.items():            
            self.gMinuit.GetParameter(p, val, err)                
            true = self.truth[key]
            true_gen = true
            if '_A' not in key and 'mass' not in key:
                true_gen = self.truth[key+'_gen']
                if self.decorrelate:
                    true = self.truth[key+'_prime']

            if '_A' in key and self.do_semianalytic:
                index = (p-1)/6*5 + (p-1)%6 - 1
                val = ROOT.Double(self.theta[index])
                err = ROOT.Double(math.sqrt(self.theta_err[index][index]))

            if key=='mass' and self.run_minos:
                errL = massL-val
                errH = massH-val
            else:
                errL = -err
                errH = err

            self.update_result_array(key=key, vals=(true, true_gen, float(val), float(err), float(errL), float(errH)) )

        if self.do_taylor_expansion:
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                for iload,load in enumerate(self.loads):
                    if iload==0:
                        continue

                    if not self.add_constant_A4:
                        for e in range(self.n_taylor):
                            index_e = (iy*(len(self.loads)-1) + (iload-1))*self.n_taylor + e
                            name = 'coeff'+str(e)+'_y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])+'_A'+str(iload-1)
                            val = ROOT.Double(self.beta[index_e])
                            err = ROOT.Double(math.sqrt(self.beta_err[index_e][index_e]))
                            errL = -err
                            errH = +err
                            self.update_result_array(key=name, vals=(0.0, 0.0, float(val), float(err), float(errL), float(errH)) )
                    else:
                        # A4
                        if iload==5:
                            for e in range(self.n_taylor+1):
                                index_e = iy*( (len(self.loads)-2)*self.n_taylor + 1*(self.n_taylor+1) ) + (len(self.loads)-2)*self.n_taylor + e
                                name = 'coeff'+str(e)+'_y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])+'_A'+str(iload-1)
                                val = ROOT.Double(self.beta[index_e])
                                err = ROOT.Double(math.sqrt(self.beta_err[index_e][index_e]))
                                errL = -err
                                errH = +err
                                self.update_result_array(key=name, vals=(0.0, 0.0, float(val), float(err), float(errL), float(errH)) )
                        # all others
                        else:
                            for e in range(self.n_taylor):
                                index_e = iy*( (len(self.loads)-2)*self.n_taylor + 1*(self.n_taylor+1) ) + (iload-1)*self.n_taylor + e
                                name = 'coeff'+str(e)+'_y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])+'_A'+str(iload-1)
                                val = ROOT.Double(self.beta[index_e])
                                err = ROOT.Double(math.sqrt(self.beta_err[index_e][index_e]))
                                errL = -err
                                errH = +err
                                self.update_result_array(key=name, vals=(0.0, 0.0, float(val), float(err), float(errL), float(errH)) )

        pvalue = -1.0
        if self.do_semianalytic:
            from scipy import stats
            pvalue = 1-stats.chi2.cdf(float(amin), self.ndof)
        self.update_result('pvalue', pvalue)

        if self.do_semianalytic:
            print 'Bins: ', self.data.size, ', number of d.o.f.: ', self.ndof, ' chi2: '+'{:4.1f}'.format(float(amin)), ' (p-value: '+'{:4.3f}'.format(pvalue)+')'
        print('Fit done in '+'{:4.1f}'.format(self.result['time'][-1])+' seconds')

    def update_result(self, var, val):
        if not self.result.has_key(var):
            self.result[var] = [val]
        else:
            self.result[var].append(val)
        print "Unfolder: ", var, " = ", val

    def update_result_array(self, key='', vals=()):
        true      = vals[0]
        true_gen  = vals[1]
        val  = vals[2]
        err  = vals[3]
        errL = vals[4]
        errH = vals[5]
        if self.verbose:
            print key, ":", true, " ==> ", val, "+/-", err, ' (', errL, ',', errH, ')'            
        if not self.result.has_key(key):
            self.result[key] = {'true' : [true], 
                                'toy' : [true_gen], 
                                'fit' : [val], 
                                'err' : [err], 
                                'errL' :  [errL],  
                                'errH' :  [errH] }
        else:
            self.result[key]['true'].append(true) 
            self.result[key]['toy'].append(true_gen) 
            self.result[key]['fit'].append(val) 
            self.result[key]['err'].append(err) 
            self.result[key]['errL'].append(errL) 
            self.result[key]['errH'].append(errH) 
            


    # save the result to a pickle file
    def save_result(self):            
        pickle.dump(self.result, self.f)
        self.f.close()
        print "Result saved in file", self.f.name
        return
