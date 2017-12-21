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

from template_parameters import accept_point, pdf_test

class Unfolder:

    def __init__(self, input_dir='../data/', params={},  mass=80.0000, num_events=1000000, fix=[], interp_deg=1, n_points=500000, job_name='TEST', verbose=True, prior_coeff=0.3, prior_xsec=0.3, strategy=0 ):

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

        self.load_files(params=params)
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

    def load_files(self, params={}):

        self.input_pt_bins = np.linspace(params['params_lep']['pt_range'][0], params['params_lep']['pt_range'][1],  params['params_lep']['pt_bins']+1)
        self.input_y_bins  = np.linspace(params['params_lep']['y_range'][0], params['params_lep']['y_range'][1],  params['params_lep']['y_bins']+1)

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
                                        if self.verbose:
                                            print "Loading file...", input_name 
                                        grid = np.load(self.input_dir+'/'+input_name+'.npy')
                                        norm = grid.sum() 
                                        #grid /= (norm if norm>0.0 else 1.0)
                                        setattr(self, input_name, grid)
                                        setattr(self, input_name+'_norm', norm)
                                        self.shape = grid.shape                                         

    # generate a rnd sample
    def load_data(self):
        self.toy_data()

    # generate rnd samples
    def toy_data(self, ntoy=0):

        self.ntoy = ntoy
        self.data = np.zeros(self.shape)
        self.truth = {}

        normalisation = 0.0
        #total_weight = 0.0
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                #total_weight += pdf_test(pt=pt_bin[0], y=y_bin[0])
                input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate([0.0, 0.0, 0.0, 0.0, 0.0]):
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
                for ic,c in enumerate([0.0, 0.0, 0.0, 0.0, 0.0]):
                    input_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                #data_rnd = np.random.poisson( getattr(self, input_name)*getattr(self, input_name+'_norm')*self.num_events/norm_acceptance )
                n_true = self.num_events*weight/normalisation
                data_rnd = np.random.poisson( getattr(self, input_name)*n_true ) #* (max(np.random.normal(1.0, self.prior_xsec), 0.0))
                self.truth[name] = n_true
                self.truth[name+'_gen'] = data_rnd.sum()/getattr(self, input_name+'_norm')
                self.truth[name+'_A0'] = 0.0
                self.truth[name+'_A1'] = 0.0
                self.truth[name+'_A2'] = 0.0
                self.truth[name+'_A3'] = 0.0
                self.truth[name+'_A4'] = 0.0
                self.data += data_rnd

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
        #self.gMinuit.mnparm( self.n_param, "mass", self.mass, 0.010, self.mass-1.000, self.mass+1.000, self.ierflg )        
        self.gMinuit.DefineParameter( self.n_param, "mass", self.mass, 0.010, self.mass-1.000, self.mass+1.000 )        
        self.map_params['mass'] = self.n_param
        if 'mass' in self.fix:
            self.gMinuit.FixParameter(0)
        self.n_param += 1

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                par_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'{:03.2f}'.format(y_bin[1])
                self.map_params[par_name] = self.n_param
                # 1/3 -- 3 with 10 step
                self.gMinuit.DefineParameter( self.n_param, par_name, self.truth[par_name], 10., self.truth[par_name]/3.0, self.truth[par_name]*3.0  )
                if 'pt_y' in self.fix: 
                    self.gMinuit.FixParameter(self.n_param)
                self.n_param += 1

                for coeff in ['A0', 'A1', 'A2', 'A3', 'A4']:
                # add one parameter per pt/y/Ai bin                    
                    par_name_coeff = par_name+'_'+coeff 
                    self.map_params[par_name_coeff] = self.n_param            
                    self.gMinuit.DefineParameter( self.n_param, par_name_coeff, 0.00, 0.001, -0.75,  0.75)
                    if coeff in self.fix:
                        self.gMinuit.FixParameter(self.n_param)
                    self.n_param += 1            

        self.loads = [ [0.0, 0.0, 0.0, 0.0, 0.0],
                       [2.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0, 0.0, 0.0],
                       [0.0, 0.0, 1.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0, 1.0, 0.0],
                       [0.0, 0.0, 0.0, 0.0, 2.0]
                       ]


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
                self.arglist[1] = self.truth[par_name] 
                self.gMinuit.mnexcm("SET PAR", self.arglist, 2, self.ierflg )

                for coeff in ['A0', 'A1', 'A2', 'A3', 'A4']:
                    par_name_coeff = par_name+'_'+coeff 
                    self.arglist[0] = self.map_params[par_name_coeff]+1
                    self.arglist[1] = self.truth[par_name_coeff] 
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

                                        
    # pdf at a given (pt,y):
    # pdf_tmp_pt_y = Norm*{(1-w)*[(1-A0-A2)*pdf(M,  pt,y,UL) + A0*pdf(M,  pt,y,L) + A1*pdf(M,  pt,y,P)] + 
    #                      w*[(1-A0-A2)*pdf(M+1,pt,y,UL) + A0*pdf(M+1,pt,y,L) + A1*pdf(M+1,pt,y,P)] } 
    def fcn(self, npar, gin, f, par, iflag ):

        #print "Evaluating fcn..."
        pdf = np.zeros(self.shape)

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

                    # normalise 
                    #pdf_tmp_coeff /= pdf_tmp_coeff.sum()
                    grid_points.append( pdf_tmp_coeff )

                # interpolate templates of different mass                
                pdf_tmp_pt_y += self.interpolate_mass(mass=mass, x=x, grid_points=grid_points, deg=self.interp_deg, shape=pdf.shape)
                # normalise bin using norm_bin
                pdf_tmp_pt_y *= norm_bin
                # add tis bin to the total pdf
                pdf += pdf_tmp_pt_y
                # increase par counter by 3 units: (pt_y, A0, A1, ..., A4)
                count_param += (1+len(self.loads[0]))

        # compute -2*log likelihood
        nll = 0.0

        # prior
        for key,p in self.map_params.items():            
            # pt_y bins
            if '_A' not in key and 'mass' not in key:                
                #true = self.truth[key+'_gen']
                true = self.truth[key]
                err = self.prior_xsec
                nll += math.pow((par[p]-true)/(err*true), 2.0)
            # coefficients
            elif '_A' in key:
                true = self.truth[key]
                err = self.prior_coeff
                nll += math.pow((par[p]-true)/err, 2.0)
                
        # prevent from any zeros in the pdf
        pdf += (np.ones(self.shape)*sys.float_info.min)
        nll += 2*(-self.data*np.log(pdf) + pdf).sum()

        #print("{:16f}".format(nll))
        #for key,p in self.map_params.items():          
        #    if '_A' not in key and 'mass' not in key:                
        #        print key, par[p]

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
        self.arglist[1] = 1.
        print("Convergence at EDM %s" % (self.arglist[1]*0.001))
        ##########################################
        self.update_result('time', -time.time())
        
        # make a global fit (DEFAULT)
        status = -1
        if self.strategy == 0:
            #self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )
            status = self.gMinuit.Migrad()

        # first fit for the cross sections, then fit for the coefficients
        elif strategy == 1:
            for key,p in self.map_params.items():            
                if 'A' in key:
                    self.gMinuit.FixParameter(p)
            self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )
            for key,p in self.map_params.items():            
                if 'A' in key:
                    self.gMinuit.Release(p)
            self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )
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

        val = ROOT.Double(0.)
        err = ROOT.Double(0.)
        for key,p in self.map_params.items():            
            self.gMinuit.GetParameter(p, val, err)
            true = self.truth[key]
            true_gen = true
            if '_A' not in key and 'mass' not in key:
                true_gen = self.truth[key+'_gen']
            if self.verbose:
                print key, p, ":", true, " ==> ", val, "+/-", err            
            if not self.result.has_key(key):
                self.result[key] = {'true' : [true], 'toy' : [true_gen], 'fit' : [float(val)], 'err' : [float(err)] }
            else:
                self.result[key]['true'].append(true) 
                self.result[key]['toy'].append(true_gen) 
                self.result[key]['fit'].append(float(val)) 
                self.result[key]['err'].append(float(err)) 
        print('Fit done in '+'{:4.1f}'.format(self.result['time'][-1])+' seconds')

    def update_result(self, var, val):
        if not self.result.has_key(var):
            self.result[var] = [val]
        else:
            self.result[var].append(val)
        print "Unfolder: ", var, " = ", val

    # save the result to a pickle file
    def save_result(self):            
        pickle.dump(self.result, self.f)
        self.f.close()
        print "Result saved in file", self.f.name
        return
