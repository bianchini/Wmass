import ROOT
import math
import sys
import csv
from pprint import pprint
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from array import array;

from template_parameters import accept_point, pdf_test

class Unfolder:

    def __init__(self, input_dir='../data/', params={},  mass=80.0000, num_events=1000000, fix=[], interp_deg=1 ):

        self.input_dir=input_dir
        self.fix=fix
        self.mass = mass
        self.interp_deg=interp_deg
        self.num_events=num_events
        self.load_files(params=params)
        self.load_data()
        self.book_parameters(params=params)

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
                                        input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(m)
                                        for ic,c in enumerate([A0,A1,A2,A3,A4]):
                                            input_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                                        print "Loading file...", input_name 
                                        grid = np.load(self.input_dir+'/'+input_name+'.npy')
                                        norm = grid.sum() 
                                        grid /= (norm if norm>0.0 else 1.0)
                                        setattr(self, input_name, grid)
                                        setattr(self, input_name+'_norm', norm)
                                        self.shape = grid.shape                                         

    def load_data(self):
        self.toy_data()

    def toy_data(self, ntoy=0):

        self.ntoy=ntoy
        self.data = np.zeros(self.shape)
        self.truth = {}

        norm_acceptance = 0.0
        #total_weight = 0.0
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                #total_weight += pdf_test(pt=pt_bin[0], y=y_bin[0])
                input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate([0.0, 0.0, 0.0, 0.0, 0.0]):
                    input_name += '_A'+str(ic)+('{:03.2f}'.format(c))                    
                norm_acceptance += getattr(self, input_name+'_norm')
        
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                #weight = pdf_test(pt=pt_bin[0], y=y_bin[0])/total_weight
                #print pt_bin, y_bin, weight
                name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])
                input_name = 'mixed_dataset_'+name+'_M'+'{:05.3f}'.format(self.mass)
                for ic,c in enumerate([0.0, 0.0, 0.0, 0.0, 0.0]):
                    input_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                data_rnd = np.random.poisson( getattr(self, input_name)*getattr(self, input_name+'_norm')*self.num_events/norm_acceptance )
                self.truth[name] = data_rnd.sum()
                self.truth[name+'_A0'] = 0.0
                self.truth[name+'_A4'] = 0.0
                self.data += data_rnd

        self.truth['mass'] = self.mass
        self.truth['num_events_toy'] = self.num_events
        self.truth['num_events_toy'] = self.data.sum()
        #pprint(self.truth)
        xx, yy = np.meshgrid(self.input_pt_bins, self.input_y_bins)        
        plt.pcolormesh(yy, xx, self.data)
        plt.show()
        plt.savefig('data_toy_'+str(ntoy)+'.png')
        #print self.data

        
    def book_parameters(self, params={}):

        self.gMinuit = ROOT.TMinuit(500)  
        self.gMinuit.SetFCN( self.fcn ) 
        self.arglist = array( 'd', 10*[0.] )
        self.ierflg = ROOT.Long(0)
        self.arglist[0] = 0.5
        self.gMinuit.mnexcm( "SET ERR", self.arglist, 1, self.ierflg )

        self.map_params={}

        # add mass
        self.n_param = 0
        self.gMinuit.mnparm( self.n_param, "mass", self.mass, 0.010, self.mass-1.000, self.mass+1.000, self.ierflg )        
        self.map_params['mass'] = self.n_param
        if 'mass' in self.fix:
            self.gMinuit.FixParameter(0)
        self.n_param += 1

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                par_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])
                self.map_params[par_name] = self.n_param
                self.gMinuit.mnparm( self.n_param, par_name, self.truth[par_name], 10.0, self.truth[par_name]/10., self.truth[par_name]*10, self.ierflg )
                if 'pt_y' in self.fix: 
                    self.gMinuit.FixParameter(self.n_param)
                self.n_param += 1

                for coeff in ['A0', 'A4']:
                # add one parameter per pt/y/Ai bin                    
                    par_name_coeff = par_name+'_'+coeff 
                    self.map_params[par_name_coeff] = self.n_param
                    self.gMinuit.mnparm( self.n_param, par_name_coeff, 0.00, 0.01, -1.50,  1.50, self.ierflg )
                    if coeff in self.fix:
                        self.gMinuit.FixParameter(self.n_param)
                    self.n_param += 1            

        self.loads = [ [0.0, 0.0, 0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0, 0.0, 1.0]
                       ]

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

        # find closest teplates for linear interpolation
        mass = par[0]
        mass_neighbours = self.find_neighbours(mass=mass)

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

                load_name = 'mixed_dataset_pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])

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
                        if iload==0:
                            pdf_tmp_coeff += (1-par[count_param+1]-par[count_param+2])*grid
                        elif iload==1:
                            pdf_tmp_coeff += par[count_param+1]*grid
                        elif iload==2:
                            pdf_tmp_coeff += par[count_param+2]*grid                    
                    grid_points.append( pdf_tmp_coeff )

                # interpolate templates of different mass
                pdf_tmp_pt_y += self.interpolate_mass(mass=mass, x=x, grid_points=grid_points, deg=self.interp_deg, shape=pdf.shape)
                # normalise bin using norm_bin
                pdf_tmp_pt_y *= norm_bin
                # add tis bin to the total pdf
                pdf += pdf_tmp_pt_y
                # increase par counter by 3 units: (pt_y, A0, A4)
                count_param += 3

        # compute the negative log likelihood
        nll = 0.0

        # prior
        for key,p in self.map_params.items():            
            true = self.truth[key]
            # pt_y bins
            if '_A' not in key and 'mass' not in key:                
                err = 0.30
                nll += math.pow( (math.log(par[p])-math.log(true))/err, 2.0)
            # coefficients
            elif '_A' in key:
                err = 0.30
                nll += math.pow((par[p]-true)/err, 2.0)
                
        nll += ((-self.data*np.log(pdf) + pdf) if np.isfinite(np.log(pdf)).all() else np.ones(self.shape)*sys.float_info.max).sum()

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


    def run(self, n_points=1000000):

        result = {}
        result['n_points'] = n_points
        result['mass_gen'] = self.mass
        result['fix'] = self.fix

        self.arglist[0] = n_points
        self.arglist[1] = 1.

        self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )

        # Print results
        amin, edm, errdef = ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.)
        nvpar, nparx, icstat = ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985)
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        self.gMinuit.mnprin( 3, amin )

        self.cov = ROOT.TMatrixDSym(self.n_param)
        self.gMinuit.mnemat(self.cov.GetMatrixArray(), self.n_param)

        val = ROOT.Double(0.)
        err = ROOT.Double(0.)
        for key,p in self.map_params.items():            
            self.gMinuit.GetParameter(p, val, err)
            true = self.truth[key]
            print key, p, ":", true, " ==> ", val, "+/-", err
            result[key]={ 'true' : true, 'fit' : val, 'fit_err' : err}
            
        f = open('mycsvfile'+str(self.ntoy)+'.txt','wb')
        #w = csv.DictWriter(f, result.keys())
        #w.writerows(result)
        f.write(result)
        f.close()
