import ROOT
import math
import sys
from pprint import pprint
import numpy as np
from array import array;

from template_parameters import accept_point, pdf_test

class Unfolder:

    def __init__(self, input_dir='../data/', params={} ):

        self.input_dir=input_dir

        self.load_files(params=params)
        self.read_data(num_events=10000000)
        self.book_parameters(params=params)

        #p_start = array( 'd' )
        #p_step = array( 'd' )
        #to_be_fixed = array('i')
        #for i in range(self.n_param):  
        #    p_start.append( 0.0 )
        #    p_step.append( 0.01 )

            #if p in to_be_fixed:
            #    self.gMinuit.FixParameter(p)

    def load_files(self, params={}):

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
                                        grid /= (grid.sum() if grid.sum()>0.0 else 1.0)
                                        setattr(self, input_name, grid)
                                        self.shape = grid.shape                                         

    def read_data(self, num_events=100000):

        self.truth = {}

        total_weight = 0.0
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                total_weight += pdf_test(pt=pt_bin[0], y=y_bin[0])

        self.data = np.zeros(self.shape)
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]

                weight = pdf_test(pt=pt_bin[0], y=y_bin[0])/total_weight
                input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(80.000)
                for ic,c in enumerate([0.0, 0.0, 0.0, 0.0, 0.0]):
                    input_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                data_rnd = np.random.poisson( getattr(self, input_name)*(num_events*weight))
                print getattr(self, input_name).sum(), num_events, weight
                self.truth['pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])] = data_rnd.sum()
                self.data += data_rnd

        self.num_events = self.data.sum()
        pprint(self.truth)
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
        self.gMinuit.mnparm( self.n_param, "mass", 80.000, 0.100,  79.000, 81.000, self.ierflg )        
        self.map_params['mass'] = self.n_param
        #self.gMinuit.FixParameter(0)
        self.n_param += 1

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                par_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])
                self.map_params[par_name] = self.n_param
                self.gMinuit.mnparm( self.n_param, par_name, self.truth[par_name], 10.0, 0.0, self.num_events, self.ierflg )
                #self.gMinuit.FixParameter(self.n_param)
                self.n_param += 1

                for coeff in ['A0', 'A4']:
                # add one parameter per pt/y/Ai bin                    
                    par_name_coeff = par_name+'_'+coeff 
                    self.map_params[par_name_coeff] = self.n_param
                    self.gMinuit.mnparm( self.n_param, par_name_coeff, 0.00, 0.01, -1.00,  1.00, self.ierflg )
                    #if coeff in ['A1', 'A2', 'A3']:
                    self.gMinuit.FixParameter(self.n_param)
                    self.n_param += 1            

        self.loads = [ [0.0, 0.0, 0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0, 0.0, 0.0],
                       [0.0, 0.0, 0.0, 0.0, 1.0]
                       ]


                                        
    def fcn(self, npar, gin, f, par, iflag ):

        #print "Evaluating fcn..."

        pdf = np.zeros(self.shape)

        # find closest teplates for interpolation
        mass = par[0]
        closest=[]
        for im in range(len(self.input_shapes_mass)-1):            
            mass_low = self.input_shapes_mass[im]
            mass_high = self.input_shapes_mass[im+1]
            if mass>=mass_low and mass<mass_high:
                closest.extend( [im,im+1, (mass-mass_low)/(mass_high-mass_low)] )
        #print mass,closest

        count_param = 1
        pdf_bin = np.zeros(self.shape)
        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]                
                norm_bin = par[count_param]
                #print "norm_bin", norm_bin

                # pdf at a given (pt,y)
                # pdf_tmp1 = Norm*{(1-w)*[(1-A0-A2)*pdf(M,  pt,y,UL) + A0*pdf(M,  pt,y,L) + A1*pdf(M,  pt,y,P)] + 
                #                      w*[(1-A0-A2)*pdf(M+1,pt,y,UL) + A0*pdf(M+1,pt,y,L) + A1*pdf(M+1,pt,y,P)] } 
                pdf_tmp1 = np.zeros(self.shape)

                load_name = 'mixed_dataset_pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])

                #print "par[count_param+1]", par[count_param+1]
                #print "par[count_param+2]", par[count_param+2]
                # mass interpolation
                for im in range(2):
                    # (1-A0-A2)*pdf(M,  pt,y,UL) + A0*pdf(M,  pt,y,L) + A1*pdf(M,  pt,y,P)]                 
                    pdf_tmp2 = np.zeros(self.shape)
                    mass_weight = 1-closest[2] if im==0 else closest[2]
                    for iload, load in enumerate(self.loads):
                        #print "mass_weight", mass_weight
                        grid_name = load_name + '_M'+'{:05.3f}'.format(self.input_shapes_mass[closest[im]]) 
                        for ic,c in enumerate( load ):
                            grid_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                        if iload==0:
                            pdf_tmp2 += (1-par[count_param+1]-par[count_param+2])*getattr(self, grid_name)
                        elif iload==1:
                            pdf_tmp2 += par[count_param+1]*getattr(self, grid_name)
                        elif iload==2:
                            pdf_tmp2 += par[count_param+2]*getattr(self, grid_name)
                    pdf_tmp2 *= mass_weight
                    pdf_tmp1 += pdf_tmp2
                
                pdf_tmp1 *= norm_bin
                pdf += pdf_tmp1
                count_param += 3


        nll = 0.0
        #for ipt in range(self.data.shape[1]):
        #    for iy in range( self.data.shape[0]):
        #        mu = pdf[iy][ipt]
        #        #print "(%s,%s) ==> %s" % (iy,ipt, mu)
        #        if mu<=0.0:
        #            continue
        #        n = self.data[iy][ipt]
        #        nll += (-math.log(mu)*n + mu)
                
        nll = ((-self.data*np.log(pdf) + pdf) if np.isfinite(np.log(pdf)).all() else np.ones(self.shape)*sys.float_info.max).sum()
        #print nll
        #nll = (mass-80.000)*(mass-80.000)/0.500/0.500
        f[0] = nll
        return


    def run(self, n_points=500000):

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
        #for p in range(self.n_param):
        #    self.gMinuit.GetParameter(p, val, err)
        #    print p, val

        #pprint(self.map_params)
        for key,p in self.map_params.items():            
            if 'mass' in key or '_A' in key:
                continue
            self.gMinuit.GetParameter(p, val, err)
            true = self.truth[key]
            print key, p, ":", true, " ==> ", val, "+/-", err
        #pprint(self.truth)
