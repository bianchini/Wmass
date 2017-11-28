import ROOT
import math
import numpy as np
from array import array;


class Unfolder:

    def __init__(self, input_dir='../data/', params={} ):

        self.input_dir=input_dir
        self.load_files(params=params)

        self.n_param = 10

        p_start = array( 'd' )
        p_step = array( 'd' )
        #to_be_fixed = array('i')

        for i in range(self.n_param):  
            p_start.append( 0.0 )
            p_step.append( 0.01 )


        self.gMinuit = ROOT.TMinuit(50)  
        self.gMinuit.SetFCN( self.fcn ) 

        self.arglist = array( 'd', 10*[0.] )
        self.ierflg = ROOT.Long(0)
        self.arglist[0] = 0.5
        self.gMinuit.mnexcm( "SET ERR", self.arglist, 1, self.ierflg )

        for p in range(self.n_param):
            p_in = p_start[p]
            delta_p = p_step[p] 
            p_low = -1.0
            p_high = +1.0
            self.gMinuit.mnparm( p, "par_"+str(p), p_in, delta_p, p_low,  p_high, self.ierflg )
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

        for ipt in range(len(self.input_shapes_pt)-1):
            pt_bin=[ self.input_shapes_pt[ipt], self.input_shapes_pt[ipt+1] ]
            for iy in range(len(self.input_shapes_y)-1):
                y_bin=[ self.input_shapes_y[iy], self.input_shapes_y[iy+1] ]
                for im,m in enumerate(self.input_shapes_mass):
                    for iA0,A0 in enumerate(self.input_shapes_A0):
                        for iA1,A1 in enumerate(self.input_shapes_A1):
                            for iA2,A2 in enumerate(self.input_shapes_A2):
                                for iA3,A3 in enumerate(self.input_shapes_A3):
                                    for iA4,A4 in enumerate(self.input_shapes_A4):
                                        input_name = 'mixed_dataset_'+'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])+'_M'+'{:05.3f}'.format(m)
                                        for ic,c in enumerate([A0,A1,A2,A3,A4]):
                                            input_name += '_A'+str(ic)+('{:03.2f}'.format(c))
                                        print "Loading file...", input_name 
                                        setattr(self, input_name, np.load(self.input_dir+'/'+input_name+'.npy') )


    def run(self, n_points=50000):

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
        

    def fcn(self, npar, gin, f, par, iflag ):
        return 1.0
