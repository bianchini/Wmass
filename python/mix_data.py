import ROOT
import math
import sys
import copy
import numpy as np
import numpy.random as ran

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from template_parameters import params_test, pdf_test , coefficients_test, accept_point

class MixData:

    def __init__(self, input_dir='../data/', output_dir='../data/'):
        print "Initialize MixData"
        self.input_dir = input_dir
        self.output_dir = output_dir
        return

    def add_shapes(self, params, symmetrise=False, make_templates=False):

        # these are the available shapes (all generated for y=0.00)
        self.input_shapes_pt = params['params_W']['pt']

        # these are the shapes to be merged
        self.output_shapes_pt = params['params_template']['pt']
        self.output_shapes_y  = params['params_template']['y']

        # these are the parameters to be varied when making templates
        self.output_shapes_A0  = params['params_template']['A0']
        self.output_shapes_A1  = params['params_template']['A1']
        self.output_shapes_A2  = params['params_template']['A2']
        self.output_shapes_A3  = params['params_template']['A3']
        self.output_shapes_A4  = params['params_template']['A4']
        self.output_shapes_mass  = params['params_template']['mass']

        # this is the function that returns the coefficients
        self.make_templates = make_templates
        self.symmetrise = symmetrise


    def mix_bin(self, y_bin=[], pt_bin=[], pdf=(lambda x : 1.0), coefficients=coefficients_test):

        # the merged sample
        mix = np.array([])

        pt_edges = np.array([])
        y_edges = np.array([])        
        
        # load a test file to read the bon edges
        test_name = 'pt{:02.1f}'.format(0.0)+'_'+'y{:03.2f}'.format(0.00)+'_M{:05.3f}'.format(80.000)
        for c in range(5):
            p =  0.0
            test_name += ('_A'+str(c)+('{:03.2f}'.format(p)))
        grid_test = np.load(self.input_dir+'/grid_lab_'+test_name+'.npy')

        # y and pt bins of the input shape
        y_edges = grid_test[1]
        pt_edges = grid_test[2]

        # must be odd!
        middle_point = (len( y_edges ) - 1)/2
        y_bin_width = (y_edges[1]-y_edges[0])
        shift_pos = range( int((y_bin[0]-y_edges[0])/y_bin_width), 
                           int((y_bin[1]-y_edges[0])/y_bin_width) )        

        # central bin (for symmetrise)        
        central = np.array([])

        for ipt,pt in enumerate(self.input_shapes_pt):
            
            if pt<pt_bin[0] or pt>=pt_bin[1]:
                continue

            for iy in shift_pos:
                y = y_edges[0]+y_bin_width*iy
                print "\t (pt, y) = (%s, %s) is in range" % (pt, y)

                in_name = 'pt{:02.1f}'.format(pt)+'_'+'y{:03.2f}'.format(0.00)
                if not self.make_templates:
                    in_name += '_M{:05.3f}'.format( self.output_shapes_mass[0] )
            
                for ic,c in enumerate( coefficients(pt_bin[0],y_bin[0]) ):
                    c = c if abs(c)>0. else 0.0
                    if self.make_templates and ic==0:                         
                        in_name      += ('_M'+('{:05.3f}'.format(c)))
                    # these are the coefficients (if make_templates)
                    elif self.make_templates and ic>0:
                        in_name      += ('_A'+str(ic-1)+('{:03.2f}'.format(c)))
                    # these are the coefficients as a functon of (pt,y) (if NOT make_templates)
                    else:
                        in_name      += ('_A'+str(ic)+('{:03.2f}'.format(c)))

                # load the file
                grid = np.load(self.input_dir+'/grid_lab_'+in_name+'.npy')

                # mix sub-samples according to a prior pdf
                weight = pdf(pt=pt,y=y)
                
                shifts = [ [iy-middle_point, grid] ]                
                for [shift,sample] in shifts:
                    sam = copy.deepcopy(sample[0]) 
                    if self.symmetrise and shift==0:
                        if central.shape!=sam.shape:
                            central = sam*weight
                        else:
                            central += sam*weight
                        continue
                    if shift>0:
                        sam[-shift:, :] = 0.
                    elif shift<0:
                        sam[:-shift, :] = 0.
                    sam_shifted = np.roll(sam, shift, axis=0 ) 
                    if not mix.shape==sam.shape:
                        mix = sam_shifted*weight
                    else:
                        mix += sam_shifted*weight

        if self.symmetrise:
            for ipt in range( mix.shape[1] ):
                for iy in range( mix.shape[0]/2 ):
                    #print iy, "+", mix.shape[0]-1-iy
                    # symmetrise varitions
                    left = mix[iy][ipt]
                    right = mix[mix.shape[0]-1-iy][ipt] 
                    mix[iy][ipt] += right
                    mix[mix.shape[0]-1-iy][ipt] += left
                    # symmetrise central bin
                    if central.shape==mix.shape:
                        leftC = central[iy][ipt]
                        rightC = central[mix.shape[0]-1-iy][ipt] 
                        central[iy][ipt] += rightC
                        central[mix.shape[0]-1-iy][ipt] += leftC
            if central.shape==mix.shape:
                mix += central*0.5
            

        out_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:03.2f}'.format(y_bin[0])+'-'+'y{:03.2f}'.format(y_bin[1])
        if not self.make_templates:
            out_name += '_M{:05.3f}'.format( self.output_shapes_mass[0] )
        else:
            for ic,c in enumerate( coefficients(0.0, 0.0) ):
                c = c if abs(c)>0. else 0.0
                if ic==0:
                    out_name += ('_M'+('{:05.3f}'.format(c)))
                else:
                    out_name += ('_A'+str(ic-1)+('{:03.2f}'.format(c)))

        # normalisation
        #mix /= (mix.sum() if mix.sum()>0.0 else 1.0)
        np.save(self.output_dir+'/mixed_dataset_'+out_name, mix)

        xx, yy = np.meshgrid(pt_edges, y_edges)        
        plt.pcolormesh(yy, xx, mix)
        plt.show()
        plt.savefig(self.output_dir+'/mixed_dataset_'+out_name+'.png')

    def mix_all_bins(self):
        for ipt in range(len(self.output_shapes_pt)-1):
            pt_bin=[ self.output_shapes_pt[ipt], self.output_shapes_pt[ipt+1] ]
            for iy in range(len(self.output_shapes_y)-1):
                y_bin=[ self.output_shapes_y[iy], self.output_shapes_y[iy+1] ]
                print "Call MixData for (%s,%s)" % (pt_bin, y_bin)
                if not self.make_templates:
                    self.mix_bin(pt_bin=pt_bin, y_bin=y_bin, pdf=pdf_test)
                else:
                    for im,m in enumerate(self.output_shapes_mass):
                        for iA0,A0 in enumerate(self.output_shapes_A0):
                            for iA1,A1 in enumerate(self.output_shapes_A1):
                                for iA2,A2 in enumerate(self.output_shapes_A2):
                                    for iA3,A3 in enumerate(self.output_shapes_A3):
                                        for iA4,A4 in enumerate(self.output_shapes_A4):
                                            if not accept_point(coeff=[A0,A1,A2,A3,A4]):
                                                continue
                                            coefficients = lambda x,y : [m,A0,A1,A2,A3,A4]
                                            print "Mixing template for", coefficients(0.0, 0.0)
                                            self.mix_bin(pt_bin=pt_bin, y_bin=y_bin, pdf=pdf_test, coefficients=coefficients) 

 


##################################
