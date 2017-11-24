import ROOT
import math
import sys
import copy
import numpy as np
import numpy.random as ran

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from template_parameters import params_test, pdf_test 

class MixData:

    def __init__(self, input_dir='../data/', output_dir='../data/'):
        print "Initialize MixData"
        self.input_dir = input_dir
        self.output_dir = output_dir
        return

    def mix(self, mass=80.000, pt_bin=[], y_bin=[], coeff=[], pdf=(lambda x : 1.0), params={}):

        pt_edges = np.array([])
        y_edges = np.array([])        
        mix = np.array([])
        for ipt,pt in enumerate(params['params_W']['pt']):
            if pt<pt_bin[0] or pt>=pt_bin[1]:
                continue
            print "pt %s is in range" % pt 
            for ic0,c0 in enumerate(params['params_W']['A0']):
                if c0!=coeff[0]:
                    continue
                for ic1,c1 in enumerate(params['params_W']['A1']):
                    if c1!=coeff[1]:
                        continue
                    for ic2,c2 in enumerate(params['params_W']['A2']):
                        if c2!=coeff[2]:
                            continue
                        for ic3,c3 in enumerate(params['params_W']['A3']):
                            if c3!=coeff[3]:
                                continue
                            for ic4,c4 in enumerate(params['params_W']['A4']):
                                if c4!=coeff[4]:
                                    continue

                                in_name = 'pt{:02.1f}'.format(pt)+'_'+'y{:02.1f}'.format(0.00)+'_M{:05.3f}'.format(mass)

                                coeff = [c0,c1,c2,c3,c4]
                                for c in range(5):
                                    in_name += ('_A'+str(c)+('{:03.2f}'.format(coeff[c])))

                                # load the file
                                sample = np.load(self.input_dir+'/grid_lab_'+in_name+'.npy')

                                y_edges = sample[1]
                                pt_edges = sample[2]

                                # must be odd!
                                #middle_point = (len(params['params_W']['y'])-1)/2
                                middle_point = (len( y_edges )-1)/2
                                y_bin_width = (y_edges[1]-y_edges[0])

                                shift_pos = range( int((y_bin[0]-y_edges[0])/y_bin_width), 
                                                   int((y_bin[1]-y_edges[0])/y_bin_width) + 1)
                                        
                                #for iy in [  middle_point-80,middle_point, middle_point+130]:
                                #for iy,y in np.ndenumerate( y_edges ):            

                                for iy in shift_pos:
                                #for iy,y in np.ndenumerate(params['params_W']['y']):            
                                    y = y_edges[0]+y_bin_width*iy
                                    print "y %s is in range" % y 

                                    weight = pdf(pt,y, c0, c1, c2, c3, c4)
                                    #weight=1.0

                                    shift = iy-middle_point
                                    #shift = int(y/y_bin_width)

                                    sam = copy.deepcopy(sample[0])
                                    if shift>0:
                                        sam[-shift:, :] = 0.
                                    elif shift<0:
                                        sam[:-shift, :] = 0.
                                    sam_shifted = np.roll(sam, shift, axis=0 ) 
                                    if not mix.shape==sam.shape:
                                        mix = sam_shifted*weight
                                    else:
                                        mix += sam_shifted*weight

        #print(mix)
        out_name = 'pt{:02.1f}'.format(pt_bin[0])+'-'+'{:02.1f}'.format(pt_bin[1])+'_'+'y{:02.1f}'.format(y_bin[0])+'-'+'y{:02.1f}'.format(y_bin[1])+'_M{:05.3f}'.format(mass)
        for c in range(5):
            out_name += ('_A'+str(c)+('{:03.2f}'.format(coeff[c])))
        np.save(self.output_dir+'/mixed_dataset_'+out_name, mix)

        xx, yy = np.meshgrid(pt_edges, y_edges)        
        plt.pcolormesh(yy, xx, mix)
        plt.show()
        plt.savefig(self.output_dir+'/mixed_dataset_'+out_name+'.png')

##################################
