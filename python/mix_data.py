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

    def mix(self, mass=80.000, pdf=(lambda x : 1.0), params={}):

        pt_edges = np.array([])
        y_edges = np.array([])        
        mix = np.array([])
        for ipt,pt in enumerate(params['params_W']['pt']):
            #for iy,y in enumerate(params['params_W']['y']):
            for ic0,c0 in enumerate(params['params_W']['A0']):
                for ic1,c1 in enumerate(params['params_W']['A1']):
                    for ic2,c2 in enumerate(params['params_W']['A2']):
                        for ic3,c3 in enumerate(params['params_W']['A3']):
                            for ic4,c4 in enumerate(params['params_W']['A4']):
                                out_name = 'pt{:02.1f}'.format(pt)+'_'+'y{:02.1f}'.format(0.00)+'_M{:05.3f}'.format(mass)
                                coeff = [c0,c1,c2,c3,c4]
                                for c in range(5):
                                    out_name += ('_A'+str(c)+('{:03.2f}'.format(coeff[c])))
                                sample = np.load(self.input_dir+'/grid_lab_'+out_name+'.npy')
                                #print sample

                                sam = copy.deepcopy(sample[0])
                                for iy,y in enumerate(sample[1]):                                    
                                    weight = pdf(pt,y, c0, c1, c2, c3, c4)
                                    #print weight

                                    if y>=0.:
                                        sam[:,-1] = 0.
                                        print iy-len(sample[1])/2 
                                        np.roll(sam, iy-len(sample[1])/2  )
                                    else:
                                        sam[:,0] = 0.
                                        np.roll(sam, iy+len(sample[1])/2 )                                        

                                    if not mix.shape==sam.shape:
                                        mix = sam*weight
                                        pt_edges = sample[2]
                                        y_edges = sample[1]
                                    else:
                                        mix += sam*weight

        #print(mix)

        xx, yy = np.meshgrid(pt_edges, y_edges)        
        plt.pcolormesh(yy, xx, mix)
        plt.show()
        plt.savefig(self.output_dir+'/mixed_dataset.png')


##################################
