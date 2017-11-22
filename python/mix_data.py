import ROOT
import math
import sys
import numpy as np
import numpy.random as ran

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


class MixData():

    def __init__(self):
        print "Initialize MixData"
        self.input_dir='../test/'
        return

    def mix(self,  coeff=[0.,0.,0.,0.,0.], mass=80.000, distribution=np.array([]) ):

        grid_tag = 'grid_lab_M{:05.3f}'.format(mass)
        for ic,c in enumerate(coeff):
            if ic<4:
                grid_tag += '_A'+str(ic)+('{:03.2f}'.format(c))

        proportions = distribution[0]
        proportions*(1./np.sum(proportions))
        pt_vals = distribution[1]
        y_vals = distribution[2]
        A4_vals = distribution[3]
        
        mix = np.array([])
        pt_edges = np.array([])
        y_edges = np.array([])
        for ipt,pt in enumerate(pt_vals):
            for iy,y in enumerate(y_vals):
                for iA4,A4 in enumerate(A4_vals):
                    bin = np.load(self.input_dir+grid_tag+'_A'+str(ic)+('{:03.2f}'.format(A4)) +'_pt{:02.1f}'.format(pt)+'_'+'y{:02.1f}'.format(y)+'_weighted{:05.3f}'.format(mass)+'.npy')
                    weight = proportions[ipt][iy][iA4]
                    if not mix.shape==bin[0].shape:
                        mix = bin[0]*weight
                        pt_edges = bin[1]
                        y_edges = bin[2]
                    else:
                        mix += bin[0]*weight 

        print(mix)

        xx, yy = np.meshgrid(pt_edges, y_edges)        
        plt.pcolormesh(xx, yy, mix.T)
        plt.show()
        plt.savefig('../test/mix.png')

##################################

distribution = np.empty((4,),dtype=object)
distribution[0] = np.zeros((4,9,2))
distribution[1] = np.array([0.0, 5.0, 10.0, 20.0])
distribution[2] = np.array([-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])
distribution[3] = np.array([-1.0, +1.0])

probs_vs_pt = [0.1, 0.4, 0.2, 0.1]
probs_vs_y = [0.1, 0.4, 0.3, 0.2, 0.2, 0.2, 0.3, 0.4, 0.1]

for ipt,pt in enumerate(distribution[1]): 
    for iy,y in enumerate(distribution[2]): 
        for iA4, A4 in enumerate(distribution[3]):
            prob_pol = 0.8 if (A4<0 and y>0) or (A4>0 and y<0) else 0.2
            prob = probs_vs_pt[ipt]*probs_vs_y[iy] * prob_pol
            distribution[0][ipt][iy][iA4] = prob

mixdata = MixData()
mixdata.mix(coeff=[0.,0.,0.,0.,-1.0], mass=80.000, distribution=distribution)
