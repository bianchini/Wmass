import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from template_parameters import params_test as params
from template_parameters import pdf_test

rebin=(1,1)
input_y_bins  = np.linspace(params['params_lep']['y_range'][0], 
                            params['params_lep']['y_range'][1],  
                            params['params_lep']['y_bins']/rebin[0]+1)
input_pt_bins = np.linspace(params['params_lep']['pt_range'][0], 
                            params['params_lep']['pt_range'][1],  
                            params['params_lep']['pt_bins']/rebin[1]+1)
output_shapes_pt = params['params_W']['pt']


dirname = 'TEST'

tot_weight = 0.0
pL = np.load('../data/'+dirname+'/mixed_dataset_pt0.0-2.0_y0.40-0.80_M80.000_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
weight = 0.0
for ipt,pt in enumerate(output_shapes_pt):
    if pt<2.0:
        weight += pdf_test(pt=pt,y=0.0)
pL *= weight
print pL.sum()
tot_weight += weight

pH = np.load('../data/'+dirname+'/mixed_dataset_pt2.0-4.0_y0.40-0.80_M80.000_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
weight = 0.0
for ipt,pt in enumerate(output_shapes_pt):
    if pt>=2.0 and pt<4.0:
        weight += pdf_test(pt=pt,y=0.0)
pH *= weight
print pH.sum()
tot_weight += weight

sumt = (pL+pH)/tot_weight
print sumt.sum()

p  = np.load('../data/'+dirname+'/mixed_dataset_pt0.0-4.0_y0.40-0.80_M80.000_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
print p.sum()

xx, yy = np.meshgrid(input_pt_bins, input_y_bins)        
plt.subplot(2, 2, 1)
plt.pcolormesh(yy, xx, pL)
plt.colorbar()

plt.subplot(2, 2, 2)
plt.pcolormesh(yy, xx, pH)
plt.colorbar()

plt.subplot(2, 2, 3)
plt.pcolormesh(yy, xx, p)
plt.colorbar()

plt.subplot(2, 2, 4)
plt.pcolormesh(yy, xx, p-sumt )
plt.colorbar()

plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
plt.savefig('overlap.png')
plt.close()

