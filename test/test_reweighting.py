import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

A4 = 1.0
ul = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
p1 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A30.00_A42.00.npy')
p2 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A30.00_A4'+'{:03.2f}'.format(A4)+'.npy')


p2a = np.zeros(ul[0].shape) 
p2a += (ul[0]*(1-A4/2) + p1[0]*A4/2)

print p2a.sum(), '<=>' , p2[0].sum()

xx, yy = np.meshgrid(ul[2], ul[1])        
plt.subplot(2, 2, 1)
plt.pcolormesh(yy, xx, p2a)
plt.colorbar()
#plt.show()
#plt.savefig('weighted.png')

plt.subplot(2, 2, 2)
plt.pcolormesh(yy, xx, p2[0])
plt.colorbar()
#plt.show()
#plt.savefig('plain.png')

diff = p2a-p2[0]
plt.subplot(2, 2, 3)
plt.pcolormesh(yy, xx, diff)
plt.colorbar()
#plt.show()
#plt.savefig('diff.png')

plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
plt.savefig('diff.png')
 
