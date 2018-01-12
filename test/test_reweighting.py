import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ul = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
p1 = np.array([])
p2 = np.array([])

A = 1.0
if False:
    p1 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A02.00_A10.00_A20.00_A30.00_A40.00.npy')
    p2 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A0'+'{:03.2f}'.format(A)+'_A10.00_A20.00_A30.00_A40.00.npy')
    A = 0.5
if False:
    p1 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A11.00_A20.00_A30.00_A40.00.npy')
    p2 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A1'+'{:03.2f}'.format(A)+'_A20.00_A30.00_A40.00.npy')
if False:
    p1 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A21.00_A30.00_A40.00.npy')
    p2 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A2'+'{:03.2f}'.format(A)+'_A30.00_A40.00.npy')
if False:
    p1 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A31.00_A40.00.npy')
    p2 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A3'+'{:03.2f}'.format(A)+'_A40.00.npy')
if True:
    p1 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A30.00_A42.00.npy')
    p2 = np.load('../root//grid_lab_pt6.0_y0.00_M80.000_A00.00_A10.00_A20.00_A30.00_A4'+'{:03.2f}'.format(A)+'.npy')
    A = 0.5

p2a = np.zeros(ul[0].shape) 
p2b = np.zeros(ul[0].shape) 
p2a += (ul[0]*(1-A) + p1[0]*A)
p2b += (-ul[0] + p1[0])*2.0

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

plt.subplot(2, 2, 4)
plt.pcolormesh(yy, xx, p2b)
plt.colorbar()

plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.show()
plt.savefig('diff.png')
 
