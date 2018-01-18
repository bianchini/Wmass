import numpy as np
np.set_printoptions(threshold=np.inf)
import copy

mass = 80.250

p0 = np.load('../data/mixed_dataset_pt0.0-20.0_y0.00-y0.40_M79.500_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
p1 = np.load('../data/mixed_dataset_pt0.0-20.0_y0.00-y0.40_M80.000_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
p2 = np.load('../data/mixed_dataset_pt0.0-20.0_y0.00-y0.40_M79.500_A00.00_A10.00_A20.00_A30.00_A40.00.npy')

x = np.array([79.500, 80.000, 80.500])

y = np.zeros( (3, p0.size) )
for i in range(p0.size):
    y[0][i] = p0.flatten()[i]
    y[1][i] = p1.flatten()[i]
    y[2][i] = p2.flatten()[i]
print p0.shape

fit = np.polyfit(x, y, deg=2)

p = np.zeros(p0.size)
for i in range(p.size):
    p[i] = fit[0][i]*mass*mass + fit[1][i]*mass + fit[2][i]
p_new = p.reshape(p0.shape)
print p_new.shape
print (p_new-p1)


