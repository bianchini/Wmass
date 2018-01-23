from sys import argv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pylab import rcParams
rcParams['figure.figsize'] = 8,7

dirname = 'TEST'
pt = argv[1]
do_deriv = int(argv[2])

y = '0.00'
mass = '80.000'

files = [ 
    ['mass', 80.500],
    ['A0', 2.0], 
    ['A1', 1.0], 
    ['A2', 1.0], 
    ['A3', 1.0], 
    ['A4', 2.0], 
    ]

ul = np.load('../root/'+dirname+'/grid_lab_pt'+pt+'_y'+y+'_M'+mass+'_A00.00_A10.00_A20.00_A30.00_A40.00.npy')
xx, yy = np.meshgrid(ul[2], ul[1])        

for iif,f in enumerate(files):
    A = f[1]
    test_mass = mass
    if f[0]=='mass' and do_deriv:
        test_mass = '{:05.3f}'.format(A)
    name = '../root/'+dirname+'/grid_lab_pt'+pt+'_y'+y+'_M'+test_mass
    for c in ['A0', 'A1', 'A2', 'A3', 'A4']:
        if c==f[0]:
            name += '_'+c+'{:03.2f}'.format(A)
        else:
            name += '_'+c+'{:03.2f}'.format(0.00)
    name += '.npy'
    p = np.load(name)
    deriv = np.zeros(p[0].shape)
    if f[0]=='mass':
        A = (f[1]-float(mass))
    if do_deriv:
        deriv += (-ul[0] + p[0])*(1./A)
    else:
        deriv += p[0]

    plt.subplot(3, 2, iif+1)
    plt.pcolormesh(yy, xx, deriv)    
    plt.xlabel('eta')
    plt.ylabel('pt')
    if do_deriv:
        plt.title('df/d('+f[0]+'), pt='+pt)
    else:
        plt.title(f[0]+', pt='+pt)
    plt.colorbar(format='%.0e')
    print f[0], ' ==> ', np.absolute(deriv).sum()

plt.subplots_adjust(wspace=0.55, hspace=0.6)
plt.show()
if do_deriv:
    plt.savefig('derivatives_pt'+pt+'_y'+y+'_M'+mass+'.png')
else:
    plt.savefig('templates_pt'+pt+'_y'+y+'_M'+mass+'.png')
 
