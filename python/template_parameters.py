import copy
import math
import numpy as np

Parameters = {

    'params_CS' : {
        'mass' : 80.000,
        'coeff' : [0.,0.,0.,0.,0.]
        },
    
    'params_lep' : {
        'pt_range' : [25.0, 50.0],
        'pt_bins' : 100,
        'y_range' : [-2.5, 2.5],
        'y_bins'  : 100,
        },
    
    'params_W' : {
        'pt' : np.linspace(0.0, 20., 2.0),
        'y'  : np.linspace(-3.0, 3.0, 0.5),
        'mass' : np.linspace(78.000, 82.000, 2),
        'A0' : np.array([0.0]),
        'A1' : np.array([0.0]),
        'A2' : np.array([0.0]),
        'A3' : np.array([0.0]),
        'A4' : np.array([0.0]),
        },    

    'params_template' : {
        'pt' : np.array([0.0, 2.0, 5.0, 10., 20.,]),
        'y'  : np.linspace(-3.0, 3.0, 31),
        'mass' : np.linspace(78.000, 82.000, 0.500),
        'A0' : np.array([0.0]),
        'A1' : np.array([0.0]),
        'A2' : np.array([0.0]),
        'A3' : np.array([0.0]),
        'A4' : np.array([0.0]),
        },    
}

def pdf_test(pt,y):

    val = 1.0

    # pt
    pt_max = 5.0
    lambdaQCD = 0.200
    val *= math.exp(-(pt+lambdaQCD)/pt_max)*(pt+lambdaQCD)/pt_max/pt_max

    # y        
    y_width = 2.5
    val *= 1./math.sqrt(2*math.pi)/(y_width)*math.exp(-0.5*y*y/y_width/y_width)
    
    return val


def coefficients_test(pt=0.0,y=0.0):
    coeff=[0.,0.,0.,0.,0.]
    return coeff


def accept_point(coeff=[], verbose=False):
    #return True
    if coeff in  [ [0.0, 0.0, 0.0, 0.0, 0.0],
                   [2.0, 0.0, 0.0, 0.0, 0.0],
                   [0.0, 1.0, 0.0, 0.0, 0.0],
                   [0.0, 0.0, 1.0, 0.0, 0.0],
                   [0.0, 0.0, 0.0, 1.0, 0.0],
                   [0.0, 0.0, 0.0, 0.0, 2.0]
                   ]:                  
        return True
    else:
        if verbose:
            print "Point", coeff, "is excluded: continue" 
        return False

def bin_ndarray(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.

    Number of output dimensions must match number of input dimensions and 
        new axes must divide old ones.

    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)

    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]

    """
    operation = operation.lower()
    if not operation in ['sum', 'mean']:
        raise ValueError("Operation not supported.")
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d,c in zip(new_shape,
                                                  ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        op = getattr(ndarray, operation)
        ndarray = op(-1*(i+1))
    return ndarray



params_test = copy.deepcopy(Parameters)
params_test['params_W']['pt'] = np.linspace(0.0, 32.0, 33)
#params_test['params_W']['pt'] = np.array([6.0])
params_test['params_W']['y'] = np.array([0.0])
#params_test['params_W']['mass'] = np.array([79.500,80.000,80.500])
params_test['params_W']['mass'] = np.array([79.500, 80.500])
params_test['params_W']['A0'] = np.array([ 0.0, 2.0 ]) 
params_test['params_W']['A1'] = np.array([ 0.0, 1.0 ]) 
params_test['params_W']['A2'] = np.array([ 0.0, 1.0 ])
params_test['params_W']['A3'] = np.array([ 0.0, 1.0 ]) 
params_test['params_W']['A4'] = np.array([ 0.0, 2.0 ])

params_test['params_template']['pt'] = np.append(np.linspace(0.0, 20.0, 6), np.array([26.0, 32.0]))
#params_test['params_template']['pt'] = np.linspace(0.0, 20.0, 2)
params_test['params_template']['y']  = np.linspace(0.0, 3.6, 10)
#params_test['params_template']['y']  = np.linspace(0.0, 3.6, 2)
params_test['params_template']['mass'] = params_test['params_W']['mass']
params_test['params_template']['A0'] = params_test['params_W']['A0']
params_test['params_template']['A1'] = params_test['params_W']['A1']
params_test['params_template']['A2'] = params_test['params_W']['A2']
params_test['params_template']['A3'] = params_test['params_W']['A3']
params_test['params_template']['A4'] = params_test['params_W']['A4']
