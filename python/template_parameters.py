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
        'pt_bins' : 50,
        'y_range' : [-2.5, 2.5],
        'y_bins' : 201,
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

def pdf_test(pt,y, c0, c1, c2, c3, c4):

    val = 1.0

    # pt
    pt_max = 5.0
    lambdaQCD = 0.200
    val *= math.exp(-(pt+lambdaQCD)/pt_max)*(pt+lambdaQCD)/pt_max/pt_max

    # y        
    y_width = 1.5
    val *= 1./math.sqrt(2*math.pi)/(y_width)*math.exp(-0.5*y*y/y_width/y_width)
    
    # coeff
    #val *= 1.0 if c4==-y/2.5 else 0.0
    #val *= (1.0 if c4*y <= 0. else 0.0)
    return val


params_test = copy.deepcopy(Parameters)
params_test['params_W']['pt'] = np.linspace(0.0, 10., 11) 
params_test['params_W']['y'] = np.array([0.0])
params_test['params_W']['mass'] = np.array([80.000])
params_test['params_W']['A0'] = np.array([0.0]) 
params_test['params_W']['A1'] = np.array([0.0]) 
params_test['params_W']['A2'] = np.array([0.0])
params_test['params_W']['A3'] = np.array([0.0]) 
params_test['params_W']['A4'] = np.array([0.0])

params_test['params_template']['pt'] = np.array([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
params_test['params_template']['y']  = np.array([-0.2, 0.0,  0.2])
params_test['params_template']['mass'] = np.array([80.000])#params_test['params_W']['mass']
params_test['params_template']['A0'] = params_test['params_W']['A0']
params_test['params_template']['A1'] = params_test['params_W']['A1']
params_test['params_template']['A2'] = params_test['params_W']['A2']
params_test['params_template']['A3'] = params_test['params_W']['A3']
params_test['params_template']['A4'] = params_test['params_W']['A4']
