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
    return val

    # pt
    if pt >= 0. and pt < 2.0:
        val *= 0.025
    elif pt >= 2.0 and pt < 4.0:
        val *= 0.050
    elif pt >= 4.0 and pt < 6.0:
        val *= 0.055
    elif pt >= 6.0 and pt < 8.0:
        val *= 0.05
    elif pt >= 8.0 and pt < 10.0:
        val *= 0.04
    elif pt >= 10.0 and pt < 12.0:
        val *= 0.035
    elif pt >= 12.0 and pt < 14.0:
        val *= 0.03
    elif pt >= 14.0 and pt < 16.0:
        val *= 0.025
    elif pt >= 14.0 and pt < 18.0:
        val *= 0.02
    elif pt >= 18.0 and pt < 20.0:
        val *= 0.015
    else:
        val *= 0.01

    # y        
    y_width = 1.5
    val *= 1./math.sqrt(2*math.pi)/(y_width)*math.exp(-0.5*y*y/y_width/y_width)
    
    # coeff
    #val *= 1.0 if c4==-y/2.5 else 0.0
    val *= (1.0 if c4*y <= 0. else 0.0)
    return val


params_test = copy.deepcopy(Parameters)
params_test['params_W']['pt'] = np.array([0.0]) 
params_test['params_W']['y'] = np.linspace(-2.5,2.5,201) 
params_test['params_W']['mass'] = np.array([80.000]) 
#params_test['params_W']['A0'] = np.array([0.0]) 
#params_test['params_W']['A1'] = np.array([0.0]) 
#params_test['params_W']['A2'] = np.array([0.0])
#params_test['params_W']['A3'] = np.array([0.0]) 
params_test['params_W']['A0'] = np.array([0.0]) 
