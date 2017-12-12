import sys
sys.path.append('./')
sys.path.append('../python/')
 
from unfolder import Unfolder

from template_parameters import params_test

unfolder = Unfolder(input_dir='../data/', params=params_test, mass=80.000, num_events=100000,  fix=['A0', 'A4', 'pt_y'], interp_deg=1)

for itoy in range(1):
    #unfolder.toy_data(ntoy=itoy)
    unfolder.run()

