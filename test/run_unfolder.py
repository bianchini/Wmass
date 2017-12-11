import sys
sys.path.append('./')
sys.path.append('../python/')
 
from unfolder import Unfolder

from template_parameters import params_test

unfolder = Unfolder(input_dir='../data/', params=params_test, num_events=1000000)
unfolder.run()

