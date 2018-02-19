import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv

import os.path
from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

from tree_utils import draw_slice 

for coeff in ['A0','A1','A2','A3','A4','A5', 'A6', 'A7']:
#for coeff in ['A0']:
    draw_slice(fname='./tree.root', var='Wdress', coeff=coeff, weight_name=0)
