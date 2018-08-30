#!/usr/bin/env python
import os

import sys
sys.path.append('./')
from sys import argv
 
from crabhelper import inputFiles
from tree_producer import TreeProducer

mass_mc = 80.419
masses = [mass_mc]

treeProducer = TreeProducer(DY='CC_FxFx', verbose=False, debug=False, filenames=inputFiles(), save_tree=False, save_histo1=False, save_histo2=False, save_histo3=False, save_histo4=True, masses=masses, postfix='test',  plot_vars_histo3=['eta','pt'])
treeProducer.run()

print "DONE"
os.system("ls -lR")
