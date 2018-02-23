#!/usr/bin/env python
import os

import sys
sys.path.append('./')
from sys import argv
 
from crabhelper import inputFiles
from tree_producer import TreeProducer

treeProducer = TreeProducer(DY='CC', verbose=False, debug=False, filenames=inputFiles(), save_tree=False, save_histo=True, postfix='test')
treeProducer.run()

print "DONE"
os.system("ls -lR")
