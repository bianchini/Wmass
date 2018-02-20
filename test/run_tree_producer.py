import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv

from tree_producer import TreeProducer

treeProducer = TreeProducer(DY='CC', verbose=False, debug=True, filenames=[], save_tree=False, save_histo=True)
treeProducer.run()

