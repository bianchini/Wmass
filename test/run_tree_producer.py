import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv
 
from tree_producer import TreeProducer

treeProducer = TreeProducer(DY='NC', verbose=False, debug=True, filenames=[])
treeProducer.run()
