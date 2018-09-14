import sys
sys.path.append('./')
sys.path.append('../python/')
from sys import argv

from tree_producer import TreeProducer

mass_mc = 80.419
masses = [mass_mc]

treeProducer = TreeProducer(DY='CC_FxFx', verbose=False, debug=True, filenames=[], save_tree=False, save_histo1=False, save_histo2=False, save_histo3=False, save_histo4=False, save_histo5=True, masses=masses, plot_vars_histo3=['eta','pt'])
treeProducer.run()

