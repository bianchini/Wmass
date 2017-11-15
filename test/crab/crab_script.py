#!/usr/bin/env python
import os
import ROOT
#from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from crabhelper import inputFiles,runsAndLumis

# here I place the actual job
inFile = ROOT.TFile.Open(os.environ['CMSSW_BASE']+"/src/Wmass/data/"+"result.root") 
inTree = inFile.Get("tree") 
print inTree.GetEntries() 
outFile = ROOT.TFile.Open("tree.root", "RECREATE")
outTree = inTree.CloneTree(10,"")
print outTree.GetEntries() 
outTree.Write() 
outFile.Close()
inFile.Close()

#for fname in inputFiles():
#    inFile = ROOT.TFile.Open(fname)
#    inTree = inFile.Get("tree")
#    print inTree.GetEntries()

print "DONE"
os.system("ls -lR")

