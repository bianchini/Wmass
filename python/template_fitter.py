import ROOT
import math
import sys
import os
import pickle
import time
from pprint import pprint
import numpy as np
import numpy.random as ran

from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from array import array;
import copy


class TemplateFitter:

    def __init__(self, DY='CC', charge='Wplus', var='WpreFSR', job_name='TEST', mc_mass=80.419):
        self.in_dir  = os.environ['CMSSW_BASE']+'/src/Wmass/data/'
        self.out_dir = os.environ['CMSSW_BASE']+'/src/Wmass/test/'
        self.job_name = job_name

        self.res_coeff = np.load(open(self.in_dir+'fit_results_'+DY+'_'+charge+'_all_A0-7.pkl', 'r'))        

        templates = np.load(self.in_dir+'template_'+charge+'_'+var+'_val.npz')
        templates_files = {'template': 0, 'masses': 1, 'bins_qt' : 2, 
                           'bins_y'  : 3, 'coeff' : 4, 'bins_eta': 5, 
                           'bins_pt' : 6, 'mc_acceptances' : 7 }
        for key,p in templates_files.items():             
            template_file = templates['arr_'+str(p)]
            print ('Adding file arr_'+str(p)+' with shape...'), template_file.shape, (' to self with key... '+key)
            setattr(self, key,template_file)
            size = -1
            # for template, save the size of (pt,eta) plane 
            if key=='template':
                size = template_file[0][0][0][0].size
            # for bins, save len(bins)-1, since bins are edges
            elif 'bins' in key:
                size = template_file.size - 1
            # for coeff, save len(coeff)-2, since coeff[-1] = MC and coeff[-2] = UL
            elif key=='coeff':
                size = template_file.size - 2
            else:
                size = template_file.size
            setattr(self, key+'_size', size)

        mc_mass_index = np.where(self.masses==mc_mass)[0][0]
        self.mc = self.template.sum(axis=(1,2))[mc_mass_index,-1]

        self.save_template_snapshot(data=self.mc, title='MC', tag='mc')
        self.fit_results = {}
        self.out_file_res = open(self.out_dir+'result_'+DY+'_'+charge+'_'+job_name+'.pkl','wb')
        pickle.dump(self.fit_results, self.out_file_res)
        self.out_file_res.close()

    def save_template_snapshot(self, data=np.array([]), title='', tag=''):
        xx,yy = np.meshgrid(self.bins_pt, self.bins_eta)        
        plt.pcolormesh(yy, xx, data)
        plt.colorbar(format='%.0e')
        plt.axis([self.bins_eta.min(), self.bins_eta.max(), self.bins_pt.min(), self.bins_pt.max()])
        plt.title(title)
        plt.show()
        plt.savefig(self.out_dir+'snapshot_'+tag+'_'+self.job_name+'.png')
        plt.close('all')
