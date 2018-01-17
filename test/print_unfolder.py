import os.path
from sys import argv

import pickle
from pprint import pprint

import numpy as np

jobs = [ 
        #['1e7_decorr_taylor2_rebin10-1', 10 ],
        #['1e7_decorr_taylor2_rebin5-1',  10 ],
        #['1e7_decorr_taylor2_rebin4-1',  10 ],
        #['1e7_decorr_taylor2_rebin2-1',  10 ],
        #['1e7_decorr_taylor2_rebin1-1',  10 ],
        #['1e7_decorr_taylor3_rebin10-1', 10 ],
        #['1e7_decorr_taylor3_rebin5-1',  10 ],
        #['1e7_decorr_taylor3_rebin4-1',  10 ],
        #['1e7_decorr_taylor3_rebin2-1',  10 ],
        #['1e7_decorr_taylor3_rebin1-1',  10 ],
        #['1e7_decorr_taylor2_rebin10-2', 10 ],
        #['1e7_decorr_taylor2_rebin5-2',  10 ],
        #['1e7_decorr_taylor2_rebin4-2',  10 ], 
        #['1e7_decorr_taylor2_rebin2-2',  10 ], 
        #['1e7_decorr_taylor2_rebin1-2',  10 ], 
        #['1e7_decorr_taylor3_rebin10-2', 10 ], 
        #['1e7_decorr_taylor3_rebin5-2',  10 ], 
        #['1e7_decorr_taylor3_rebin4-2',  10 ], 
        #['1e7_decorr_taylor3_rebin2-2',  10 ], 
        #['1e7_decorr_taylor3_rebin1-2',  10 ], 
        ]

for ijob,job in enumerate(jobs):
    print "............................................"
    print "Fit: ", job

    nrun = 0
    for j in range(job[1]):
        if not os.path.isfile('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl'):
            continue
        nrun += 1

    time = []
    ntime = np.zeros((nrun))

    edm = []
    nedm = np.zeros((nrun))

    amin = []
    namin = np.zeros((nrun))

    status = []
    nstatus = np.zeros((nrun))

    pvalue = []
    npvalue = np.zeros((nrun))

    mass = []
    nmass = np.zeros((nrun))

    mass_err = []
    nmass_err = np.zeros((nrun))

    nrun = 0
    for j in range(job[1]):
        if not os.path.isfile('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl'):
            continue
        f = open('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl', 'r')
        res = pickle.load(f)
        for toy in range(res['ntoys']):        
            time.append( "{:02.1f}".format(res['time'][toy]/3600.) )
            ntime[nrun] = res['time'][toy]/3600.

            edm.append( "{:03.3f}".format(res['edm'][toy] / 1e-04) )
            nedm[nrun] = res['edm'][toy] / 1e-04

            amin.append( "{:03.0f}".format(res['amin'][toy]))
            namin[nrun] = res['amin'][toy]
            
            status.append( "{:01.0f}".format(res['status'][toy]) )
            nstatus[nrun] = res['status'][toy]

            mass.append( "{:05.3f}".format(res['mass']['fit'][toy]))
            nmass[nrun] = res['mass']['fit'][toy]

            mass_err.append( "{:05.3f}".format(res['mass']['err'][toy]))
            nmass_err[nrun] = res['mass']['err'][toy]

            pvalue.append("{:02.2f}".format(res['pvalue'][toy]))
            npvalue[nrun] = res['pvalue'][toy]
        nrun += 1
        f.close()

    print "\t - status   = ", nstatus.mean(), "    ", status
    print "\t - edm*1e-4 = ", "{:03.3f}".format(nedm.mean()),  "  ", edm
    print "\t - amin     = ", "{:03.0f}".format(namin.mean()),  "   ", amin
    print "\t - pvalue   = ", "{:02.2f}".format(npvalue.mean()),  "   ", pvalue
    print "\t - time [h] = ", "{:02.1f}".format(ntime.mean()),  "    ", time
    print "\t - mass     = ", "{:05.3f}".format(nmass.mean()),  " ", mass
    print "\t - mass_err = ", "{:05.3f}".format(nmass_err.mean()),  "  ", mass_err

