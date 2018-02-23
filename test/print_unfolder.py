import os.path
from sys import argv

import pickle
from pprint import pprint

import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


save_dir = 'plots/'

jobs = [ 
    #['1e7_decorr_quad_rebin4-1', 10 ],
    #['1e7_decorr_cube_rebin4-1',  10 ],
    #['1e7_corr_quad_rebin4-1', 10 ],
    #['1e7_corr_cube_rebin4-1',  10 ],
    #['1e7_corr_cube_rebin4-1_morebins',   10],
    #['1e7_corr_tetra_rebin4-1', 10]
    #['1e7_decorr_quad_rebin4-1_fixmass',  10],
    #['1e7_decorr_cube_rebin4-1_fixmass',  10],    

    #['1e7_newtemp_pt0-2_corr_quad_rebin4-1',    10],
    #['1e7_newtemp_pt0-2_corr_cube_rebin4-1',    10],
    #['1e7_newtemp_pt0-2_corr_tetra_rebin4-1',   10],
    #['1e7_newtemp_pt0-2_corr_quad_rebin4-2',    10],
    #['1e7_newtemp_pt0-2_corr_cube_rebin4-2',    10],
    #['1e7_newtemp_pt0-2_corr_tetra_rebin4-2',   10],
    #['1e7_newtemp_pt0-2_corr_quad_rebin4-4',    10],
    #['1e7_newtemp_pt0-2_corr_cube_rebin4-4',    10],
    #['1e7_newtemp_pt0-2_corr_tetra_rebin4-4',   10],
    
    #['1e7_newtemp_pt0-2_corr_p2_rebin4-4',       100],
    #['1e7_newtemp_pt0-2_corr_p2p4_rebin4-4',     100],
    ['1e7_newtemp_pt0-2_corr_p2_rebin4-4_cov_fixPower',     100],
    ['1e7_newtemp_pt0-2_corr_p2p4_rebin4-4_cov_fixPower',   100],
    ['1e7_newtemp_pt0-2_corr_p2p4p6_rebin4-4_cov_fixPower', 100],
    #['1e7_newtemp_pt0-2_corr_p2p4p8_rebin4-4',   100],
    #['1e7_newtemp_pt0-2_corr_p2_rebin4-4_syst',       100],
    #['1e7_newtemp_pt0-2_corr_p2p4_rebin4-4_syst',     100],
    #['1e7_newtemp_pt0-2_corr_p2p4p8_rebin4-4_syst',   10],

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
    print "Fit: ", job[0]

    nrun = 0
    for j in range(job[1]):
        if not os.path.isfile('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl'):
            continue
        f = open('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl', 'r')
        res = pickle.load(f)
        if res['status'][0]==0:
            nrun += 1

    time = []
    ntime = np.zeros((nrun))

    edm = []
    nedm = np.zeros((nrun))

    amin = []
    namin = np.zeros((nrun))

    status = []
    nstatus = np.zeros((nrun))

    ndof = []
    nndof = np.zeros((nrun))

    pvalue = []
    npvalue = np.zeros((nrun))

    mass = []
    nmass = np.zeros((nrun))

    mass_err = []
    nmass_err = np.zeros((nrun))

    mass_pull = []
    nmass_pull = np.zeros((nrun))

    nrun = 0
    for j in range(job[1]):
        if not os.path.isfile('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl'):
            continue
        f = open('crab/crab_unfolder_'+job[0]+'/results/result_'+job[0]+'_'+str(j)+'.pkl', 'r')
        res = pickle.load(f)
        if res['status'][0]>0:
            continue
        for toy in range(res['ntoys']):        
            time.append( "{:02.1f}".format(res['time'][toy]/3600.) )
            ntime[nrun] = res['time'][toy]/3600.

            edm.append( "{:03.3f}".format(res['edm'][toy] / 1e-04) )
            nedm[nrun] = res['edm'][toy] / 1e-04

            amin.append( "{:03.0f}".format(res['amin'][toy]))
            namin[nrun] = res['amin'][toy]
            
            status.append( "{:01.0f}".format(res['status'][toy]) )
            nstatus[nrun] = res['status'][toy]

            ndof.append( "{:02.0f}".format(res['ndof'][toy]) )
            nndof[nrun] = res['ndof'][toy]

            mass.append( "{:05.3f}".format(res['mass']['fit'][toy]))
            nmass[nrun] = res['mass']['fit'][toy]

            mass_err.append( "{:05.3f}".format(res['mass']['err'][toy]))
            nmass_err[nrun] = res['mass']['err'][toy]

            mass_pull.append( "{:05.3f}".format( (res['mass']['fit'][toy]-80.000)/res['mass']['err'][toy] ))
            nmass_pull[nrun] =  (res['mass']['fit'][toy]-80.000)/res['mass']['err'][toy]

            pvalue.append("{:02.2f}".format(res['pvalue'][toy]))
            npvalue[nrun] = res['pvalue'][toy]
        nrun += 1
        f.close()

    print "\t - njobs    = ", nstatus.size
    print "\t - status   = ", nstatus.mean()#, "    ", status
    print "\t - edm*1e-4 = ", "{:03.3f}".format(nedm.mean())#,  "  ", edm
    print "\t - chi2     = ", "{:03.0f}".format(namin.mean())#,  "   ", amin
    print "\t - ndof     = ", "{:02.0f}".format(nndof.mean())#,  "   ", ndof
    print "\t - pvalue   = ", "{:02.2f}".format(npvalue.mean())#,  "   ", pvalue
    print "\t - time [h] = ", "{:02.1f}".format(ntime.mean())#,  "    ", time
    print "\t - mass     = ", "{:05.3f}".format(nmass.mean())#,  " ", mass
    print "\t - mass_err = ", "{:05.3f}".format(nmass_err.mean())#,  "  ", mass_err
    print "\t - mass_pull= ", "{:05.3f}".format(nmass_pull.mean()),  " +/- "+"{:05.3f}".format( np.std(nmass_pull)/math.sqrt(nmass_pull.size) ), " sigma="+"{:05.3f}".format( np.std(nmass_pull) )#, "  ", mass_pull

    n, bins, patches = plt.hist(nmass_pull, 20, normed=1, facecolor='g', alpha=0.75)    
    plt.xlabel('Pull')
    plt.ylabel('Toys')
    plt.title(job[0])
    plt.text(-3.0, 0.9, r'$\mu='+'{:05.3f}'.format(nmass_pull.mean())+' \pm  '+'{:05.3f}'.format( np.std(nmass_pull)/math.sqrt(nmass_pull.size) )+' ,\ \sigma='+'{:05.3f}'.format( np.std(nmass_pull) )+'$', fontsize=20)
    plt.axis([-4, 4, 0, 1])
    plt.grid(True)
    plt.show()
    plt.savefig(save_dir+'pull_mass_'+job[0]+'.png')
    plt.close()

    n, bins, patches = plt.hist(npvalue, 10, normed=1, facecolor='g', alpha=0.75)    
    plt.xlabel('p-value')
    plt.ylabel('Toys')
    plt.title(job[0])
    plt.text(-3.0, 0.9, r'p-value')
    plt.axis([0, 1, 0, 5])
    plt.grid(True)
    plt.show()
    plt.savefig(save_dir+'p-value_'+job[0]+'.png')
    plt.close()

    n, bins, patches = plt.hist(nmass, 20, normed=1, facecolor='g', alpha=0.75)    
    plt.xlabel('mass')
    plt.ylabel('Toys')
    plt.title(job[0])
    plt.text(79.950, 40, r'$\mu='+'{:05.3f}'.format( nmass.mean() )+'\;\;\; \mathrm{RMS}='+'{:05.3f}'.format( np.std(nmass) )+'$', fontsize=20)
    plt.axis([79.900, 80.100, 0, 50])
    plt.grid(True)
    plt.show()
    plt.savefig(save_dir+'mass_'+job[0]+'.png')
    plt.close()

