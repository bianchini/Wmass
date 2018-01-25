import os
from sys import argv

debug = False

if argv[1]=='syst':
    if not os.path.isdir(os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST_syst_p1'):
        print 'Replace data/'
        if os.path.isdir(os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST'):
            print 'mv '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST ~/data/mixed_TEST_backup'
        if not debug:
            os.system('mv '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST ~/data/mixed_TEST_backup')
        print 'mv ~/data/mixed_TEST_nominal '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST'
        if not debug:
            os.system('mv ~/data/mixed_TEST_nominal '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST')
        print 'mv ~/data/mixed_TEST_syst_p1 '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST_syst_p1'
        if not debug:
            os.system('mv ~/data/mixed_TEST_syst_p1 '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST_syst_p1')

if argv[1]=='restore':
    if os.path.isdir(os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST_syst_p1'):
        print 'Restore data/'
        print 'mv '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST ~/data/mixed_TEST_nominal'
        if not debug:
            os.system('mv '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST ~/data/mixed_TEST_nominal' )
        print 'mv '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST_syst_p1 ~/data/mixed_TEST_syst_p1' 
        if not debug:
            os.system('mv '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST_syst_p1 ~/data/mixed_TEST_syst_p1')
        print 'mv ~/data/mixed_TEST_backup '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST'
        if not debug:
            os.system('mv ~/data/mixed_TEST_backup '+os.environ['CMSSW_BASE']+'/src/Wmass/data/TEST' )
