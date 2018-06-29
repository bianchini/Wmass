import os

path = '/gpfs/ddn/srm/cms/store/user/bianchi/CRAB_PrivateMC/TEST/'
dirs = []

# Asimov
#dirs.extend(['180626_153513', '180626_153606', '180626_153634', '180626_153729', '180626_153757', '180626_153934', '180626_154003', '180626_154035', '180626_154104', '180626_154134', '180626_154203', '180626_154237', '180626_154308', '180627_084120', '180627_095334', '180627_095521', '180627_095627', '180627_095728', '180627_095834'])
#dirs.extend(['180626_154035', '180627_095521'])

# Random
#dirs.extend(['180627_140453'],
#'180627_140519', '180627_140552']
#)

#dirs.extend(['180628_082425', '180628_082452'])
#dirs.extend(['180628_134051', '180628_134122'])

#dirs.extend([#'180628_152329',
             #'180628_152403',
             #'180628_152440',
             #'180628_152514',
             #'180628_152550',
             #'180628_152625',
             #'180628_152703',
             #'180628_152746',
             #'180628_154106',
             #'180628_154140'
            #'180629_102823'
             #])

dirs.extend(['180629_133044'])

for d in dirs:
    print 'Directory '+d
    root_files = [f for f in os.listdir(path+'/'+d+'/0000/') if os.path.isfile(os.path.join(path+'/'+d+'/0000/', f)) and '.root' in f]
    if len(root_files)>0:
        root_name_split = root_files[0].split('_')
        idx = root_files[0].find( root_name_split[-1] )
        name = root_files[0][:(idx-1)]
        print '\tHadding ROOT file(s) '+name
        os.system('hadd -f '+name+'.root '+path+'/'+d+'/0000/'+name+'_*.root')
    print '\tCopy plots to local'
    for ext in ['png', 'C']:
        os.system('cp '+path+'/'+d+'/0000/*.'+ext+' .')
