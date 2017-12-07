import ROOT
import math
import sys
import numpy as np
import numpy.random as ran

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from template_parameters import accept_point

# Class that computes grid (pt,y) in the lab
class TemplateMaker:

    # initialise
    def __init__(self, out_dir='../test/'):
        print "Initialize TemplateMaker"
        self.out_dir = out_dir
        self.save_plots = True
        self.width = 2.0
        return

    # specify the mesh size
    def set_CS_grid_mesh(self, grid_px_cos=500, grid_px_phi=10):
        self.grid_px_cos = grid_px_cos
        self.grid_px_phi = grid_px_phi        

        # the uniform grid in cos and phi (CS)
        self.grid_cos = np.linspace(-1., 1.,      grid_px_cos, dtype=np.float)
        self.grid_phi = np.linspace( 0., 2*np.pi, grid_px_phi, dtype=np.float)
        self.grid_indexes = np.arange(0, len(self.grid_cos)*len(self.grid_phi), 1)
        
    # compute dsigma/dphidcos in the CS frame
    def angular_pdf(self, x, y, coeff=[]):
        UL = (1.0 + x*x)
        L = 0.5*(1-3*x*x)
        T = 2*x*np.sqrt(1-x*x)*np.cos(y) 
        I = 0.5*(1-x*x)*np.cos(2*y)
        A = np.sqrt(1-x*x)*np.cos(y)
        P = x
        return 3./16./math.pi * ( UL + coeff[0]*L + coeff[1]*T + coeff[2]*I + coeff[3]*A + coeff[4]*P)

    # NR Breit-Wigner
    def BW(self, x, mass):
        bw = 1.0/math.pi/self.width/(1+math.pow((x-mass)/self.width,2.0))
        return bw

    # sample mass from BW
    def sample_BW(self, mass, truncate=True):
        x = np.random.standard_cauchy()
        check = lambda x : (x<3.0 and x>-3.0) 
        while truncate and not check(x): 
            x = np.random.standard_cauchy()
        return (x*self.width + mass)

    # (cos,phi,E) --> (E,px,py,pz)
    def boost_to_lab(self, x=np.array([0.,0.,0.]), y=0.0, pt=0.0):       
        cos = x[0]
        phi = x[1]
        E = x[2]

        px = E*math.sqrt(1-cos*cos)*math.cos(phi)  
        py = E*math.sqrt(1-cos*cos)*math.sin(phi)  
        pz = E*cos

        beta_z = math.tanh(y)
        M = 2.0*E
        beta_t = math.sqrt( pt*pt/(M*M + pt*pt) * (1-beta_z*beta_z) )
        gamma = 1./math.sqrt(1 - beta_z*beta_z - beta_t*beta_t)
        
        boost = np.array([[gamma, -beta_t*gamma , 0.0, -beta_z*gamma],
                          [-beta_t*gamma*gamma/math.sqrt(1+gamma*gamma*beta_t*beta_t), math.sqrt(1+gamma*gamma*beta_t*beta_t), 0.0, beta_t*beta_z*gamma*gamma/math.sqrt(1+gamma*gamma*beta_t*beta_t)],
                          [0., 0., 1., 0.],
                          [-beta_z*gamma/math.sqrt(1+gamma*gamma*beta_t*beta_t), 0, 0, gamma/math.sqrt(1+gamma*gamma*beta_t*beta_t)]
                          ] )
        x_CS = np.array([ E, px, py, pz ])
        x_lab = np.linalg.tensorsolve(boost, x_CS)
        return x_lab

    # fill (cos,phi) grid in CS frame
    def make_grid_CS(self, coeff=[0., 0., 0., 0., 0.], mass=80., ntoys=1000):
        # number of toys used to fill the template
        self.ntoys=ntoys

        # array that stores the rnd values
        rnd_grid = np.zeros((ntoys,3))
        xx, yy = np.meshgrid( self.grid_cos, self.grid_phi )        

        print "Evaluating the pdf over the grid..."
        pdf_evals = self.angular_pdf(xx,yy, coeff=coeff).flatten()
        pdf_evals_norm = np.sum(pdf_evals)    
        pdf_evals *= (1./pdf_evals_norm if pdf_evals_norm>0. else 1.0)        

        # these are the values chosen by the pdf
        rnd = np.random.choice( self.grid_indexes, size=ntoys, p=pdf_evals )

        print "Make the grid in the CS frame"
        # fill the rnd_grid
        for i,idx in enumerate(rnd):
            idx_cos =  idx%self.grid_px_cos
            idx_phi =  idx/self.grid_px_cos
            # cos
            rnd_grid[i][0] = self.grid_cos[idx_cos]
            # phi
            rnd_grid[i][1] = self.grid_phi[idx_phi]
            # E
            rnd_grid[i][2] = self.sample_BW(mass)*0.5
        return rnd_grid

    # fill (pt,y) template and save output as .npy
    def make_grid_lab(self, 
                      params = {},
                      boost_vecs=[], weights_mass=[], 
                      coeff=[0.,0.,0.,0.,0.], mass=80., ntoys=100000):

        # read from params how many templates we need
        boost_vecs = []
        for ipt,pt in enumerate(params['params_W']['pt']):
            for iy,y in enumerate(params['params_W']['y']):
                boost_vecs.append([pt,y])
        n_weights_mass = len(params['params_W']['mass'])
        n_weights_A0 = len(params['params_W']['A0'])
        n_weights_A1 = len(params['params_W']['A1'])
        n_weights_A2 = len(params['params_W']['A2'])
        n_weights_A3 = len(params['params_W']['A3'])
        n_weights_A4 = len(params['params_W']['A4'])
        n_weights = n_weights_mass*n_weights_A0*n_weights_A1*n_weights_A2*n_weights_A3*n_weights_A4

        grid_tag = 'M{:05.3f}'.format(mass)
        for ic,c in enumerate(coeff):
            grid_tag += '_A'+str(ic)+('{:03.2f}'.format(c))

        rnd_grid_CS = self.make_grid_CS( coeff=coeff, mass=mass, ntoys=ntoys )
        np.save(self.out_dir+'/grid_CS_'+grid_tag, rnd_grid_CS)

        if self.save_plots:
            self.plot(data_x=rnd_grid_CS[:,0], label_x='cos(theta)',              
                      data_y=rnd_grid_CS[:,1], label_y='phi', nbins=30, 
                      weights=np.ones(ntoys),
                      title=( 'CS frame, M={:05.3f}, A=[{:03.2f},{:03.2f},{:03.2f},{:03.2f},{:03.2f}]'.format(mass, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4])),     
                      name='grid_CS_'+grid_tag)

        # loop over pt,y values
        for boost_vec in boost_vecs:
            pt, y = boost_vec

            print ("Make the grid in the lab for (pt,y)=(%s,%s)" % (pt,y))            
            boost_tag = 'pt{:02.1f}'.format(pt)+'_'+'y{:03.2f}'.format(y)

            # fill the template
            setattr(self, 'rnd_grid_lab_'+boost_tag,  np.zeros((self.ntoys, 2+n_weights)) )

            # needed to reweight sample to various Ai,M
            weight_tensor = np.zeros((n_weights_mass,n_weights_A0,n_weights_A1,n_weights_A2,n_weights_A3,n_weights_A4), dtype=int)

            # fill the templates toy-by-toy
            for itoy in range(self.ntoys):
                toy_CS = rnd_grid_CS[itoy]
                angle_pdf = self.angular_pdf(toy_CS[0],toy_CS[1], coeff=coeff )

                x_lab = self.boost_to_lab( toy_CS, y=y, pt=pt)                        
                this_toy = getattr(self, 'rnd_grid_lab_'+boost_tag)[itoy] 
                this_toy[0] = math.sqrt(x_lab[0]*x_lab[0] - x_lab[3]*x_lab[3]) if x_lab[0]>abs( x_lab[3] ) else 0. # pt_lab
                this_toy[1] = 0.5*math.log((1+x_lab[3]/x_lab[0])/(1-x_lab[3]/x_lab[0])) if  x_lab[0]>abs( x_lab[3] ) else 10. # eta lab

                index = 0
                for im,m in enumerate(params['params_W']['mass']):
                    # mass weight is ratio of two BW
                    weight_mass = self.BW(2.0*toy_CS[2], m)/self.BW(2.0*toy_CS[2], mass)
                    for iA0,A0 in enumerate(params['params_W']['A0']):
                        for iA1,A1 in enumerate(params['params_W']['A1']):
                            for iA2,A2 in enumerate(params['params_W']['A2']):
                                for iA3,A3 in enumerate(params['params_W']['A3']):
                                    for iA4,A4 in enumerate(params['params_W']['A4']):

                                        if not accept_point(coeff=[A0,A1,A2,A3,A4]):
                                            continue
                                        # weight angle is the ratio between angular pdfs
                                        weight_angle = self.angular_pdf(toy_CS[0], toy_CS[1], coeff=[A0, A1, A2, A3, A4] )/angle_pdf
                                        this_toy[2+index] = weight_mass*weight_angle
                                        #print im, iA0, iA1, iA2, iA3, iA4, index, weight_angle
                                        weight_tensor[im][iA0][iA1][iA2][iA3][iA4] = index
                                        index += 1


            # save the histogram
            for im,m in enumerate(params['params_W']['mass']):
                for iA0,A0 in enumerate(params['params_W']['A0']):
                    for iA1,A1 in enumerate(params['params_W']['A1']):
                        for iA2,A2 in enumerate(params['params_W']['A2']):
                            for iA3,A3 in enumerate(params['params_W']['A3']):
                                for iA4,A4 in enumerate(params['params_W']['A4']):

                                    if not accept_point(coeff=[A0,A1,A2,A3,A4]):
                                        continue
                                    template = np.histogram2d( getattr(self, 'rnd_grid_lab_'+boost_tag)[:,1], getattr(self, 'rnd_grid_lab_'+boost_tag)[:,0],
                                                               weights=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,2+weight_tensor[im][iA0][iA1][iA2][iA3][iA4]],
                                                               normed=True,
                                                               bins=[params['params_lep']['y_bins'], params['params_lep']['pt_bins']],
                                                               range=[ params['params_lep']['y_range'], params['params_lep']['pt_range']]
                                                               )

                                    # save to disk
                                    out_name = boost_tag+'_M{:05.3f}'.format(m)
                                    out_name += ('_A0{:03.2f}'.format(A0)+'_A1{:03.2f}'.format(A1)+'_A2{:03.2f}'.format(A2)+'_A3{:03.2f}'.format(A3)+'_A4{:03.2f}'.format(A4))
                                    print "Saving histogram for "+out_name
                                    np.save(self.out_dir+'/grid_lab_'+out_name, template)

                                    if self.save_plots:            
                                        self.plot(data_x=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,1], label_x='eta',              
                                                  data_y=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,0], label_y='pt',
                                                  weights=getattr(self, 'rnd_grid_lab_'+boost_tag)[:, 2+weight_tensor[im][iA0][iA1][iA2][iA3][iA4] ],
                                                  title=( 'Lab frame, M={:05.3f}, A=[{:03.2f},{:03.2f},{:03.2f},{:03.2f},{:03.2f}] (pt,y)=({:02.1f},{:02.1f})'.format(m, A0, A1, A2, A3, A4, pt, y)), 
                                                  xmin=params['params_lep']['y_range'][0], xmax=params['params_lep']['y_range'][1], ymin=params['params_lep']['pt_range'][0], ymax=params['params_lep']['pt_range'][1], nbins=[params['params_lep']['y_bins'], params['params_lep']['pt_bins']],
                                                  name='/grid_lab_'+out_name)

    # plot (pt,y) for each .npy
    def plot(self, data_x=np.array([]), label_x='eta', data_y=np.array([]), label_y='pT (GeV)', weights=np.array([]),
             xmin=-1.,xmax=1., ymin=0, ymax=2*math.pi, nbins=[], title='', name='grid_CS'):
        plt.title(title)
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        plt.hist2d(data_x, data_y, bins=nbins, weights=weights, normed=True) 
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        plt.savefig(self.out_dir+'/'+name+'.png')


#################################################
