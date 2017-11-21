import ROOT
import math
import sys
import numpy as np
import numpy.random as ran

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class TemplateMaker():

    def __init__(self):
        print "Initialize TemplateMaker"
        return
        
    def pdf(self, x, y, coeff=[]):
        return 3./16./np.pi * ((1.0 + x*x) + 0.5*coeff[0]*(1-3*x*x) + coeff[1]*2*x*np.sqrt(1-x*x)*np.cos(y) + 0.5*coeff[2]*(1-x*x)*np.cos(2*y) + coeff[3]*np.sqrt(1-x*x)*np.cos(y) + coeff[4]*x)
        #return 3./8.*(1.0 + x*x)

    def sample_BW(self, mass, width):
        x = max(min(np.random.standard_cauchy(),-3.0),+3.0)
        return (x*width + mass)

    def make_grid_CS(self, grid_px_cos=500, grid_px_phi=10, coeff=[0., 0., 0., 0., 0.], mass=80., width=2.0, ntoys=1000):

        self.ntoys=ntoys

        self.grid_px_cos = grid_px_cos
        self.grid_px_phi = grid_px_phi        

        # the uniform grid in cos and phi (CS)
        self.grid_cos = np.linspace(-1., 1.,      grid_px_cos, dtype=np.float)
        self.grid_phi = np.linspace( 0., 2*np.pi, grid_px_phi, dtype=np.float)
        self.grid_indexes = np.arange(0, len(self.grid_cos)*len(self.grid_phi), 1)

        # array that stores the rnd values
        self.rnd_grid = np.zeros((ntoys,3))

        xx, yy = np.meshgrid( self.grid_cos, self.grid_phi )        

        print "Evaluating the pdf over the grid..."
        pdf_evals = self.pdf(xx,yy, coeff=coeff).flatten()
        pdf_evals_norm = np.sum(pdf_evals)    
        pdf_evals *= (1./pdf_evals_norm if pdf_evals_norm>0. else 1.0)        

        # these are the values chosen by the pdf
        rnd = np.random.choice( self.grid_indexes, size=ntoys, p=pdf_evals )

        # fill the rnd_grid
        for i,idx in enumerate(rnd):
            idx_cos =  idx%self.grid_px_cos
            idx_phi =  idx/self.grid_px_cos
            # cos
            self.rnd_grid[i][0] = self.grid_cos[idx_cos]
            # phi
            self.rnd_grid[i][1] = self.grid_phi[idx_phi]
            # E
            self.rnd_grid[i][2] = self.sample_BW(mass, width)*0.5

        #print self.rnd_grid

        
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
        

    def make_grid_lab(self, y=0.0, pt=0.0):

        print ("Make the grid in the lab for (pt,y)=(%s,%s)" % (pt,y))
        self.rnd_grid_lab = np.zeros((self.ntoys,2))

        for itoy in range(self.ntoys):
            toy_CS = self.rnd_grid[itoy]
            x_lab = self.boost_to_lab( toy_CS, y=y, pt=pt)
                        
            # pt_lab
            self.rnd_grid_lab[itoy][0] = math.sqrt(x_lab[0]*x_lab[0] - x_lab[3]*x_lab[3]) if x_lab[0]>abs( x_lab[3] ) else 0.
            # eta_lab
            self.rnd_grid_lab[itoy][1] = 0.5*math.log((1+x_lab[3]/x_lab[0])/(1-x_lab[3]/x_lab[0])) if  x_lab[0]>abs( x_lab[3] ) else 10.
            #if  self.rnd_grid_lab[itoy][0]>50.:
            #    print toy_CS, x_lab, "pt,eta", self.rnd_grid_lab[itoy][0], self.rnd_grid_lab[itoy][1] 

        print "Done!"

    def plot(self, data_x=np.array([]), data_y=np.array([]), nbins=10, name='grid_CS'):
        plt.hist2d(data_x, data_y, bins=nbins)
        plt.savefig(name+'.png')



##############################

template = TemplateMaker()
template.make_grid_CS( grid_px_cos=200, grid_px_phi=200, coeff=[0.,0.,0.,0.,0.], mass=80., width=2.0, ntoys=100000 )
template.plot(data_x=getattr(template,'rnd_grid')[:,0], data_y=getattr(template,'rnd_grid')[:,1],nbins=30, name='grid_CS')

#template.make_grid_lab( y=1.0, pt=0.0 )
#template.plot(data_y=getattr(template,'rnd_grid_lab')[:,0], data_x=getattr(template,'rnd_grid_lab')[:,1], nbins=100, name='grid_lab')

#template.make_grid_lab( y=1.0, pt=2.0 )
#template.plot(data_y=getattr(template,'rnd_grid_lab')[:,0], data_x=getattr(template,'rnd_grid_lab')[:,1], nbins=100, name='grid_lab_pt2')

#template.make_grid_lab( y=1.0, pt=4.0 )
#template.plot(data_y=getattr(template,'rnd_grid_lab')[:,0], data_x=getattr(template,'rnd_grid_lab')[:,1], nbins=100, name='grid_lab_pt4')

template.make_grid_lab( y=1.0, pt=10.0 )
template.plot(data_y=getattr(template,'rnd_grid_lab')[:,0], data_x=getattr(template,'rnd_grid_lab')[:,1], nbins=100, name='grid_lab_pt10')
