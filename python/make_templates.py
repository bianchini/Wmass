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
        self.save_plots = True
        self.width = 2.0
        return


    def set_CS_grid_mesh(self, grid_px_cos=500, grid_px_phi=10):
        self.grid_px_cos = grid_px_cos
        self.grid_px_phi = grid_px_phi        

        # the uniform grid in cos and phi (CS)
        self.grid_cos = np.linspace(-1., 1.,      grid_px_cos, dtype=np.float)
        self.grid_phi = np.linspace( 0., 2*np.pi, grid_px_phi, dtype=np.float)
        self.grid_indexes = np.arange(0, len(self.grid_cos)*len(self.grid_phi), 1)
        
    def pdf(self, x, y, coeff=[]):
        return 3./16./np.pi * ((1.0 + x*x) + 0.5*coeff[0]*(1-3*x*x) + coeff[1]*2*x*np.sqrt(1-x*x)*np.cos(y) + 0.5*coeff[2]*(1-x*x)*np.cos(2*y) + coeff[3]*np.sqrt(1-x*x)*np.cos(y) + coeff[4]*x)
        #return 3./8.*(1.0 + x*x)

    def BW(self, x, mass):
        bw = 1.0/math.pi/self.width/(1+math.pow((x-mass)/self.width,2.0))
        return bw

    def sample_BW(self, mass, truncate=True):
        x = np.random.standard_cauchy()
        check = lambda x : (x<3.0 and x>-3.0) 
        while truncate and not check(x): 
            x = np.random.standard_cauchy()
        return (x*self.width + mass)

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


    def make_grid_CS(self, coeff=[0., 0., 0., 0., 0.], mass=80., ntoys=1000):

        self.ntoys=ntoys

        # array that stores the rnd values
        rnd_grid = np.zeros((ntoys,3))

        xx, yy = np.meshgrid( self.grid_cos, self.grid_phi )        

        print "Evaluating the pdf over the grid..."
        pdf_evals = self.pdf(xx,yy, coeff=coeff).flatten()
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

    def make_grid_lab(self, boost_vecs=[], weights_mass=[], coeff=[0.,0.,0.,0.,0.], mass=80., ntoys=100000):

        grid_tag = 'M{:05.3f}'.format(mass)
        for ic,c in enumerate(coeff):
            grid_tag += '_A'+str(ic)+('{:03.2f}'.format(c))

        rnd_grid_CS = self.make_grid_CS( coeff=coeff, mass=mass, ntoys=ntoys )
        np.save('../test/grid_CS_'+grid_tag, rnd_grid_CS)

        if self.save_plots:
            self.plot(data_x=rnd_grid_CS[:,0], label_x='cos(theta)',              
                      data_y=rnd_grid_CS[:,1], label_y='phi', nbins=30, 
                      weights=np.ones(ntoys),
                      title=('CS frame, M=%s, A=[%s,%s,%s,%s,%s]' % (mass, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4] )), 
                      name='grid_CS_'+grid_tag)

        
        for boost_vec in boost_vecs:
            pt, y = boost_vec
            print ("Make the grid in the lab for (pt,y)=(%s,%s)" % (pt,y))
            
            boost_tag = grid_tag+'_pt{:02.1f}'.format(pt)+'_'+'y{:02.1f}'.format(y)

            # fill the template
            setattr(self, 'rnd_grid_lab_'+boost_tag,  np.zeros((self.ntoys, 2 + len(weights_mass))) )
            for itoy in range(self.ntoys):
                toy_CS = rnd_grid_CS[itoy]
                x_lab = self.boost_to_lab( toy_CS, y=y, pt=pt)                        
                this_toy = getattr(self, 'rnd_grid_lab_'+boost_tag)[itoy] 
                this_toy[0] = math.sqrt(x_lab[0]*x_lab[0] - x_lab[3]*x_lab[3]) if x_lab[0]>abs( x_lab[3] ) else 0. # pt_lab
                this_toy[1] = 0.5*math.log((1+x_lab[3]/x_lab[0])/(1-x_lab[3]/x_lab[0])) if  x_lab[0]>abs( x_lab[3] ) else 10. # eta lab
                for im,m in enumerate(weights_mass):
                    weight = self.BW(2.0*toy_CS[2], m)/self.BW(2.0*toy_CS[2], mass)
                    #print mass, 2.0*toy_CS[2], " ==> " , m, weight
                    this_toy[2+im] = weight
            
            for im,m in enumerate(weights_mass):
                template = np.histogram2d( getattr(self, 'rnd_grid_lab_'+boost_tag)[:,1], getattr(self, 'rnd_grid_lab_'+boost_tag)[:,0],
                                           weights=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,2+im],
                                           normed=True,
                                           range=[[-2.5, +2.5], [15.0, 60.0]],
                                           bins=[30,30])
                # save to disk
                np.save('../test/grid_lab_'+boost_tag+'_weighted{:05.3f}'.format(m), template)

                if self.save_plots:            
                    self.plot(data_x=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,1], label_x='eta',              
                              data_y=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,0], label_y='pt', nbins=30, 
                              weights=getattr(self, 'rnd_grid_lab_'+boost_tag)[:,2+im],
                              title=( 'Lab frame, M=%s, A=[%s,%s,%s,%s,%s] (pt,y)=(%s,%s)' % (m, coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], pt, y)), 
                              xmin=-3.0,xmax=3.0, ymin=20, ymax=60,
                              name='grid_lab_'+boost_tag+'_weighted{:05.3f}'.format(m))

    def plot(self, data_x=np.array([]), label_x='eta', data_y=np.array([]), label_y='pT (GeV)', weights=np.array([]),
             xmin=-1.,xmax=1., ymin=0, ymax=2*math.pi, nbins=10, title='', name='grid_CS'):
        plt.title(title)
        plt.xlabel(label_x)
        plt.ylabel(label_y)
        plt.hist2d(data_x, data_y, bins=nbins, weights=weights, normed=True) 
        #plt.xlim((xmin,xmax))
        #if ymax>0.:
        #    plt.ylim((ymin,ymax))
        plt.savefig('../test/'+name+'.png')



##############################

boost_vecs=[]
weights_mass=[80.000]
for pt in [0.0, 5.0, 10.0, 20.0]:
    for y in [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]:
        boost_vecs.append([pt,y])
coeff=[0.,0.,0.,0.,+1.0]

template = TemplateMaker()
template.set_CS_grid_mesh(grid_px_cos=200, grid_px_phi=200)
template.make_grid_lab( boost_vecs=boost_vecs, weights_mass=weights_mass, coeff=coeff, mass=80.000, ntoys=10000 )
