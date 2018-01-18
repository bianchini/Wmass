import ROOT
import numpy as np
import math
from array import array;

np.random.seed(0)

class Bias:
    def __init__(self, mu=100, prior=0.1):
        self.mu = mu
        self.mu_tmp = mu
        self.prior = prior
        self.gMinuit = ROOT.TMinuit(10)
        self.gMinuit.SetPrintLevel(-1)
        self.gMinuit.SetFCN( self.fcn ) 
        self.arglist = array( 'd', 10*[0.] )
        self.ierflg = ROOT.Long(0)
        self.arglist[0] = 1.0
        self.gMinuit.mnexcm( "SET ERR", self.arglist, 1, self.ierflg )
        self.gMinuit.DefineParameter( 0, "par", mu-10, 0.1, mu/2, mu*2)
        return

    def get_random(self):
        self.n = np.random.poisson(self.mu_tmp, 1)

    def sample_mu(self):
        a = -1
        while a<0:
            a = np.random.normal(self.mu, self.prior)
        self.mu_tmp = a

    def fcn(self, npar, gin, f, par, iflag ):
        nll = 0.0
        nll += math.pow( (par[0] - self.mu_tmp)/(self.prior*self.mu_tmp), 2.0)
        nll += 2*(par[0] - self.n*math.log(par[0]))        
        #nll += math.pow( (par[0] - self.n)/(math.sqrt(self.mu_tmp)), 2.0)
        f[0] = nll

    def run(self, ntoys=100):
        
        pull = 0.0
        ratio = 0.0
        for i in range(ntoys):
            self.sample_mu()
            self.get_random()
            self.arglist[0] = 50000
            self.arglist[1] = 1.
            self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )
            val = ROOT.Double(0.)
            err = ROOT.Double(0.)
            self.gMinuit.GetParameter(0, val, err)
            pull += (val-self.mu_tmp)/err
            ratio += val/self.mu_tmp
            self.arglist[0] = 1
            self.arglist[1] = self.mu
            self.gMinuit.mnexcm("SET PAR", self.arglist, 2, self.ierflg )
        
        print('%s +/- %s' % (pull/ntoys, 1./math.sqrt(ntoys)))
        print('%s +/- %s' % (ratio/ntoys, 1./math.sqrt(ntoys)))

######################
bias = Bias(mu=10, prior=0.1)
bias.run(10000)
