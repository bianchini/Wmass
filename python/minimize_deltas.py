import ROOT
import math
from array import array;

# class that fits for the shifts
class Minimze_deltas:

    def __init__(self, deltas=[], nbins=1, step=0.01, fix_params=[]):
        
        self.nbins = nbins
        self.deltas = deltas
        self.nmu = len(self.deltas)-self.nbins

        # initial values for the deltas
        mu_start = array( 'd' )
        mu_step = array( 'd' )

        mu_start.append( deltas[0] )
        mu_start.append( deltas[1] )
        for i in range(nbins):
            mu_start.append( deltas[2+i] )
        for i in range(nbins):
            mu_start.append( deltas[2+2*nbins+i] )
        for i in range(nbins):
            mu_start.append( deltas[2+3*nbins+i] )
        mu_start.append( deltas[2+4*nbins] )
        mu_start.append( deltas[2+4*nbins+1] )

        for i in range(self.nmu):  
            mu_step.append( step )

        self.gMinuit = ROOT.TMinuit(50)  
        self.gMinuit.SetFCN( self.fcn ) 

        self.arglist = array( 'd', 10*[0.] )
        self.ierflg = ROOT.Long(0)
        self.arglist[0] = 0.5
        self.gMinuit.mnexcm( "SET ERR", self.arglist, 1, self.ierflg )

        to_be_fixed = array('i')
        to_be_fixed.extend( fix_params )
        #to_be_fixed.extend( [1, 2+3*self.nbins] )
        #to_be_fixed.extend( range(self.nbins+2, 2*self.nbins+2) )        
        #to_be_fixed = array('i')

        for p in range(self.nmu):
            p_start = max(mu_start[p], 1.0) if p not in to_be_fixed else 0.0
            p_low = -1.0
            p_high = mu_start[p] + max(5.*math.sqrt( mu_start[p] ), 20. ) 
            self.gMinuit.mnparm( p, "par_"+str(p), p_start, mu_step[p], p_low,  p_high, self.ierflg )
            if p in to_be_fixed:
                self.gMinuit.FixParameter(p)




    def run(self, n_points=50000):
        
        delta_n_u = array('f', self.nbins*[0.] )
        delta_n_d = array('f', self.nbins*[0.] )
        delta_m_u = array('f', self.nbins*[0.] )
        delta_m_d = array('f', self.nbins*[0.] )
        var_delta_m = array('f', self.nbins*[0.] )

        self.arglist[0] = n_points
        self.arglist[1] = 1.
        self.gMinuit.mnexcm( "MIGRAD", self.arglist, 2, self.ierflg )

        # Print results
        amin, edm, errdef = ROOT.Double(0.), ROOT.Double(0.), ROOT.Double(0.)
        nvpar, nparx, icstat = ROOT.Long(1983), ROOT.Long(1984), ROOT.Long(1985)
        self.gMinuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat )
        self.gMinuit.mnprin( 3, amin )

        self.cov = ROOT.TMatrixDSym(self.nmu)
        self.gMinuit.mnemat(self.cov.GetMatrixArray(), self.nmu)

        self.res_val = array('d', self.nmu*[0.])
        self.res_err = array('d', self.nmu*[0.])
        for mu in range(self.nmu):
            val = ROOT.Double(0.)
            err = ROOT.Double(0.)
            self.gMinuit.GetParameter(mu, val, err)
            self.res_val[mu] = val
            self.res_err[mu] = err

        for i in range(self.nbins):
            delta_n_u[i] = self.deltas[2+i-1 if i>0 else 0] + self.deltas[2+2*self.nbins+i+1 if i<self.nbins-1 else 2+4*self.nbins] - self.deltas[2+i] - self.deltas[2+2*self.nbins+i] 
            delta_n_d[i] = self.deltas[2+self.nbins+i-1 if i>0 else 1] + self.deltas[2+3*self.nbins+i+1 if i<self.nbins-1 else 2+4*self.nbins+1] - self.deltas[2+self.nbins+i] - self.deltas[2+3*self.nbins+i] 
            delta_m_u[i] = self.res_val[2+i-1 if i>0 else 0 ] + self.res_val[2+self.nbins+i+1 if i<self.nbins-1  else 2+3*self.nbins] - self.res_val[ 2+i ] - self.res_val[ 2+self.nbins+i ]
            delta_m_d[i] = -delta_m_u[i]
            var_delta_m[i] = 0.
            indexes = [2+i-1 if i>0 else 0, 2+self.nbins+i+1 if i<self.nbins-1  else 2+3*self.nbins, 2+i,  2+self.nbins+i]
            for ind in range(4):
                add = math.pow(self.res_err[indexes[ind]], 2.)
                var_delta_m[i] += add
            for ind1 in range(4):
                for ind2 in range(4):
                    if ind1==ind2: 
                        continue
                    add = self.res_err[indexes[ind1]]*self.res_err[indexes[ind2]]*self.cov(indexes[ind1],indexes[ind2]) * (-1 if (ind1>1 and ind2<2) or (ind1>1 and ind2<2) else +1)
                    var_delta_m[i] += 2*add
    
        self.gMinuit.IsA().Destructor(self.gMinuit)
        return [delta_n_u, delta_n_d, delta_m_u, delta_m_d, var_delta_m]

    def lf(self, m, n ):
        value = 0.
        if n>0:
            value = (m - n*math.log(m)) if m >0. else 1e+10
        else:
            #value = abs(math.pow(m,2.))
            value = math.pow(abs(m),2.0)
        return value

 
    def fcn(self, npar, gin, f, par, iflag ):

        nll = 0.
        nll += self.lf(par[0], self.deltas[0])
        nll += self.lf(par[1], self.deltas[1])
        
        for i in range(self.nbins):      
            nll += self.lf(par[2+i], self.deltas[2+i])
            mu = -par[2+i] + par[2+self.nbins+(i+1) if i<self.nbins-1 else 2+3*self.nbins] + par[2+2*self.nbins+(i+1) if i<self.nbins-1 else 2+3*self.nbins+1] - par[2+self.nbins] - par[2+2*self.nbins] + par[0] + par[1]
            nll += self.lf( mu , self.deltas[2+self.nbins+i])
            nll += self.lf(par[2+self.nbins+i], self.deltas[2+2*self.nbins+i])
            nll += self.lf(par[2+2*self.nbins+i], self.deltas[2+3*self.nbins+i])

        nll += self.lf(par[2+3*self.nbins], self.deltas[2+4*self.nbins])
        nll += self.lf(par[2+3*self.nbins+1], self.deltas[2+4*self.nbins+1])
        f[0] = nll


################################

