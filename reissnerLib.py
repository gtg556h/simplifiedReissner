import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
from scipy.special import i1
from scipy.special import k1
plt.rcParams.update({'font.size':11})


# State: Reasonably comfortable that I coded the solution correctly
# However, still a SIGNIFICANT dispartity between direct formula for 'beta' and the relation between g and beta.  There is a significant variation in shape, not just magnitude, so I suspect that the culprit may be a typo in the argument of the Bessel functions.


class reissner(object):

    def __init__(self, params):
        print('Initializing system')

        # Unpack dictionary
        self.loading = params['loading']
        self.E = params['E']
        self.a = params['a']
        self.nu = params['nu']
        self.dr = params['dr']
        self.h = params['h']
        self.eps0 = params['eps0']

        if self.loading == 'pressure':
            self.p = params['p']             # pressure, [Pa]
            self.N = self.p*self.a**2/2      # Stress resultant, [N]
        elif self.loading == 'point':
            self.P = params['P']             # Point load, [N]
            self.N = self.P/(2*np.pi)        # Stress resultant, [N]
   
        # Basic parameters:
        self.D = (self.E*self.h**3)/(12*(1-self.nu**2))     # Bending stiffness
        self.A = 1/(self.E*self.h)                     # Stretching compliance
        self.k = 1/(1-self.nu)
        self.eps = ((self.A*self.D)/self.a**2)**(0.5)
        self.sigma = self.A*self.N/self.a

        # Nondimensionalization
        self.alpha = np.log(self.eps0)/np.log(self.eps)
        self.gamma = np.log(self.sigma)/np.log(self.eps)
        self.r = np.arange(self.dr,1,self.dr)
        if self.alpha > 2 and self.gamma > 3:
            print('Region 1')
            self.region = 1
        elif self.alpha < 2 and self.gamma > 3*self.alpha/2.0:
            print('Region 2')
            self.region = 2
        elif self.gamma < 3 and self.alpha > 2*self.gamma/3.0:
            print('Region 3')
            self.region = 3
        else:
            print('My mind is blown')


    def solveReissner(self):
        if self.region == 2 and self.loading == 'pressure':
            self.region2_pressure()
        else:
            print('Code more solutions you bum')


    def region2_pressure(self):
        self.kappa = np.sqrt(12*self.eps0*(1+self.nu)*(self.a/self.h)**2)
        self.delta = self.gamma - self. alpha
        self.g = -self.eps0/self.eps**2/self.kappa**2*(self.r-i1(self.kappa*self.r)/i1(self.kappa))
        self.beta = -6*(1-self.nu**2)/self.kappa**2*(self.p*self.a**3)/(self.E*self.h**3)*(self.r-i1(self.k*self.r)/i1(self.k))
        self.numIntegrateBeta()

        # Making ends meet: compare caclulated beta from relation to g with direct formula.  Is this appropriate??
        if 1:
            self.beta = self.eps**self.delta*self.g
            self.numIntegrateBeta()


######## VALIDATE THIS SCRIPT!!!!
    def numIntegrateBeta(self):
        self.w = np.zeros_like(self.beta)
        for i in range(0,self.beta.shape[0]):
            self.w[i] = self.dr*np.sum(np.sin(self.beta[i:self.r.shape[0]-1]))
            
        self.maxDisp = np.max(self.w) - np.min(self.w)


###################################################
###################################################

def cellPressure(C, h, rho):
    p = C * h * (rho-1.00)*9810
    return p
