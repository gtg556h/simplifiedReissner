import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg
from scipy.special import i1
from scipy.special import k1
plt.rcParams.update({'font.size':11})



class reissner(object):

    def __init__(self, params):
        print('Initializing system')

        # Unpack dictionary
        self.loading = params['loading']
        if self.loading == 'pressure':
            self.p = params['p']
        elif self.loading == 'point':
            self.P = params['P']
        self.E = params['E']
        self.a = params['a']
        self.nu = params['nu']
        self.dr = params['dr']
        self.h = params['h']
            
        # Basic parameters:
        self.D = (E*h**3)/(12*(1-nu**2))     # Bending stiffness
        self.A = 1/(E*h)                     # Stretching compliance
        self.k = 1/(1-nu)

        N = p*a**2/2    # WHAT THE FUCK IS THIS???
        
        # Nondimensionalization:

    def function(self):
        print('blah')




###################################################
###################################################

def cellPressure(C, h, rho):
    p = C * h_cell * (rho_cell-1.00)*9810
