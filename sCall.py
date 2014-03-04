import numpy as np
import reissnerLib as SR
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11})


# Film parameters:
loading = 'pressure'  # Specify type of film loading
E = 1.32E6            # Modulus of elasticity, [Pa]
nu = 0.5              # Poisson's ratio
a = 7.5E-3            # Radius of film, [m]
dr = .001             # Discretization for nondimensional r
h = 18E-6             # Film thickness, [m]
quiet = 0


# Calculate prestrain due to thermal strain:
dT = 23                     # Temperature shift, [K]
CTE_pdms = 3.04E-4          # Coefficient of Thermal Expansion, [K^{-1}]
eps0 = dT * CTE_pdms        # Prestrain, []

# Calculate pressure due to cell loading, [Pa]:
h_cell = 3E-6                           # Height of cell, [m]
rho_cell = 1.05                          # Density of cell, [g/mL]
C = 1.0                                 # Percent confluency
p = SR.cellPressure(C,h_cell,rho_cell)  # Distributed pressure, [Pa]


# Build dictionary
params = {'loading':loading, 'E':E, 'a':a, 'nu':nu, 'dr':dr, 'h':h,'eps0':eps0,'quiet':quiet} 

if params['loading'] == 'pressure':
    params['p'] = p
elif params['loading'] == 'point':
    params['P'] = P


# Initialize object
plate = SR.reissner(params)
plate.solveReissner()

