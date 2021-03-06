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
quiet = 1


# Calculate pressure due to cell loading, [Pa]:
h_cell = 3E-6                           # Height of cell, [m]
rho_cell = 1.05                          # Density of cell, [g/mL]
C = 1.0                                 # Percent confluency
p = SR.cellPressure(C,h_cell,rho_cell)  # Distributed pressure, [Pa]


# Build dictionary
params = {'loading':loading, 'E':E, 'a':a, 'nu':nu, 'dr':dr, 'h':h,'quiet':quiet}#,'eps0':eps0} 

if params['loading'] == 'pressure':
    params['p'] = p
elif params['loading'] == 'point':
    params['P'] = P


# Main Loop:
tempVec = np.arange(1,30,.25)
displacement = np.zeros(tempVec.shape[0])
regionVec = np.zeros(tempVec.shape[0])
for i in range(0,tempVec.shape[0]):
    if np.mod(i,10)==0:
        print(i/tempVec.shape[0])
# Calculate prestrain due to thermal strain:
    dT = tempVec[i]                     # Temperature shift, [K]
    CTE_pdms = 3.04E-4          # Coefficient of Thermal Expansion, [K^{-1}]
    eps0 = dT * CTE_pdms        # Prestrain, []
    params['eps0'] = eps0


    # Initialize object
    plate = SR.reissner(params)
    plate.solveReissner()
    displacement[i] = plate.maxDisp
    regionVec[i] = plate.region

plt.plot(tempVec,displacement)
plt.show()

