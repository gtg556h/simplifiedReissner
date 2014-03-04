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
#C = 1.0                                 # Percent confluency
dT = 10                     # Temperature shift, [K]
CTE_pdms = 3.04E-4          # Coefficient of Thermal Expansion, [K^{-1}]
eps0 = dT * CTE_pdms        # Prestrain, []



# Build dictionary
params = {'loading':loading, 'E':E, 'a':a, 'nu':nu, 'dr':dr, 'h':h,'quiet':quiet,'eps0':eps0}#,'eps0':eps0} 



# Main Loop:
loadVec = np.arange(0,1,.01)
displacement = np.zeros(loadVec.shape[0])
regionVec = np.zeros(loadVec.shape[0])
for i in range(0,loadVec.shape[0]):
    if np.mod(i,10)==0:
        print(i/loadVec.shape[0])
# Calculate prestrain due to thermal strain:
    C = loadVec[i]
    p = SR.cellPressure(C,h_cell, rho_cell)
    params['p']=p

    # Initialize object
    plate = SR.reissner(params)
    plate.solveReissner()
    displacement[i] = plate.maxDisp
    regionVec[i] = plate.region

plt.plot(loadVec,displacement*10**6)
plt.show()

