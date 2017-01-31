import numpy as np

# Define some constants.
D = 1.          # The diffusion coefficient. 
omega = 1.      # The monomer volume.
gamma = 1.      # The jump frequency.
T = 273.        # The temperature of the system. 
d_jump = 1.     # The jump distance for monomer to attach to cluster. 
k_B = 1.38e-23  # The Boltzmann constant.
#C_0 = N_0/N_s   # Fraction of accessible sites. 
C_0 = 1.        # Temporary value. 

# -----------------------------------------------------------------------------
# Function to calculate the condensation & evaporation rates.
# Will return to this later, once some sort of matrix system is set up.
# -----------------------------------------------------------------------------

def beta(n, C_1):  # C_1 is monomer concentration, n is the class of clusters. 

    R_n = ((3*n*omega)/(4*np.pi))**(1./3.)
    kappa = D/(gamma*d_jump)

    return 4*np.pi*(R_n**2/(R_n + kappa))*(D/omega)*C_1

def sigma(n):   # Calculate surface free energy barrier. 

    return 1.

def alpha(n):  # Calculate evaporation rate from (n+1) clusters.

    exponent = (36*np.pi*omega**2)**(1./3.) * \
               (((n+1)**(2./3.)*sigma(n+1)  - \
               n**(2./3.)*sigma(n) - sigma(1))) / (k_B*T)

    #return beta(n, 1.)*C_0*np.exp(exponent)
    return 6.

#n_class = 6
#M = update_matrix(n_class=n_class)
