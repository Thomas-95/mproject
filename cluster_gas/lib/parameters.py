import numpy as np

# -----------------------------------------------------------------------------
# Module for calculating the condensation & evaporation rates, given the
# parameters with which the system is set up. 
# FIXME Most of these parameters are yet to be properly calculated. 
# -----------------------------------------------------------------------------

# Define some constants.
D = 2.32e-9            # The mass diffusion coefficient (298.16K, 1atm).
mon_vol = 2.992e-29    # The monomer volume, taken as vol of one water molecule.
k_B = 1.38e-23         # The Boltzmann constant. 
T = 273.15             # The temperature of the system. 
jump_freq = 1.         # The jump frequency.
d_jump = 1.            # The jump distance for monomer to attach to cluster.
C_0 = 1.               # All sites are available for homogeneous nucleation. 
C_1 = 1.               # Initial condition for C_1, conc. of monomers.


def beta(n, C_1=C_1):  # n is the class of clusters. 

    R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)
    kappa = D/(jump_freq*d_jump)

    return 4*np.pi*(R_n**2/(R_n + kappa))*(D/mon_vol)*C_1

def sigma(n):          # Surface free energy in generalised capillary approx.

    sigma_const = (28.0 + (T-273.15)/4.)*1e-4      # -36 < T < 0 Celcius.
    n_0 = ((32e-30*np.pi)/(3*mon_vol))**(1./3.)

    return sigma_const*(1 + (n_0/float(n))**(1./3.))**-2

def alpha(n):  # Calculate evaporation rate from (n+1) clusters.

    exponent = (36*np.pi*mon_vol**2)**(1./3.) * \
               (((n+1)**(2./3.)*sigma(n+1)    - \
               n**(2./3.)*sigma(n) - sigma(1))) / (k_B*T)

    return beta(n, 1.)*C_0*np.exp(exponent)
