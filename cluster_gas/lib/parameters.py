import numpy as np

# -----------------------------------------------------------------------------
# Module for calculating the condensation & evaporation rates, given the
# parameters with which the system is set up. 
# -----------------------------------------------------------------------------

# Define some constants.
D = 2.32e-9            # Mass diffusion coefficient (298.16K, 1atm).
mon_vol = 2.992e-29    # Monomer volume, taken as vol of one water molecule.
k_B = 1.38e-23         # Boltzmann constant. 
T = 272.85             # Temperature of the system. 
jump_freq = 1.e13      # Jump freq. taken as vibration freq of water molecule.
d_jump = 2.5e-10       # Jump distance for monomer to attach to cluster.
                       # Taken as length of Hydrogen bond.
#jump_freq = 1.        # COULD THESE BE FITTED EMPIRICALLY?
#d_jump = 1.
C_0 = 1.               # All sites are available for homogeneous nucleation.
sigma_const = (28.0 + (T-273.15)/4.)*1e-4      # -36 < T < 0 Celcius.


def beta(n, C_1):  # n is the class of clusters. 

    R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)
    kappa = D/(jump_freq*d_jump)

    return 4*np.pi*(R_n**2/(R_n + kappa))*(D/mon_vol)*C_1

def sigma(n):          # Surface free energy in generalised capillary approx.

    n_0 = ((32e-30*np.pi)/(3*mon_vol))**(1./3.)    # Tolman, 1949.

    return sigma_const*(1 + (n_0/float(n))**(1./3.))**-2

def alpha(n):  # Calculate evaporation rate from (n+1) clusters.

    exponent = (36*np.pi*mon_vol**2)**(1./3.) * \
               (((n+1)**(2./3.)*sigma(n+1)    - \
               n**(2./3.)*sigma(n) - sigma(1))) / (k_B*T)

    return beta(n, 1.)*C_0*np.exp(exponent)

def critical_radius(T=T):
    assert(T < 273.15), "System not supercooled: please lower the temperature."

    L = 333.55
    mon_mass = 2.992e-23  # mass of water molecule
    T_m = 273.15
    G_nuc = L*mon_mass*((T-T_m)/T_m)

    return -(2./3.)*(36*np.pi*mon_vol**2)**(1./3.)*(sigma_const/G_nuc)

# Implement a binary search algorithm to see where alpha = beta to check above.
