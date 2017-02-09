import numpy as np

# -----------------------------------------------------------------------------
# Module for calculating the condensation & evaporation rates, given the
# parameters with which the system is set up. 
# -----------------------------------------------------------------------------

# Could jump_freq, d_jump be fitted empirically?

sigma_const = (28.0 + (T-273.15)/4.)*1e-4      # -36 < T < 0 Celcius.


def critical_radius(T=T):
    assert(T < 273.15), "System not supercooled: please lower the temperature."

    L = 333.55
    mon_mass = 2.992e-23  # mass of water molecule
    T_m = 273.15
    G_nuc = L*mon_mass*((T-T_m)/T_m)

    return -(2./3.)*(36*np.pi*mon_vol**2)**(1./3.)*(sigma_const/G_nuc)

# Implement a binary search algorithm to see where alpha = beta to check above.


class ClusterSystem():

    k_B = 1.38e-23                  # Boltzmann constant. 
    
    def __init__(self, n_class, T, D=2.32e-9, mon_vol=2.992e-29, \
                       d_jump=2.5e-10, jump_freq=1.e13):
                       
        self.n_class = n_class      # Number of cluster sizes to sample.
        self.T = T                  # Temperature of the system. 
        self.D = D                  # Mass diffusion coeff't (298.16K, 1atm).
        self.mon_vol = mon_vol      # Monomer volume = one water molecule vol.
        self.C_0 = 1.               # All sites poss for homogeneous nucleation.
        self.d_jump = d_jump        # Jump distance for mon to attach to cluster.
        self.jump_freq = jump_freq  # Taken as vibration freq of water molecule.
        
        
    @staticmethod
    def beta(n, C_1):  # n = the class of clusters, C_1 is monomer conc.

        R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)
        kappa = D/(jump_freq*d_jump)

    return 4*np.pi*(R_n**2/(R_n + kappa))*(D/mon_vol)*C_1  
    
          
    @staticmethod      
    def sigma(n):      # Surface free energy in generalised capillary approx.

        n_0 = ((32e-30*np.pi)/(3*mon_vol))**(1./3.)    # Tolman, 1949.

    return sigma_const*(1 + (n_0/float(n))**(1./3.))**-2
          
          
    @staticmethod
    def alpha(n):      # Calculate evaporation rate from (n+1) clusters.

        exponent = (36*np.pi*mon_vol**2)**(1./3.) * \
                   (((n+1)**(2./3.)*sigma(n+1)    - \
                   n**(2./3.)*sigma(n) - sigma(1))) / (k_B*T)
               
        kappa = D/(jump_freq*d_jump)
        R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)

    return 4*np.pi*(R_n**2/(R_n + kappa))*(D/mon_vol)*C_0*np.exp(exponent)
    
    
    

