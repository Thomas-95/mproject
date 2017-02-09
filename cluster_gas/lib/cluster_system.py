import numpy as np, ODE_int

# To be done:

# Could jump_freq, d_jump be fitted empirically?
# Implement a binary search algorithm to see where alpha = beta.

'''def critical_radius(T=T):
    assert(T < 273.15), "System not supercooled: please lower the temperature."

    L = 333.55
    mon_mass = 2.992e-23  # mass of water molecule
    T_m = 273.15
    G_nuc = L*mon_mass*((T-T_m)/T_m)

    return -(2./3.)*(36*np.pi*mon_vol**2)**(1./3.)*(sigma_const/G_nuc)'''


class ClusterSystem():

    k_B = 1.38e-23                  # Boltzmann constant. 
    
    def __init__(self, n_class, T, D=2.32e-9, mon_vol=2.992e-29, \
                       d_jump=2.5e-10, jump_freq=1.e13):
                       
        self.n_class = n_class      # Number of cluster sizes to sample.
        self.T = T                  # Temperature of the system. 
        self.D = D                  # Mass diffusion coeff't (298.16K, 1atm).
        self.mon_vol = mon_vol      # Monomer volume = one water molecule vol.
        self.C_0 = 1.               # All sites poss for homogeneous nucleation.
        self.d_jump = d_jump        # Jump distance for mon to attach w/cluster.
        self.jump_freq = jump_freq  # Taken as vibration freq of water molecule.
        self.sigma_const = (28.0 + (T-273.15)/4.)*1e-4    # -36 < T < 0 Celcius.
        self.C_init = np.zeros(shape=(n_class+1))
        self.C_init[0] = 1.0        # Set initial condition - all monomers.
        
        
    def beta(self, n, C_1):  # n = the class of clusters, C_1 is monomer conc.

        R_n = ((3*n*self.mon_vol)/(4*np.pi))**(1./3.)
        kappa = D/(self.jump_freq*self.d_jump)

        return 4*np.pi*(R_n**2/(R_n + kappa))*(self.D/self.mon_vol)*C_1  
      
      
    def sigma(self, n):      # Surface free energy in generalised capillary approx.

        n_0 = ((32e-30*np.pi)/(3*mon_vol))**(1./3.)    # Tolman, 1949.

        return sigma_const*(1 + (n_0/float(n))**(1./3.))**-2
          
          
    def alpha(self, n):      # Calculate evaporation rate from (n+1) clusters.

        exponent = (36*np.pi*mon_vol**2)**(1./3.) * \
                   (((n+1)**(2./3.)*sigma(n+1)    - \
                   n**(2./3.)*sigma(n) - sigma(1))) / (k_B*T)
               
        kappa = D/(jump_freq*d_jump)
        R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)

        return 4*np.pi*(R_n**2/(R_n + kappa))*(D/mon_vol)*C_0*np.exp(exponent)
    
    
    @staticmethod
    def tridiag(a, b, c, k1=-1, k2=0, k3=1):
        return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


    def generate__update_matrix(self):
        '''Function generating the matrix that describes our linear system. 
        Work within the small cluster regime, for now.'''
    
        diags = [0 for i in range(self.n_class)]
        lower_diags = [self.beta(i+1, C_1) for i in range(self.n_class-1)]
        upper_diags = [self.alpha(i+1) for i in range(self.n_class-1)]

        diags[-1] = -self.alpha(self.n_class)
        for i in range(1, n_class-1):
            diags[i] = -(self.beta(i+1, C_1) + self.alpha(i+1))

        M = self.tridiag(lower_diags, diags, upper_diags)

        M[0][1:-1] += self.alpha(i) - self.beta(i, C_1)
        M[0,0]      = -2*self.beta(1, C_1)
        M[0, -1]    = self.alpha(self.n_class)
    
        return M
        
        
    @staticmethod
    def calc_mass(x):
        return sum((index+1)*value for index, value in enumerate(x))


    def solve_system(self, N_ITER):
        soln = [C_init]
        times = np.linspace(0, 0.8e-12, N_ITER)
    
        for t in times:
            M = self.generate_update_matrix(n_class=n_class, C_1=soln[-1][0])
            x_new = ODE_int.RK4(M, soln[-1][0:-1], t, 1)[0]
            #x_new[0] += (1-self.calc_mass(x_new))  # Add line to conserve mass?
            soln.append(x_new)

        return np.asarray(soln)
