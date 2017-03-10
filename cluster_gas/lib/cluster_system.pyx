import numpy as np, ODE_int
cimport numpy as np


cdef class ClusterSystem(object):
    
    cdef int n_class, soln
    cdef double T, D, mon_vol, sigma_const,  mon_mass, k_B
    cdef np.ndarray C_init
    #cdef double C_init[]
    
    def __init__(self, nclass, temp):
        self.n_class = nclass      # Number of cluster sizes to sample.
        self.T = temp                 # Temperature of the system. 
        self.D = 2.32e-9                  # Mass diffusion coeff't (298.16K, 1atm).
        self.mon_vol = 2.992e-29      # Monomer volume = one water molecule vol.
        self.sigma_const = (28.0 + (self.T-273.15)/4.)*1e-4    # -36 < T < 0 Celcius.
        self.C_init = np.zeros(shape=(self.n_class+1))
        self.C_init[0] = 3.35      # Set initial condition - all monomers.
        self.mon_mass = 2.992e-26
        self.k_B = 1.38e-23                  # Boltzmann constant.
        
        
    cdef double R_n(self, n):
        return ((3*n*self.mon_vol)/(4*np.pi))**(1./3.)

    cdef double beta(self, n, C_1):  # n = the class of clusters, C_1 is monomer conc.

        cdef double R_n
        R_n = ((3*n*self.mon_vol)/(4*np.pi))**(1./3.)

        return 4*np.pi*(R_n**2/(R_n))*(self.D/self.mon_vol)*C_1
        #return 4*np.pi*(R_n**2/(R_n + self.kappa))*(self.D/self.mon_vol)*C_1 
        # Assume diffusive regime - neglect kappa.
     
      
    cdef double sigma(self, n):      # Surface free energy, gen. capillary approx.

        cdef double n_0
        n_0 = ((32e-30*np.pi)/(3*self.mon_vol))**(1./3.)    # Tolman, 1949.

        return self.sigma_const*(1 + (n_0/float(n))**(1./3.))**-2
          
          
    cdef double alpha(self, n):      # Calculate evaporation rate from (n+1) clusters.

        cdef double exponent, R_n
        
        exponent = (36*np.pi*self.mon_vol**2)**(1./3.)        * \
                   (((n+1)**(2./3.)*self.sigma(n+1)           - \
                   n**(2./3.)*self.sigma(n) - self.sigma(1))) / \
                   (self.k_B*self.T)
               
        R_n = ((3*n*self.mon_vol)/(4*np.pi))**(1./3.)

        #return 4*np.pi*(R_n**2/(R_n + self.kappa)) * \
               #(self.D/self.mon_vol)*np.exp(exponent)
        return 4*np.pi*(R_n**2/(R_n)) * (self.D/self.mon_vol)*np.exp(exponent)
               
    
    
    @staticmethod
    def tridiag(a, b, c, k1=-1, k2=0, k3=1):
        return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


    cdef np.ndarray generate_update_matrix(self, C_1):
        '''Function generating the matrix that describes our linear system. 
        Work within the small cluster regime, for now.'''
        
        cdef int i
        cdef np.ndarray diags, lower_diags, upper_diags
        
        diags = np.zeros(self.n_class)
        lower_diags = np.zeros(self.n_class - 1)
        upper_diags = np.zeros(self.n_class - 1)
        
        for i in xrange(self.n_class-1):
            diags[i]       = -(self.beta(i+1, C_1) + self.alpha(i+1))
            lower_diags[i] = self.beta(i+1, C_1)
            upper_diags[i] = self.alpha(i+2)

        diags[-1] = -self.alpha(self.n_class)

        M = ClusterSystem.tridiag(lower_diags, diags, upper_diags)
        
        for i in range(1, self.n_class-1):
            M[0][i] += self.alpha(i+1) - self.beta(i+1, C_1)
        
        M[0,0]      = -2*self.beta(1, C_1)
        M[0, -1]    = self.alpha(self.n_class)
    
        return M
        
        
    @staticmethod
    def calc_mass(x):
        return sum((index+1)*value for index, value in enumerate(x))


    cdef np.ndarray solve_system(self, h, N_ITER):
        '''Solve the system, given C_init. fast=True will return ONLY the final
           distributions, not a distribution for each timestep.'''
        soln = [self.C_init]
        times = np.linspace(0, h*N_ITER, N_ITER)


        for t in times:
            M = self.generate_update_matrix(C_1=soln[-1][0])
            x_new = ODE_int.RK4(M, soln[-1][0:-1], t, 1)[0]
            soln.append(x_new)
        
        return np.asarray(soln)
        
    def wrap_solve(self, h, N_ITER):
        return self.solve_system(h, N_ITER)
