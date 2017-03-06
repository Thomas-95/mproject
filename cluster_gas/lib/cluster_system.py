import numpy as np, ODE_int


class ClusterSystem():

    k_B = 1.38e-23                  # Boltzmann constant. 
    
    def __init__(self, n_class, T, E=0., D=2.32e-9, mon_vol=2.992e-29, \
                       d_jump=2.5e-10, jump_freq=1.e13):
                       
        self.n_class = n_class      # Number of cluster sizes to sample.
        self.T = T                  # Temperature of the system. 
        self.mon_vol = mon_vol      # Monomer volume = one water molecule vol.
        self.C_0 = 1.               # All sites poss for homogeneous nucleation.
        self.d_jump = d_jump        # Jump distance for mon to attach w/cluster.
        self.jump_freq = jump_freq  # Taken as vibration freq of water molecule.
        self.sigma_const = (28.0 + (T-273.15)/4.)*1e-4    # -36 < T < 0 Celcius.
        self.C_init = np.zeros(shape=(n_class+1))
        self.C_init[0] = 1.0        # Set initial condition - all monomers.
        self.mon_mass = 2.992e-26
        self.E = E                  # Electric field strength.
        if E == 0.:
            self.D = D
        else:
            r      = 2.5e-20
            self.D = D*(1-np.exp(-(self.k_B*T)/(np.pi*r*(self.E**0.5)))) 
      
      
    def sigma(self, n):      # Surface free energy, gen. capillary approx.

        n_0 = ((32e-30*np.pi)/(3*self.mon_vol))**(1./3.)    # Tolman, 1949.

        return self.sigma_const*(1 + (n_0/float(n))**(1./3.))**-2

               
    def alpha(self, n):     
        k_B = 1.35e-23 
        mfp = 2.5e-10
        exponent = (36*np.pi*self.mon_vol**2)**(1./3.)        * \
                   (((n+1)**(2./3.)*self.sigma(n+1)           - \
                   n**(2./3.)*self.sigma(n) - self.sigma(1))) / \
                   (k_B*self.T)
               
        R_n = ((3*n*self.mon_vol)/(4*np.pi))**(1./3.)

        return 3*exponent*(np.pi)*self.D*R_n**2/(mfp*self.mon_vol)
        
        
    def beta(self, n, C_1):
        k_B = 1.35e-23 
        mfp = 2.5e-10
        R_n = ((3*n*self.mon_vol)/(4*np.pi))**(1./3.)
        return 3*C_1*(np.pi)*self.D*R_n**2/(mfp*self.mon_vol)
    
    
    @staticmethod
    def tridiag(a, b, c, k1=-1, k2=0, k3=1):
        return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)


    def generate_update_matrix(self, C_1):
        '''Function generating the matrix that describes our linear system. 
        Work within the small cluster regime, for now.'''
    
        diags = [0 for i in range(self.n_class)]
        lower_diags = [self.beta(i+1, C_1) for i in range(self.n_class-1)]
        upper_diags = [self.alpha(i+2) for i in range(self.n_class-1)]

        diags[-1] = -self.alpha(self.n_class)
        for i in range(1, self.n_class-1):
            diags[i] = -(self.beta(i+1, C_1) + self.alpha(i+1))

        M = self.tridiag(lower_diags, diags, upper_diags)
        
        for i in range(1, self.n_class-1):
            M[0][i] += self.alpha(i+1) - self.beta(i+1, C_1)
        
        M[0,0]      = -2*self.beta(1, C_1)
        M[0, -1]    = self.alpha(self.n_class)
    
        return M
        
        
    @staticmethod
    def calc_mass(x):
        return sum((index+1)*value for index, value in enumerate(x))


    def solve_system(self, h, N_ITER):
        soln = [self.C_init]
        #times = np.linspace(0, self.kappa, N_ITER)
        times = np.linspace(0, h*N_ITER, N_ITER)
    
        for t in times:
            M = self.generate_update_matrix(C_1=soln[-1][0])
            x_new = ODE_int.RK4(M, soln[-1][0:-1], t, 1)[0]
            soln.append(x_new)

        self.soln = np.asarray(soln)
