from lib import cluster_system as cs, plot_system

if __name__ == "__main__":

    #from results import cluster_distribution
    
    # Check if system is in steady state:
    import numpy as np
    
    n_class = 40
    T = 200.
    
    #NEED AN ANALYTIC SOLUTION.
    
    system = cs.ClusterSystem(n_class=n_class, T=T)
    system.C_init[0] = 3.35e28  # Number of molecules in metre cubed of water
    
    import matplotlib.pyplot as plt
    from lib import plot_system
    #plot_system.plot_clusters(system)
    #plot_system.plot_mass(system)
    #num_iters = np.linspace(1e2, 1e5, num=5)
    num_iters = [5e4]
    print num_iters
    for n in num_iters:
        system.solve_system(n)
        soln = system.soln
        plt.semilogy(range(1, n_class+1), soln[-1][0:-1], \
                     label=r"$N_{ITER}$=%i" %n)
    plt.legend()
    plt.show()
    
    '''D = system.D
    d_jump = system.d_jump
    jump_freq = system.jump_freq
    mon_vol = system.mon_vol
    k_B = system.k_B
    
    def b_N(n=n_class):
        R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)
        kappa = 0#D/(jump_freq*d_jump)
        return 4*np.pi*(R_n**2/(R_n + kappa))*(D/mon_vol) # From beta() funct.
    
    a = 36*np.pi*mon_vol**2
    b = (n_class)**(2./3.)*sigma(n_class)
    c = (2./3.)*sigma(2.) - (n_class)*sigma(1.)
    d = k_B * T
    
    exponent = (a*b*c)/d
    
    print soln[-1][0]**n_class / np.exp(exponent)
    print soln[-1][-2]'''
