from lib import cluster_system as cs, plot_system

if __name__ == "__main__":

    #from results import cluster_distribution
    
    n_class = 4
    T = 200.
  
  
    system = cs.ClusterSystem(n_class=n_class, T=T)
    system.C_init[0] = 3.35e28  # Number of molecules in metre cubed of water
    
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import linregress as linr
    
    
    #plot_system.plot_clusters(system)
    #plot_system.plot_mass(system)
    #num_iters = np.linspace(1e2, 1e5, num=5)
    
    
    num_iters = [5]#, 1e2, 5e2]#, 1e3, 5e3, 1e4, 5e4]
    #h         = 1e-46  # Close to unstable.
    h = 1e-60
    
    nums = range(1, n_class + 1)
    r_vals = []
    
    for n in num_iters:
        system.solve_system(h=h, N_ITER=n)
        #print "Time simulated (s):", h*n
        soln = system.soln
        log_soln = np.log(soln[-1][0:-1])
        #print np.gradient(log_soln)
        plt.figure(1)
        plt.title(r"Variation in Size Distribution with $N_{ITER}$")
        plt.xlabel("Cluster size class $n$")
        plt.ylabel("Size class concentration $C_n$")
        plt.semilogy(range(1, n_class+1), soln[-1][0:-1], \
                     label=r"$N_{ITER}$=%i" %n)
        r_vals.append(linr(nums, log_soln)[2])
    
    
    plt.legend(loc=3)
    plt.figure(2)
    plt.title(r"Variation in linear fit with $N_{iter}$")
    plt.xlabel(r"$N_{iter}$")
    plt.ylabel(r"Regression $r$ value")
    print np.log(num_iters)
    print r_vals
    plt.semilogx(num_iters, r_vals)
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
