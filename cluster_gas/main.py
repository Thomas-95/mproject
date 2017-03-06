from lib import cluster_system as cs, plot_system
from datetime import datetime
startTime = datetime.now()

if __name__ == "__main__":

    n_class = 4
    T = 200.
  
    system = cs.ClusterSystem(n_class=n_class, T=T)
    system.C_init[0] = 3.35e28  # Number of molecules in metre cubed of water
    
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import linregress as linr
    
    
    num_iters = [5, 1e2, 5e2, 1e3, 5e3, 1e4, 3e5]
    #h         = 1e-46  # Close to unstable.
    h = 1e-60
    
    nums = range(1, n_class + 1)
    r_vals = []
    
    for n in num_iters:
        system.solve_system(h=h, N_ITER=n, fast=False)
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
    #print np.log(num_iters)
    #print r_vals
    plt.semilogx(num_iters, r_vals)
    #plt.show()
    
    print datetime.now() - startTime
