from lib import cluster_system as cs, plot_system
from datetime import datetime

if __name__ == "__main__":

    h = 1e-50
    N_ITER = 1e6

    b = cs.ClusterSystem(5,5)
    import time
    start_time = time.time()
    soln = b.wrap_solve(h, N_ITER)
    print time.time() - start_time
    
    n_class = 50
    T = 100.


    '''system = cs.ClusterSystem(n_class=n_class, T=T)
    system.C_init[0] = 3.35e28  # Number of molecules in metre cubed of water
    
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import linregress as linr
    
    num_iters = [3e2]
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
        plt.semillogy(range(1, n_class+1), soln[-1][0:-1], \
                     label=r"$N_{ITER}$=%i" %n)
        line = linr(nums, log_soln)
        print line
        r_vals.append(line[2])
    
    
    plt.legend(loc=3)'''
    '''plt.figure(2)
    plt.title(r"Variation in linear fit with $N_{iter}$")
    plt.xlabel(r"$N_{iter}$")
    plt.ylabel(r"Regression $r$ value")
    #print np.log(num_iters)
    #print r_vals
    plt.semilogx(num_iters, r_vals)'''
    #plt.show()
