from lib import cluster_system as cs, plot_system

if __name__ == "__main__":

    #from results import log_distributions
    #from results import scaled_system
    #from results import constant_update_matrix
    
    A = 6.022e23
    h = 1e-50 * A
    n_class = 10
    nums = range(1, n_class+1)
    T = 250
    N_ITER = 2e8
    
    
    import time
    start_time = time.time()
    system = cs.ClusterSystem(n_class, T, 5.6e4)
    soln = system.wrap_solve(h, N_ITER)
    print "Time (s) for simulation:", time.time() - start_time
    print "Time (s) at simulation termination:", h * N_ITER
    import sys
    print "Memory (GB) taken by solution:", float(sys.getsizeof(soln))/1073741824.
    
    import matplotlib.pyplot as plt
    plt.semilogy(nums, soln[-1]*A)
    
    import numpy as np
    from scipy.stats import linregress as linr  
    print linr(nums, np.log(soln[-1]*A))
    
    #plt.figure(2)
    #times = np.linspace(0, h*N_ITER, N_ITER)
    #plt.plot(times, soln[:,-1])
    plt.show()
