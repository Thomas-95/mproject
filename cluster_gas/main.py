from lib import cluster_system as cs, plot_system
import numpy as np

if __name__ == "__main__":

    #from results import log_distributions
    #from results import scaled_system
    #from results import constant_update_matrix
    
    A = 6.022e23
    h = 8e-40 * A
    n_class = 50
    nums = range(1, n_class+1)
    T = 250
    N_ITER = 3e8
    N_RUNS = 10
    
        
    import matplotlib.pyplot as plt
    plt.figure(1)
    
    
    import time
    start_time = time.time()
    soln = np.zeros(n_class)
    soln[0] = 5.6e4
    
    from progressbar import ProgressBar
    
    pbar = ProgressBar()
    for n in pbar(xrange(N_RUNS)):
        system = cs.ClusterSystem(n_class, T, soln)
        soln = system.wrap_solve(h, N_ITER)[-1]
        if n%2 == 0:
            plt.semilogy(nums, soln*A, label="Run %i"%(n+1))
        
    print "Time (s) to run simulation:", time.time() - start_time
    print "Time (s) at simulation termination:", h * N_ITER * N_RUNS
    #import sys
    #print "Memory (GB) taken by solution:", float(sys.getsizeof(soln))/1073741824.
    
    #print soln[-1]
    
    
    
    import numpy as np
    from scipy.stats import linregress as linr  
    print linr(nums, np.log(soln*A))
    
    plt.legend()
    plt.show()
