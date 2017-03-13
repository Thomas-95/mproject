import numpy as np, matplotlib.pyplot as plt
from lib import cluster_system as cs
from scipy.stats import linregress as linr    
    
h = 1e-60
N_ITER = 1e7
n_class = 10
T = 5.
nums = range(1, n_class+1)
    
sys1 = cs.ClusterSystem(n_class, T)
import time
start_time = time.time()
soln1 = sys1.wrap_solve(h, N_ITER)
print time.time() - start_time
    
sys2 = cs.ClusterSystem(n_class, T, 5.6e4)
soln2 = sys2.wrap_solve(h=h*6.022e23, N_ITER=N_ITER)

print soln1[-1]
print soln2[-1]*6.022e23
    

#plt.semilogy(nums, soln1[-1], label="system 1")
#plt.semilogy(nums, soln2[-1]*6.022e23, label="system 2")
#plt.legend()
#plt.show()

  
log_soln1 = np.log(soln1[-1])
log_soln2 = np.log(soln2[-1])
print linr(nums, log_soln1)
print linr(nums, log_soln2)       


# Results for slope within 1% of each other. Could just use a scaled system.
# The Intercept, however, is some way off. 
# This holds for a constant matrix system.
# Need to see the effect of non-constant matrix.
