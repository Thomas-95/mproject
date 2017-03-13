import numpy as np, matplotlib.pyplot as plt, time
from scipy.stats import linregress as linr
from lib import cluster_system as cs
    
    
h = 1e-50
T = 250.
n_class = 10
N_ITER = 1e5
nums = range(1, n_class + 1)
    
    
system = cs.ClusterSystem(n_class, T)
start_time = time.time()
const_M_soln = system.wrap_solve(h, N_ITER, const_M=1)
time1 = time.time() - start_time
start_time = time.time()
vary_M_soln  = system.wrap_solve(h, N_ITER, const_M=0)
time2 = time.time() - start_time
    
    
print "Time (s) to run with constant M:", time1
print "Time (s) to run with varying M:", time2
print "Constant M is quicker by factor of:", time2/time1
print const_M_soln[-1]
print vary_M_soln[-1]
    
    
print linr(nums, np.log(const_M_soln[-1]))
print linr(nums, np.log(vary_M_soln[-1]))
    
plt.semilogy(nums, const_M_soln[-1], label=r"$M$ constant")
plt.semilogy(nums, vary_M_soln[-1], label=r"$M$ varying")
plt.legend()
plt.show()
