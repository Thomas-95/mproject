import numpy as np, matplotlib.pyplot as plt, time
from lib import cluster_system as cs
from progressbar import ProgressBar
from scipy.stats import linregress as linr 


A = 6.022e23
h = 8e-42 * A
n_class = 50
nums = range(1, n_class+1)
T = 250
N_ITER = 5e7
N_RUNS = 11
start_time = time.time()
soln = np.zeros(n_class)
soln[0] = 5.6e4


plt.figure(1)
pbar = ProgressBar()
for n in pbar(xrange(N_RUNS)):
    system = cs.ClusterSystem(n_class, T, soln)
    soln = system.wrap_solve(h, N_ITER)[-1]
        
    if n%2 == 0:
        t = (n+1)*h*N_ITER
        a, b = '{:.5e}'.format(t).split('e')
        b = int(b)
        plt.semilogy(nums, soln*A, label=r"$t =  %s \times 10^{%s}s$"%(a,b))
        
print "Time (s) to run simulation:", time.time() - start_time
print "Time (s) at simulation termination:", h * N_ITER * N_RUNS
#import sys
#print "Memory (GB) taken by solution:", float(sys.getsizeof(soln))/1073741824.
print linr(nums, np.log(soln*A))
    
plt.title("Convergence of system over time")
plt.xlabel(r"Cluster class size $n$")
plt.ylabel(r"Cluster class concentration $C_n$")
plt.legend(loc=3)
plt.xlim([1, n_class])
plt.show()
