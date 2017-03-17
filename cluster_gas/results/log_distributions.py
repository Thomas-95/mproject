import numpy as np, matplotlib.pyplot as plt
from lib import cluster_system as cs
from scipy.stats import linregress as linr

n_class = 10
T = 250
h = 1e-60
system = cs.ClusterSystem(n_class, T)


num_iters = [3e2, 6e3, 9e4]
nums = range(1, n_class + 1)
r_vals = []
    
for n in num_iters:
    soln = system.wrap_solve(h, n)
    #print "Time simulated (s):", h*n
    log_soln = np.log(soln[-1])
    #print np.gradient(log_soln)
    plt.figure(1)
    plt.title(r"Variation in Size Distribution with $N_{ITER}$")
    plt.xlabel("Cluster size class $n$")
    plt.ylabel("Size class concentration $C_n$")
    plt.semilogy(range(1, n_class+1), soln[-1], \
                 label=r"$N_{ITER}$=%i" %n)
    line = linr(nums, log_soln)
    print line
    r_vals.append(line[2])
    
    
plt.legend(loc=3)
plt.figure(2)
plt.title(r"Variation in linear fit with $N_{iter}$")
plt.xlabel(r"$N_{iter}$")
plt.ylabel(r"Regression $r$ value")
#print np.log(num_iters)
#print r_vals
plt.semilogx(num_iters, r_vals)
#plt.show()
