import env, matplotlib.pyplot as plt, numpy as np
from lib import cluster_system as cs

n_class = 15
T = 200
N_ITER = 500

system = cs.ClusterSystem(n_class, T)
system.solve_system(N_ITER)

times = np.linspace(0, system.kappa, N_ITER + 1)
plt.figure(1)
plt.title("Chnage in cluster class sizes over time.")
for i in range(n_class):
    plt.plot(system.soln[:, i])
    
plt.figure(2)
plt.title("Total mass of system over time.")
plt.ylabel("Total mass of system.")
plt.xlabel("Time (s)")
plt.plot(times, system.calc_mass([system.soln[:,i] for i in range(n_class)]))

plt.figure(3)
plt.title("Final cluster size distribution.")
plt.plot(range(1, n_class+1), system.soln[-1][:-1])

plt.show()
