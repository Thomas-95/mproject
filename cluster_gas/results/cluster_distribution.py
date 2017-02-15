import env, matplotlib.pyplot as plt, numpy as np
from lib import cluster_system as cs, plot_system

n_class = 15
T = 200
N_ITER = 500

system = cs.ClusterSystem(n_class, T)
system.solve_system(N_ITER)

plot_system.plot_mass(system)
plot_system.plot_distribution(system)
plot_system.plot_clusters(system)
