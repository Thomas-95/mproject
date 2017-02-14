import env, matplotlib.pyplot as plt
from lib import cluster_system as cs

system = cs.ClusterSystem(n_class=10, T=200)
N_LIMIT = 1000
n_vals = range(1, N_LIMIT)
alpha_vals = [system.alpha(n) for n in n_vals]
beta_vals  = [system.beta(n, 1.) for n in n_vals]
plt.figure(1)
plt.plot(n_vals, beta_vals, 'b', label=r'$\beta$ with $C_1=1$')
plt.plot(n_vals, alpha_vals, 'r', label=r'$\alpha$')
plt.legend(loc=2)
plt.show()
