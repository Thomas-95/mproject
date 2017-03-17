import env, matplotlib.pyplot as plt
from lib import cluster_system as cs

system = cs.ClusterSystem(n_class=10, T=200)
N_LIMIT = 1000
n_vals = range(1, N_LIMIT)
alpha_vals  = [system.alpha(n) for n in n_vals]
beta_vals1  = [system.beta(n, 1.)  for n in n_vals]
beta_vals2  = [system.beta(n, 0.9) for n in n_vals]
beta_vals3  = [system.beta(n, 0.8) for n in n_vals]
beta_vals4  = [system.beta(n, 0.7) for n in n_vals]
plt.figure(1)
plt.plot(n_vals, beta_vals1, label=r'$\beta$ with $C_1=1$')
plt.plot(n_vals, beta_vals2, label=r'$\beta$ with $C_1=0.9$')
plt.plot(n_vals, beta_vals3, label=r'$\beta$ with $C_1=0.8$')
plt.plot(n_vals, beta_vals4, label=r'$\beta$ with $C_1=0.7$')
plt.plot(n_vals, alpha_vals, '--', label=r'$\alpha$')
plt.legend(loc=2)
plt.show()
