import numpy as np
from lib import cluster_system
import matplotlib.pyplot as plt

n_class = 5
T = 276
N_ITER = 500

if __name__ == "__main__":

    sys1 = cluster_system.ClusterSystem(n_class, T)
    sys1.solve_system(N_ITER)
    
    print sys1.soln[-1][:-1]

    times = np.linspace(0, 0.8e-12, N_ITER + 1)
    plt.figure(1)
    for i in range(n_class):
        plt.plot(sys1.soln[:, i])
    plt.figure(2)
    plt.title("Total mass of system over time.")
    plt.ylabel("Total mass of clusters.")
    plt.xlabel("Time (s)")
    print 
    plt.plot(times, sys1.calc_mass([sys1.soln[:,i] for i in range(n_class)]))
    plt.show()
    
    # Check if system is in steady state:
    
    n=n_class
    R_n = ((3*n*pr.mon_vol)/(4*np.pi))**(1./3.)
    kappa = pr.D/(pr.jump_freq*pr.d_jump)
    b_N = 4*np.pi*(R_n**2/(R_n + kappa))*(pr.D/pr.mon_vol) # From beta() funct.
    
    exponent = (36*np.pi*pr.mon_vol**2)**(1./3.) * \
               ((n_class)**(2./3.)*pr.sigma(n_class) - \
               (2./3.)*pr.sigma(2.) - n_class*pr.sigma(1.))/(pr.k_B * pr.T)

    print (b_N/np.exp(exponent)) * soln[-1][0]**n_class
    print soln[-1][-2]
    
    #print b_1 * soln[-1][0]**2, pr.alpha(2) * soln[-1][1]
    #print b_1 * soln[-1][0] * soln[-1][1], pr.alpha(3) * soln[-1][2]
