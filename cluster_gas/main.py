import numpy as np
from lib import ODE_int, generate_matrix, parameters as pr
import matplotlib.pyplot as plt

N_ITER = 2500
n_class = 3
C_init = np.zeros(shape=(n_class+1))
C_init[0] = 1.

def calc_mass(x):
    return sum((index+1)*value for index, value in enumerate(x))

def solve_system(n_class=n_class):
    soln = [C_init]
    times = np.linspace(0, 0.8e-12, N_ITER)
    
    for t in times:
        M = generate_matrix.update_matrix(n_class=n_class, C_1=soln[-1][0])
        #print M
        x_new = ODE_int.RK4(M, soln[-1][0:-1], t, 1)[0]
        #print x_new
        #x_new[0] += (1-calc_mass(x_new))  # Add this line to conserve mass?
        soln.append(x_new)
    
    return np.asarray(soln)

if __name__ == "__main__":
    soln = solve_system()
    print soln[-1][:-1]

    times = np.linspace(0, 0.8e-12, N_ITER + 1)
    plt.figure(1)
    for i in range(n_class):
        plt.plot(soln[:, i])
    plt.figure(2)
    plt.title("Total mass of clusters over time.")
    plt.ylabel("Total mass of clusters.")
    plt.xlabel("Time (s)")
    plt.plot(times, calc_mass([soln[:,i] for i in range(n_class)]))
    #plt.show()
    
    # Check if system is in steady state:
    
    n=1
    R_n = ((3*n*pr.mon_vol)/(4*np.pi))**(1./3.)
    kappa = pr.D/(pr.jump_freq*pr.d_jump)
    b_1 = 4*np.pi*(R_n**2/(R_n + kappa))*(pr.D/pr.mon_vol) # From beta() funct.
    
    exponent = ((n_class + 1)**(2./3.)*pr.sigma(n_class + 1) - \
               (2./3.)*pr.sigma(2.) - n_class*pr.sigma(1.))/(pr.k_B * pr.T)
               
    print (n_class + 1)**(2./3.)*pr.sigma(n_class + 1) 
    print (2./3.)*pr.sigma(2.)  
    print n_class*pr.sigma(1.)
    print pr.k_B * pr.T
    print (b_1/np.exp(exponent)) * soln[-1][0]**n_class
