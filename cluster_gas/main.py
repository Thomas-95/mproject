import numpy as np
from lib import ODE_int, generate_matrix, parameters
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #Example of RK4 for a simple pendulum.
    '''D_hat = 0.2
    theta_0 = 0.1
    y = np.array([0,theta_0])
    L = np.array([[-D_hat,-1],[1,0]])

    soln = ODE_int.RK4(L, y, h=0.02, n_iter=2500)
    plt.plot(soln[:,2],soln[:,1])
    plt.show()'''

    '''temps = np.linspace(263.15, 273.10, 1000)
    radii = [parameters.critical_radius(T) for T in temps]
    print radii
    def reciprocal(t):
        return 100./(-t+273.15)
    plt.plot(temps, radii, 'g', label="n*(T)")
    plt.plot(temps, reciprocal(temps), 'r--', label="f(T)=1/T")
    plt.axvline(x=273.15, label="T=273.15K")
    plt.xlabel("Temperature of the system (K)")
    plt.ylabel("Critical Cluster size")
    plt.legend(loc=2)
    #plt.show()'''
    
    n_class = 5
    M = generate_matrix.update_matrix(n_class=n_class)
    C_init = np.zeros(shape=(n_class))
    C_init[0] = 1.
    timestep = 0.2e-13
    n_iter = 8000
    soln = ODE_int.RK4(M, C_init, timestep, n_iter)
    print soln
    plt.figure(1)
    plt.plot(soln[:,0])
    plt.plot(soln[:,1])
    plt.plot(soln[:,2])
    plt.plot(soln[:,3])
    plt.plot(soln[:,4])
    plt.figure(2)
    plt.title("Total mass of clusters over time.")
    plt.plot(soln[:,0]+2*soln[:,1]+3*soln[:,2]+4*soln[:,3]+5*soln[:,4])
    plt.show()
