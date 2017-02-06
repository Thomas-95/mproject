import numpy as np
from lib import ODE_int, generate_matrix, parameters
import matplotlib.pyplot as plt

n_class = 10
C_init = np.zeros(shape=(n_class+1))
C_init[0] = 1.

def calc_mass(x):
    mass = 0.
    for index, value in enumerate(x):
        mass += (index+1)*value
    return mass

def solve_system(n_class=n_class):
    soln = [C_init]
    times = np.linspace(0, 0.8e-12, 8000)
    
    for t in times:
        M = generate_matrix.update_matrix(n_class=n_class, C_1=soln[-1][0])
        x_new = ODE_int.RK4(M, soln[-1][0:-1], t, 1)[0]
        #x_new[0] += (1-calc_mass(x_new))  # Add this line to conserve mass?
        soln.append(x_new)
    
    return np.asarray(soln)

if __name__ == "__main__":
    soln = solve_system()

    times = np.linspace(0, 0.8e-12, 8001)
    plt.figure(1)
    for i in range(n_class-1):
        plt.plot(soln[:, i])
    plt.figure(2)
    plt.figure(2)
    plt.title("Total mass of clusters over time.")
    plt.ylabel("Total mass of clusters.")
    plt.xlabel("Time (s)")
    plt.plot(times, soln[:,0]+2*soln[:,1]+3*soln[:,2]+4*soln[:,3]+5*soln[:,4])
    #plt.plot(times, calc_mass(soln[:,]))
    plt.show()
