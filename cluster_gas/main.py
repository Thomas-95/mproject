import numpy as np
from lib import ODE_int, generate_matrix, parameters
import matplotlib.pyplot as plt

n_class = 5
C_init = np.zeros(shape=(n_class+1))
C_init[0] = 1.

def solve_system(n_class=n_class):

    soln = [C_init]
    
    times = np.linspace(0, 0.8e-12, 800)
    for t in times:
        M = generate_matrix.update_matrix(n_class=n_class)
        soln.append(ODE_int.RK4(M, soln[-1][0:-1], t, 1)[0])
    
    return np.asarray(soln)

if __name__ == "__main__":
    soln = solve_system()

    times = np.linspace(0, 0.8e-12, 801)
    plt.figure(1)
    plt.plot(soln[:,0])
    plt.plot(soln[:,1])
    plt.plot(soln[:,2])
    plt.plot(soln[:,3])
    plt.plot(soln[:,4])
    plt.figure(2)
    plt.figure(2)
    plt.title("Total mass of clusters over time.")
    plt.ylabel("Total mass of clusters.")
    plt.xlabel("Time (s)")
    plt.plot(times, soln[:,0]+2*soln[:,1]+3*soln[:,2]+4*soln[:,3]+5*soln[:,4])
    plt.show()
    
