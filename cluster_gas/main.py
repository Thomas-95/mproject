import numpy as np
from lib import ODE_int, generate_matrix
import matplotlib.pyplot as plt

if __name__ == "__main__":
    #Example of RK4 for a simple pendulum.
    D_hat = 0.2
    theta_0 = 0.1
    y = np.array([0,theta_0])
    L = np.array([[-D_hat,-1],[1,0]])

    soln = ODE_int.RK4(L, y, h=0.02, n_iter=2500)
    plt.plot(soln[:,2],soln[:,1])
    #plt.show()

    M = generate_matrix.update_matrix(n_class=5)
    print M
