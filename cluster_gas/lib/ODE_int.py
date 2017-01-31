import numpy as np

# -----------------------------------------------------------------------------
# Module written for solving linear systems in a matrix format: dx/dt = A.x 
# using the methods supplied by linear_solver.
# -----------------------------------------------------------------------------

def RK4(L, x, h, n_iter):
    ''' Arguments: Update matrix for system to be solved, L.
        Numpy array of initial conditions, x.
        Step size increment in time, h.
        No. of iteraions to perform, n_iter.'''

    soln = []

    for n in range(n_iter):
        k_1 = h*np.dot(L,x)                                                   
        k_2 = h*np.dot(L, x + 0.5*k_1)
        k_3 = h*np.dot(L, x + 0.5*k_2)
        k_4 = h*np.dot(L, x + k_3)

        x += (1./6.)*k_1 + (1./3.)*k_2 + (1./3.)*k_3 + (1./6.)*k_4

        soln.append(np.append(x, n*h))

    return np.asarray(soln)
