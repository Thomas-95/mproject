import numpy as np
import scipy

# -----------------------------------------------------------------------------
# Module written for solving linear systems in a matrix format: A.x = b using
# the Jacobi method.
# Solution taken from Wiki entry: https://en.wikipedia.org/wiki/Jacobi_method
# -----------------------------------------------------------------------------


def Jacobi_solve(A, x, b):
    '''Insert form of below here for use in cluster_gas.py
    Arguments: matrix, x vector of unknowns, b vector of values. ''' 

    pass

''' It is then passed into an intergrator which solves the rate values for time,
before being passed back into the Jacobi_solve function. This generates the 
change in concentration gradients over time, which can be plotted. 
scipy.integrate.odeint(C, C_init, t)'''


ITERATION_LIMIT = 100

# initialize the matrix
A = np.array([[10., -1., 2., 0.],
              [-1., 11., -1., 3.],
              [2., -1., 10., -1.],
              [0.0, 3., -1., 8.]])
# initialize the RHS vector
b = np.array([6., 25., -11., 15.])

# prints the system
print("System:")
for i in range(A.shape[0]):
    row = ["{}*x{}".format(A[i, j], j + 1) for j in range(A.shape[1])]
    print(" + ".join(row), "=", b[i])
print()

x = np.zeros_like(b)
for it_count in range(ITERATION_LIMIT):
    print("Current solution:", x)
    x_new = np.zeros_like(x)

    for i in range(A.shape[0]):
        s1 = np.dot(A[i, :i], x[:i])
        s2 = np.dot(A[i, i + 1:], x[i + 1:])
        x_new[i] = (b[i] - s1 - s2) / A[i, i]

    if np.allclose(x, x_new, atol=1e-10):
        break

    x = x_new

print("Solution:")
print(x)
error = np.dot(A, x) - b
print("Error:")
print(error)
