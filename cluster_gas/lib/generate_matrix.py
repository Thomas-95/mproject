import numpy as np
from parameters import *

# -----------------------------------------------------------------------------
# Module written for generating the matrix that describes our linear system. 
# Work within the small cluster regime, for now.
# -----------------------------------------------------------------------------

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

def update_matrix(n_class):

    '''We set the diagonal values, then the lower & upper diagonals, then the 
       first row.'''
    # These coefficients depend on C_1, so must be updated at each iteration. 
    # FIXME This capacity does not currently exist.

    # Rather than recalculating the same values, it would be better to calculate
    # them all at once then access the values - I am recalculating the same
    # values, at times.

    diags = [0 for i in range(n_class)]
    lower_diags = [beta(i+1) for i in range(n_class-1)]   # n+1 for 0 indexing.
    upper_diags = [alpha(i+1) for i in range(n_class-1)]

    diags[-1] = -alpha(n_class)
    for i in range(1, n_class-1):
        diags[i] = -(beta(i+1) + alpha(i+1)) 

    M = tridiag(lower_diags, diags, upper_diags)

    M[0][1:-1] += alpha(i) - beta(i)
    M[0,0] = -2*beta(i)
    M[0, -1] = alpha(n_class)

    return M
