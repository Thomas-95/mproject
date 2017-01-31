import numpy as np
from parameters import *

# -----------------------------------------------------------------------------
# Module written for generating the matrix that describes our linear system. 
# Work within the small cluster regime, for now.
# -----------------------------------------------------------------------------

def update_matrix(n_class):

    '''We iterate across the rows, filling in coefficients, as required.'''
    # These coefficients depend on C_1, so must be updated at each iteration. 
    # FIXME This capacity does not currently exist.

    M = np.zeros(shape=(n_class, n_class))

    for i in range(1, n_class-1):
        M[i-1, i-2] = beta(i-1, C_1)             # Lower diagonal components.
        M[i, i]     = -(beta(i, C_1) + alpha(i)) # Diagonal components.
        M[i, i+1]   = alpha(i)                   # Upper diagonal components.

    M[n_class-1, n_class-1] = -alpha(n_class-1)     # Final diagonal.
    M[n_class-1, n_class-2] = beta(n_class-2, C_1)  # Final lower diagonal.
    M[0,0] = -2*beta(1, C_1)
    
    M[0,1] = 2*alpha(2) - beta(2, C_0)           # Sort out first row of matrix.
    for i in range(2, n_class-1):
        M[0, i] = i

    M[0, n_class-1] = alpha(4)

    return M
