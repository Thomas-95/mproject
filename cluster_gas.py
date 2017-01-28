import numpy as np

# Define some constants.
D = 1.          # The diffusion coefficient. 
omega = 1.      # The monomer volume.
gamma = 1.      # The jump frequency.
T = 273.        # The temperature of the system. 
d_jump = 1.     # The jump distance for monomer to attach to cluster. 
k_B = 1.38e-23  # The Boltzmann constant.
#C_0 = N_0/N_s   # Fraction of accessible sites. 
C_0 = 1.        # Temporary value. 

# -----------------------------------------------------------------------------
# Function to calculate the condenstation & evaporation rates.
# Will return to this later, once some sort of matrix system is set up.
# -----------------------------------------------------------------------------

def beta(n, C_1):  # C_1 is monomer concentration, n is the class of clusters. 

    R_n = ((3*n*omega)/(4*np.pi))**(1./3.)
    kappa = D/(gamma*d_jump)

    return 4*np.pi*(R_n**2/(R_n + kappa))*(D/omega)*C_1

def sigma(n):   # Calculate surface free energy barrier. 

    return 1.

def alpha(n):  # Calculate evaporation rate from (n+1) clusters.

    exponent = (36*np.pi*omega**2)**(1./3.) * \
               (((n+1)**(2./3.)*sigma(n+1)  - \
               n**(2./3.)*sigma(n) - sigma(1))) / (k_B*T)

    #return beta(n, 1.)*C_0*np.exp(exponent)
    return 6.

# -----------------------------------------------------------------------------
# Setup the matrix for use in matrix equation.
# -----------------------------------------------------------------------------

# Work within the small cluster regime, for now.

n_class = 6

M = np.zeros(shape=(n_class, n_class))

# We need to iterate across the rows, filling in coefficients, as required.

C_1 = 0.1
for i in range(1, n_class-1):
    M[i-1, i-2] = beta(i-1, C_1)           # The lower diagonal components.
    M[i, i] = -(beta(i, C_1) + alpha(i))   # The diagonal components.
    M[i, i+1] = alpha(i)                   # The upper diagonal components.

M[n_class-1, n_class-1] = -alpha(n_class-1)     # Final diagonal.
M[n_class-1, n_class-2] = beta(n_class-2, C_1)  # Final lower diagonal.
M[0,0] = -2*beta(1, C_1)

M[0,1] = 2*alpha(2) - beta(2, C_0)         # Sort out first row of the matrix.
for i in range(2, n_class-1):
    M[0, i] = i

print M
