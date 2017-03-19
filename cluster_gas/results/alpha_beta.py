import env, matplotlib.pyplot as plt, numpy as np

mon_vol = 2.992e-29
T       = 270.
sigma_const = (28.0 + (T-273.15)/4.)*1e-4
k_B     = 1.38e-23
mfp     = 2.5e-10
D       = 2.32e-9


def calc_D(E):
    if E != 0.:
        r = 2.5e-20
        return D*(1-np.exp(-(k_B*T)/(np.pi*r*(E**0.5)))) 
    else: return D


def sigma(n):

    n_0 = ((32e-30*np.pi)/(3*mon_vol))**(1./3.)

    return sigma_const*(1 + (n_0/float(n))**(1./3.))**-2

               
def alpha(n):
    
    exponent = (36*np.pi*mon_vol**2)**(1./3.)        * \
               (((n+1)**(2./3.)*sigma(n+1)           - \
               n**(2./3.)*sigma(n) - sigma(1))) / \
               (k_B*T)
               
    R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)

    return 3*exponent*(np.pi)*D*R_n**2/(mfp*mon_vol)
    

def beta(n, C_1):

    R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)

    return 4*np.pi*R_n*(D/mon_vol)*C_1
        
        
def izzy_beta(n, C_1):
        
    R_n = ((3*n*mon_vol)/(4*np.pi))**(1./3.)
        
    return 3*C_1*(np.pi)*D*R_n**2/(mfp*mon_vol)



n_vals = np.linspace(1., 100., num=1000)
E_vals = np.linspace(0., 1.e4, num=1e5)
D_vals = [calc_D(E) for E in E_vals]

plt.figure(1)
plt.title(r"Diffusion coefficient $D$ for varying field strengths $|E|$")
plt.xlabel(r"Field strengths $|E|$ $(Vm^{-1})$")
plt.ylabel(r"Diffusion coefficient $D$ $(m^2 s^{-1})$")
plt.semilogy(E_vals, D_vals)

E_vals = [0., 1.e1, 1.e2, 1.e3, 1.e4, 1.e5]#, 1.e6]

plt.figure(2)
plt.xlim([1, 100])
plt.title(r"Condensation coefficient $\beta_n$ for varying field strengths $|E|$")
plt.xlabel(r"Cluster size $n$")
plt.ylabel(r"$\beta_n (m^{-3}s^{-1})$")
for E in E_vals:
    D = calc_D(E)
    plt.semilogy(n_vals, beta(n_vals, C_1=1), label="$|E| = %i Vm^{-1}$"%E)
plt.legend(loc=4)

plt.figure(3)
plt.xlim([1, 100])
plt.title(r"Evaporation coefficient $\alpha_n$ for varying field strengths $|E|$")
plt.xlabel(r"Cluster size $n$")
plt.ylabel(r"$\alpha_n (m^{-3}s^{-1})$")
for E in E_vals:
    D = calc_D(E)
    alpha_vals = [-alpha(n) for n in n_vals]
    plt.semilogy(n_vals, alpha_vals, label="$|E| = %i Vm^{-1}$"%E)
plt.legend(loc=4)
plt.show()

