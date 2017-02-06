import parameters

# Will use code to find critical cluster size at different temps here.

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
