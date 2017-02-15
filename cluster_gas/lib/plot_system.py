import matplotlib.pyplot as plt


def plot_mass(system):

    plt.title(r"Total mass of system over iterations")
    plt.ylabel(r"Total mass of system.")
    plt.xlabel(r"Number of iterations")
    plt.plot(system.calc_mass([system.soln[:,i] \
                               for i in range(system.n_class)]))
    
    return plt.show()
    
    
def plot_distribution(system):

    plt.title(r"Final cluster size distribution.")
    plt.plot(range(1, system.n_class+1), system.soln[-1][:-1])
    
    return plt.show()
    
    
def plot_clusters(system):

    plt.title(r"Change in cluster class sizes over iterations.")
    for i in range(system.n_class):
        plt.plot(system.soln[:, i])
        
    return plt.show()
