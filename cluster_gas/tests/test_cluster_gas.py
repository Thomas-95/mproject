import env, numpy as np
from lib import cluster_system as cs


def test_mass_conservation():

    test_system = cs.ClusterSystem(n_class=200, T=200.)
    test_system.solve_system(N_ITER=300)
    
    x = cs.ClusterSystem.calc_mass(test_system.soln[0])
    y = cs.ClusterSystem.calc_mass(test_system.soln[-1])

    assert(np.allclose(x, y, rtol=1e-10))
