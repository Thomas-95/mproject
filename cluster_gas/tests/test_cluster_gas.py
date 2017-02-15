import env, numpy as np
from lib import cluster_system as cs


def test_mass_conservation():

    test_system = cs.ClusterSystem(n_class=200, T=200.)
    test_system.solve_system(N_ITER=300)
    
    x = cs.ClusterSystem.calc_mass(test_system.soln[0])
    y = cs.ClusterSystem.calc_mass(test_system.soln[-1])

    assert(np.allclose(x, y, rtol=1e-10))
    
    
def test_steady_state():

    test_system = cs.ClusterSystem(n_class=4, T=200.)
    test_system.solve_system(N_ITER=3000)
    soln = test_system.soln
    beta = test_system.beta
    alpha = test_system.alpha
    
    assert(np.allclose(soln[-1][0]*beta(1, C_1=soln[-1][0]),
                       soln[-1][1]*alpha(2),                
                       rtol = 1e-5))
    assert(np.allclose(beta(2, C_1=soln[-1][0])*soln[-1][1]/ alpha(3),
                       soln[-1][2],                  
                       rtol = 1e-5))
    assert(np.allclose(beta(3, C_1=soln[-1][0])*soln[-1][2]/ alpha(4),
                       soln[-1][3],     
                       rtol = 1e-5))
