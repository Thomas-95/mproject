import pytest, env, numpy as np
from lib import cluster_system as cs

def beta_test(*args):  return 0.1
def alpha_test(*args): return 1.0
    
test_system = cs.ClusterSystem(n_class=5, T=1.)
test_system.beta = beta_test
test_system.alpha = alpha_test

def test_generate_update_matrix():

    x = test_system.generate_update_matrix(C_1=1.)
    M = np.array([[-0.2,  1.9,  0.9,  0.9,  1. ],
                  [ 0.1, -1.1,  1. ,  0. ,  0. ],
                  [ 0. ,  0.1, -1.1,  1. ,  0. ],
                  [ 0. ,  0. ,  0.1, -1.1,  1. ],
                  [ 0. ,  0. ,  0. ,  0.1, -1. ]])
                
    assert (x == M).all()
