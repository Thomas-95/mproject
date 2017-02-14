import pytest, env, numpy as np
from lib import cluster_system as cs

def beta_test(n, *args):  return 0.1*n
def alpha_test(n, *args): return 1.0*n
    
test_system = cs.ClusterSystem(n_class=5, T=1.)
test_system.beta  = beta_test
test_system.alpha = alpha_test

def test_generate_update_matrix():

    x = test_system.generate_update_matrix(C_1=1.)
    M = np.array([[-0.2,  3.8,  2.7,  3.6,  5.0 ],
                  [ 0.1, -2.2,  3. ,  0. ,  0.  ],
                  [ 0. ,  0.2, -3.3,  4.0,  0.  ],
                  [ 0. ,  0. ,  0.3, -4.4,  5.  ],
                  [ 0. ,  0. ,  0. ,  0.4, -5.  ]])
                  
    assert (x == M).all()
