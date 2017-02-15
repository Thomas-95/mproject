import pytest, env, numpy as np
from lib import cluster_system as cs

def beta_test(n, *args):  return 1*n
def alpha_test(n, *args): return 10*n
    
test_system = cs.ClusterSystem(n_class=5, T=1.)
test_system.beta  = beta_test
test_system.alpha = alpha_test

def test_generate_update_matrix():

    x = test_system.generate_update_matrix(C_1=1)                 
    M = np.array([[-2,  38,  27,  36,  50],
                  [ 1, -22,  30,   0,   0],
                  [ 0,   2, -33,  40,   0],
                  [ 0,   0,   3, -44,  50],
                  [ 0,   0,   0,   4, -50 ]])
             
    assert (x == M).all()
