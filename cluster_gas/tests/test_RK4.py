import pytest, env
import numpy as np
from lib import ODE_int

def test_RK4():                      # Example of RK4 for a simple pendulum.

    D_hat = 0.                       # Scaled damping constant.
    theta_0 = 0.1                    # Initial displacement. v_0 = 0. 
    y = np.array([0,theta_0])
    L = np.array([[-D_hat,-1],[1,0]])
    soln = ODE_int.RK4(L, y, h=0.02, n_iter=2500)
    assert np.allclose(soln[-1][1],  0.1*np.cos(50), rtol=1e-10)
    assert np.allclose(soln[-1][0], -0.1*np.sin(50), rtol=1e-10)
