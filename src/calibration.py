import numpy as np
from scipy.optimize import minimize_scalar

from robustness import mismatched_petz_fidelity

def optimal_petz_calibration(gamma_true):

    func = lambda gam0 : -mismatched_petz_fidelity(gamma_true, gam0)
    results = minimize_scalar(func, bounds= (1e-3*gamma_true, 1-1e-4), method="bounded", options={"xatol": 1e-7})
    assert results.success
    return results.x, -results.fun