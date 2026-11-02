import numpy as np
from scipy.optimize import minimize_scalar

from robustness import mismatched_petz_fidelity

def optimal_petz_calibration(gamma_true):

    func = lambda gam0 : -mismatched_petz_fidelity(gamma_true, gam0)
    results = minimize_scalar(func, bounds= (1e-3*gamma_true, gamma_true), method="bounded", options={"xatol": 1e-7})
    return results.x, -results.fun