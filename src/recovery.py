import numpy as np
from scipy.linalg import orth

from codes import code_projector, leung_isometry

def baseline_recovery_kraus():

    # returns the Kraus operators for the baseline recovery channel

    V = leung_isometry()
    Q = np.eye(16) - code_projector(V)
    q = orth(Q)
    base = [np.array([1,0]), np.array([0,1])]
    kraus = [np.outer(b, q[:, k].conj())/np.sqrt(2) for b in base for k in range(q.shape[1])]
    return [V.conj().T] + kraus