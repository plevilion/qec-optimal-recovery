import numpy as np

from fidelity import entanglement_fidelity
from petz import completed_petz_kraus
from channels import noisy_kraus_operator
from codes import leung_isometry

def mismatched_petz_fidelity(gamma_true, gamma_assumed):
    
    # Entanglement fidelity of Petz recovery constructed for gamma_assumed and applied to amplitude damping gamma_true.
    # The reference state is sigma = I_2 / 2.

    sigma = np.eye(2)/2
    R = completed_petz_kraus(sigma, gamma_assumed)
    bits = [format(i, "04b") for i in range(16)]
    V = leung_isometry()
    E = [noisy_kraus_operator(bit, gamma_true)@V for bit in bits]
    K = [r@e for r in R for e in E]

    return entanglement_fidelity(K)