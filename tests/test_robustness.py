import numpy as np
import pytest

from robustness import mismatched_petz_fidelity
from codes import leung_isometry
from channels import noisy_kraus_operator
from petz import completed_petz_kraus
from fidelity import entanglement_fidelity
from sdp import optimal_recovery


GAMMAS =[0.0, 0.001, 0.1, 0.18, 0.32, 0.67, 1.0]
GAMMAS_0 = [0.0, 0.0015, 0.09, 0.17, 0.34, 0.66, 1.0]


@pytest.mark.parametrize("gamma", GAMMAS)
def test_match_recovery(gamma):
    
    R = completed_petz_kraus(np.eye(2)/2, gamma)
    V = leung_isometry()
    bits = [format(i, "04b") for i in range(16)]
    E = np.array([noisy_kraus_operator(bit, gamma)@V for bit in bits])
    K = [r@e for r in R for e in E]
   
    assert np.isclose(mismatched_petz_fidelity(gamma, gamma), entanglement_fidelity(K), atol=1e-12)


@pytest.mark.parametrize("gamma_true", GAMMAS)
@pytest.mark.parametrize("gamma_assumed", GAMMAS_0)
def test_mismatched_fidelity_is_suboptimal(gamma_true, gamma_assumed):
    f_opt = optimal_recovery(gamma_true)[0]
    assert mismatched_petz_fidelity(gamma_true, gamma_assumed) <= f_opt + 1e-6