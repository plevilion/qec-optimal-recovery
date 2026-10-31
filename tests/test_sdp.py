import numpy as np 
import cvxpy as cp
import pytest

from sdp import objective_matrix, optimal_recovery, choi_petz_from_kraus
from fidelity import entanglement_fidelity
from channels import noisy_kraus_operator
from petz import completed_petz_kraus
from codes import leung_isometry
from recovery import baseline_recovery_kraus


GAMMAS =[0.0, 0.001, 0.1, 0.32, 0.67, 1.0]

@pytest.mark.parametrize("gamma", GAMMAS)
def test_hermitian_obj_mat(gamma):
    C = objective_matrix(gamma)
    np.testing.assert_allclose(C , C.conj().T, atol=1e-12)
    assert C.shape == (32, 32)


GAMMA = [0.1]
@pytest.mark.parametrize("gamma", GAMMA)
def test_optimized_choi(gamma):
    
    J = optimal_recovery(gamma)[1]
    pt = np.trace(J.reshape(16, 2, 16, 2), axis1=1, axis2=3) 
    np.testing.assert_allclose(J , J.conj().T, atol=1e-5)                   # J is hermitian
    assert np.min(np.linalg.eigvalsh(J)) >= -1e-4                           # J is positive semi definite
    np.testing.assert_allclose(pt, np.eye(16), atol = 1e-4)                 # Partial trace of J is identity


@pytest.mark.parametrize("gamma", GAMMA)
def test_fidelity_from_choi(gamma):
    V = leung_isometry()
    bits = [format(i, '04b') for i in range(16)]
    N = [noisy_kraus_operator(bit, gamma) for bit in bits]
    R = completed_petz_kraus(np.eye(2)/2, gamma)
    kraus = [r @ n @ V for r in R for n in N]
    C = objective_matrix(gamma)
    assert np.isclose(np.real(np.trace(C @ choi_petz_from_kraus(gamma))), entanglement_fidelity(kraus))


@pytest.mark.parametrize("gamma", GAMMA)
def test_sdp_is_optimal(gamma):
    V = leung_isometry()
    bits = [format(i, '04b') for i in range(16)]
    noisy_kraus = [noisy_kraus_operator(bit, gamma) for bit in bits]
    noisy_encoding_kraus = [N@V for N in noisy_kraus]
    R = baseline_recovery_kraus()
    R_til = completed_petz_kraus(np.eye(2)/2 , gamma)
    baseline_recov = [r@E for r in R for E in noisy_encoding_kraus]
    petz_recove = [r_til@E for r_til in R_til for E in noisy_encoding_kraus]

    fidelity_baseline = entanglement_fidelity(baseline_recov)
    fidelity_recov = entanglement_fidelity(petz_recove)
    fidelity_sdp = optimal_recovery(gamma)[0]

    assert fidelity_sdp >= fidelity_recov - 1e-6 and fidelity_sdp >= fidelity_baseline - 1e-6


def test_sdp_decodes_perfectly_without_noise():
    assert np.isclose(optimal_recovery(0.0)[0], 1.0, atol=1e-5)