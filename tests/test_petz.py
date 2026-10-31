import numpy as np
import pytest

from petz import petz_kraus, support_projector, completed_petz_kraus
from channels import encoded_noisy_channel, noisy_kraus_operator
from codes import rand_rho, leung_isometry
from fidelity import entanglement_fidelity



GAMMAS =[0.0, 0.001, 0.1, 0.32, 0.67, 1.0]
RHOS = [rand_rho(2, rng = np.random.default_rng(123)), rand_rho(2, rng = np.random.default_rng(231)), rand_rho(2, rng = np.random.default_rng(312))]


@pytest.mark.parametrize("gamma", GAMMAS)
def test_dim_kraus_petz(gamma):
   
    # making sure the Petz's map kraus operators have correct dimension

    R = petz_kraus(np.eye(2)/2, gamma)
    for k in range(len(R)):
        np.testing.assert_allclose(np.shape(R[k]), [2, 16], atol=1e-12)

@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("rho", RHOS)
def test_recovery_preserves_trace(rho, gamma):

    # testing if the recovery channel is trace preserving on the support of phi(rho)

    R = petz_kraus(rho, gamma)
    np.testing.assert_allclose(sum((R[k].conj().T)@R[k] for k in range(len(R))), support_projector(encoded_noisy_channel(rho, gamma)), atol=1e-8)



@pytest.mark.parametrize("gamma", GAMMAS)
def test_recovery(gamma):

    # testing if the recovery of the transmitted encoded mixed state is still the same

    sigma = np.eye(2)/2
    tau = encoded_noisy_channel(sigma, gamma)
    R = petz_kraus(sigma,gamma)
    recov = sum(R[k]@tau@R[k].conj().T for k in range(len(R)))
    np.testing.assert_allclose(recov, sigma, atol = 1e-12)




@pytest.mark.parametrize("rho", RHOS)
def test_noiseless_recovery(rho):

    # testing if the completed Petz maps recovers initial when there is no noise, i.e. gamma = 0

    R = completed_petz_kraus(np.eye(2)/2, 0.0)
    V = leung_isometry()
    recov = sum(r @ V @ rho @ V.conj().T @ r.conj().T for r in R)
    np.testing.assert_allclose(recov, rho, atol=1e-12)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_recovery_fidelity(gamma):

    # testing if the values for the fidelity of recovery are coherent 

    R = petz_kraus(np.eye(2)/2, gamma)
    V = leung_isometry()
    bits = [format(i, "04b") for i in range(16)]
    E = np.array([noisy_kraus_operator(bit, gamma)@V for bit in bits])
    kraus = [r@e for r in R for e in E]
    F = entanglement_fidelity(kraus)
    
    assert 0.25 - 1e-12 <= F <= 1 + 1e-12            # always in [1/4, 1]
    if gamma == 0.0:
        assert np.isclose(F, 1.0, atol=1e-12)
    if gamma == 1.0:
        assert np.isclose(F, 0.25, atol=1e-12)




@pytest.mark.parametrize("gamma", GAMMAS)
def test_completed_petz_is_trace_preserving(gamma):
    R = completed_petz_kraus(np.eye(2)/2, gamma)
    np.testing.assert_allclose(sum(K.conj().T @ K for K in R), np.eye(16), atol=1e-8)