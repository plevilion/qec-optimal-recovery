import numpy as np
import pytest

from channels import noisy_kraus_operator, encoded_noisy_channel
from codes import basis_states, leung_isometry, code_projector, ket_to_rho


GAMMAS =[0.0, 0.001, 0.1, 0.32, 0.67, 1.0]

@pytest.mark.parametrize("gamma", GAMMAS)
def test_noisy_kraus_preserves_trace(gamma):
    bits = [format(i, "04b") for i in range(16)]
    ks= [noisy_kraus_operator(s, gamma) for s in bits]
    np.testing.assert_allclose(sum(K.conj().T @ K for K in ks), np.eye(16), atol=1e-12)


def rand_rho(dim=2, rng=None):
    # creates a random density matrix 
    if rng is None:
        rng = np.random.default_rng(1234)

    M = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))

    rho = M @ M.conj().T
    return rho / np.trace(rho)


RHOS = [ket_to_rho(np.array([1,0])), ket_to_rho(np.array([0,1])), ket_to_rho(np.array([1/np.sqrt(2), 1/np.sqrt(2)])), np.eye(2)/2, rand_rho()]

@pytest.mark.parametrize("rho", RHOS)
@pytest.mark.parametrize("gamma", GAMMAS)
def test_noisy_encoding_channel_preserves_density(rho, gamma, atol=1e-12):
    transmitted_rho = encoded_noisy_channel(rho, gamma)
    assert np.isclose(np.trace(transmitted_rho), 1.0, atol=atol)
    np.testing.assert_allclose(transmitted_rho.conj().T, transmitted_rho, atol=atol)
    assert np.linalg.eigvalsh(transmitted_rho).min() >= -atol
    


@pytest.mark.parametrize("rho", RHOS)
def test_no_noise_is_just_encoding(rho):
    V = leung_isometry()
    np.testing.assert_allclose(encoded_noisy_channel(rho, 0.0), V @ rho @ V.conj().T, atol=1e-12)


@pytest.mark.parametrize("rho", RHOS)
def test_full_damping_sends_everything_to_0000(rho):
    zero = basis_states("0000")
    np.testing.assert_allclose(encoded_noisy_channel(rho, 1.0), np.outer(zero, zero), atol=1e-12)