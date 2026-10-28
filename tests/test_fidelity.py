import numpy as np
import pytest

from fidelity import entanglement_fidelity, entanglement_fidelity_2
from channels import kraus_operators


GAMMAS =[0.0, 0.001, 0.1, 0.32, 0.67, 1.0]

@pytest.mark.parametrize("gamma", GAMMAS)
def test_damping_channel_fidelity(gamma):
    # damping channel fidelity
    kraus = [kraus_operators(gamma, 0), kraus_operators(gamma, 1)]
    expected = ((1+np.sqrt(1-gamma))**2)/4
    assert np.isclose(entanglement_fidelity(kraus), expected, rtol = 1e-12)
    assert np.isclose(entanglement_fidelity_2(kraus), expected, rtol = 1e-12)
    


def test_identity_channel_fidelity():
    # identiy channel fidelity 
    kraus = [np.eye(2)]
    assert np.isclose(entanglement_fidelity(kraus), 1.0, rtol=1e-12)
    assert np.isclose(entanglement_fidelity_2(kraus), 1.0, rtol = 1e-12)

PROBA = [0, 0.001, 0.5, 0.67, 1]
@pytest.mark.parametrize("p", PROBA)
def test_bit_flip_fidelity(p):
    # bit flip channel with probability p
    K0 = np.array([[np.sqrt(1-p), 0], [0, np.sqrt(1-p)]])
    K1 = np.array([[0,np.sqrt(p)], [np.sqrt(p), 0]])
    kraus = [K0, K1]
    assert np.isclose(entanglement_fidelity(kraus), 1-p, rtol=1e-12)
    assert np.isclose(entanglement_fidelity_2(kraus), 1-p, rtol=1e-12)



THETAS = [0.0, 0.3, 1.2, np.pi / 2, np.pi]

@pytest.mark.parametrize("theta", THETAS)
def test_z_rotation_fidelity(theta):
    # rotation exp(-i theta Z / 2) = diag(e^{-i theta/2}, e^{i theta/2}), truly complex
    U = np.diag([np.exp(-1j * theta / 2), np.exp(1j * theta / 2)])
    expected = np.cos(theta / 2) ** 2
    assert np.isclose(entanglement_fidelity([U]), expected, atol=1e-12)
    assert np.isclose(entanglement_fidelity_2([U]), expected, atol=1e-12)