import numpy as np
import pytest

from recovery import baseline_recovery_kraus
from codes import rand_rho, leung_isometry


def test_dim_base_kraus():

    # testing if the we have the correct number and the correct dimension for the kraus operators of the baseline recovery channel

    R = np.array(baseline_recovery_kraus())
    assert R.shape == (29, 2, 16)



def test_trace_preservation_base_recov():
    
    # testing if the baseline recovery channel is trace preserving 

    R = baseline_recovery_kraus()
    np.testing.assert_allclose(sum(K.conj().T @ K for K in R), np.eye(16), atol = 1e-12)


RHOS = [rand_rho(2, rng = np.random.default_rng(123)), rand_rho(2, rng = np.random.default_rng(231)), rand_rho(2, rng = np.random.default_rng(312))]

@pytest.mark.parametrize('rho', RHOS)
def test_noiseless_decoding(rho):

    # testing if the baseline recovery channel decodes perfeclty when there is no noise

    V = leung_isometry()
    kraus = baseline_recovery_kraus()
    encoded = V@rho@V.conj().T
    decoded = sum(R @ encoded @ R.conj().T for R in kraus)
    np.testing.assert_allclose(decoded, rho, atol=1e-12)