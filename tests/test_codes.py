import numpy as np
import pytest

from qec.codes import basis_states, leung_isometry


def test_bit_flip_code_corrects_single_error():
    # 3-qubit repetition code: parity checks Z1Z2 and Z2Z3
    H = np.array([[1, 1, 0], [0, 1, 1]])
    for flipped in range(3):
        error = np.zeros(3, dtype=int)
        error[flipped] = 1
        syndrome = H @ error % 2
        lookup = {(1, 0): 0, (1, 1): 1, (0, 1): 2}
        assert lookup[tuple(syndrome)] == flipped


def test_basis_state_matches_kron_product():
    ket = {"0": np.array([1, 0]), "1": np.array([0, 1])}
    for i in range(16):
        bits = format(i, "04b")
        expected = ket[bits[0]]
        for b in bits[1:]:
            expected = np.kron(expected, ket[b])
        np.testing.assert_array_equal(basis_states(bits), expected)


def test_basis_states_are_orthonormal():
    states = np.array([basis_states(format(i, "04b")) for i in range(16)])
    np.testing.assert_array_equal(states.conj() @ states.T, np.eye(16))


@pytest.mark.parametrize("bad", ["", "010", "01010", "01a1"])
def test_rejects_invalid_bit_strings(bad):
    with pytest.raises(ValueError):
        basis_states(bad)


def test_leung_isometry_codewords():
    V = leung_isometry()
    assert V.shape == (16, 2)
    s = 1 / np.sqrt(2)
    zero_L, one_L = np.zeros(16), np.zeros(16)
    zero_L[[0b0000, 0b1111]] = s
    one_L[[0b0011, 0b1100]] = s
    np.testing.assert_allclose(V[:, 0], zero_L)
    np.testing.assert_allclose(V[:, 1], one_L)


def test_leung_isometry_is_isometry():
    V = leung_isometry()
    np.testing.assert_allclose(V.conj().T @ V, np.eye(2), atol=1e-12)
    P = V @ V.conj().T
    np.testing.assert_allclose(P @ P, P, atol=1e-12)
    assert np.isclose(np.trace(P), 2)
