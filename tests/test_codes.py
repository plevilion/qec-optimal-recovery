import numpy as np
import pytest

from qec.codes import basis_states


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
