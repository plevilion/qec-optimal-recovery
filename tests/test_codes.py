import numpy as np


def test_bit_flip_code_corrects_single_error():
    # 3-qubit repetition code: parity checks Z1Z2 and Z2Z3
    H = np.array([[1, 1, 0], [0, 1, 1]])
    for flipped in range(3):
        error = np.zeros(3, dtype=int)
        error[flipped] = 1
        syndrome = H @ error % 2
        lookup = {(1, 0): 0, (1, 1): 1, (0, 1): 2}
        assert lookup[tuple(syndrome)] == flipped
