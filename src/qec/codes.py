"""Quantum codes on 4 qubits."""

import numpy as np

N_QUBITS = 4
DIM = 2**N_QUBITS


def basis_states(bits: str) -> np.ndarray:
    """Return the computational basis state |bits> in C^16.

    The first character is the most significant qubit, so "0011" maps to
    index 3, matching the ordering of np.kron(|0>, |0>, |1>, |1>).
    """
    if len(bits) != N_QUBITS or set(bits) - {"0", "1"}:
        raise ValueError(f"expected a {N_QUBITS}-character bit string, got {bits!r}")
    state = np.zeros(DIM, dtype=complex)
    state[int(bits, 2)] = 1
    return state
