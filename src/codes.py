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


def leung_isometry() -> np.ndarray:
    """Return the 16x2 encoding isometry V = |0_L><0| + |1_L><1| of the Leung code.

    The codewords are |0_L> = (|0000> + |1111>)/sqrt(2) and
    |1_L> = (|0011> + |1100>)/sqrt(2), so the columns of V are |0_L> and |1_L>.
    """
    zero_L = basis_states("0000") + basis_states("1111")
    one_L = basis_states("0011") + basis_states("1100")
    zero_L /= np.linalg.norm(zero_L)
    one_L /= np.linalg.norm(one_L)
    return np.column_stack([zero_L, one_L])

def code_projector(V: np.ndarray) -> np.ndarray:
    return V @ V.conj().T
