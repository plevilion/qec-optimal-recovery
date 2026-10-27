import numpy as np
from functools import reduce
from codes import leung_isometry

def tensor(tns):
    return reduce(np.kron, tns)

def kraus_operators(gamma, index):
    if index == 0:
        return np.array([[1, 0], [0, np.sqrt(1-gamma)]])

    elif index == 1:
        return np.array([[0, np.sqrt(gamma)], [0, 0]])

def amplitude_damping_channel(rho, gamma):
    K0 = kraus_operators(gamma, 0)
    K1 = kraus_operators(gamma, 1)
    return K0 @ rho @ K0.conj().T + K1 @ rho @ K1.conj().T

def noisy_kraus_operator(bits, gamma):
    return tensor(kraus_operators(gamma, int(bit)) for bit in bits)

def noisy_channel(rho, gamma):
    bit = [format(i,"04b") for i in range (16)]
    return sum(noisy_kraus_operator(bits, gamma)@rho@noisy_kraus_operator(bits,gamma).conj().T for bits in bit )

def encoded_noisy_channel(rho, gamma):
    bit = [format(i,"04b") for i in range (16)]
    V = leung_isometry()
    return sum(noisy_kraus_operator(bits, gamma)@V@rho@V.conj().T@noisy_kraus_operator(bits,gamma).conj().T for bits in bit )