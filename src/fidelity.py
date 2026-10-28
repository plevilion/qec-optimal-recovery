import numpy as np
from codes import ket_to_rho

def entanglement_fidelity(kraus):
    # computes the entanglement fidelity of the mixed state rho = I/2 and a channel with Kraus operators kraus 
    d = np.shape(kraus[0])[1]
    return float(sum(abs(np.trace(kn))**2 for kn in kraus)/d**2)

def entanglement_fidelity_2(kraus):
    # same fidelity but with different formula 
    # only works for 1 qubit 
    phi = (np.kron([1,0], [1,0]) + np.kron([0,1], [0,1]))/np.sqrt(2)
    rho = ket_to_rho(phi)
    rho_out = sum(np.kron(np.eye(2), ks)@rho@(np.kron(np.eye(2), ks)).conj().T for ks in kraus)
    return float(np.real(np.vdot(phi, rho_out@phi)))