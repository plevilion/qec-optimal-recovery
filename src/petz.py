import numpy as np

from channels import encoded_noisy_channel, noisy_kraus_operator
from codes import leung_isometry


def pseudo_inv_sqrt(A, tol = 1e-12):
    
    #computes the pseudo inverse square root of a positive semidefinite symmetric/hermitian matrix

    egvalues, U = np.linalg.eigh(A)
    inv_sqrt = np.zeros_like(egvalues)
    for i in range(len(egvalues)):
        if egvalues[i] > tol:
            inv_sqrt[i] = 1/np.sqrt(egvalues[i])
        else : 
            inv_sqrt[i] = 0
    return U@np.diag(inv_sqrt)@U.conj().T

def psd_sqrt(A):
    w, U = np.linalg.eigh(A)
    return U @ np.diag(np.sqrt(np.clip(w, 0, None))) @ U.conj().T
    

def support_projector(A, tol = 1e-12):
    eigenvalues, U = np.linalg.eigh(A)
    mask = eigenvalues > tol
    U_support = U[:, mask]
    return U_support @ U_support.conj().T


def petz_kraus(sigma, gamma):

    # returns the kraus operators of the Petz recovery map for sigma and the noisy encoding channel 
    # note that the Petz map is not trace preserving outside the support of tau 

    V = leung_isometry()
    bits = [format(i,"04b") for i in range(16)]
    E = np.array([noisy_kraus_operator(bit, gamma)@V for bit in bits])
    tau = encoded_noisy_channel(sigma, gamma)
    sqrt_sigma = psd_sqrt(sigma)
    inv_sqrt_tau = pseudo_inv_sqrt(tau)
    return np.array([sqrt_sigma @ E[i].conj().T @ inv_sqrt_tau for i in range(len(E))])