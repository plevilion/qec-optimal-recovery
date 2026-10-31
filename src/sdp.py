import numpy as np
import cvxpy as cp

from channels import encoded_noisy_channel
from petz import completed_petz_kraus



def logic_op(i, j):
    
    op = np.zeros((2,2), dtype=complex)
    op[i][j] = 1
    return op


def objective_matrix(gamma):

    B = [[encoded_noisy_channel(logic_op(0, j), gamma) for j in range(2)], [encoded_noisy_channel(logic_op(1, j), gamma) for j in range(2)]]
    return sum(np.kron(B[i][j].T, logic_op(i, j).T) for i in range(2) for j in range (2))/4

   
def optimal_recovery(gamma):

    C = objective_matrix(gamma)
    J = cp.Variable((32,32), hermitian = True)
    constraints = [J >> 0, cp.partial_trace(J, (16, 2), axis=1) == np.eye(16)]
    objective = cp.Maximize(cp.real(cp.trace(C @ J)))
    problem = cp.Problem(objective, constraints)
    fidelity = problem.solve(solver=cp.SCS, eps=1e-6, max_iters=30000)
    return float(fidelity), J.value



def vec_map(A):
    d = A.shape[1]
    omega = np.eye(d).flatten()
    return np.kron(np.eye(d), A) @ omega



def choi_petz_from_kraus(gamma):

    R = completed_petz_kraus(np.eye(2)/2, gamma)
    vec_R = [vec_map(r) for r in R]
    return sum(np.outer(vec, vec.conj()) for vec in vec_R )