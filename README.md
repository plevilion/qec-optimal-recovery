# Quantum Error Correction: Petz vs. SDP-Optimal Recovery

A numerical study of approximate quantum error correction under amplitude-damping noise.

## Overview

This project investigates quantum error recovery for the four-qubit Leung–Nielsen–Chuang–Yamamoto (LNCY) code under amplitude-damping noise.

The main objective is to implement and compare different recovery strategies, particularly the Petz (transpose) recovery map and the optimal recovery channel obtained through semidefinite programming (SDP).

The project draws on established results in approximate quantum error correction and aims to reproduce selected numerical comparisons from the literature, particularly the recent work of Laura Burri (2026).

An additional objective is to explore how recovery performance changes when the actual noise model differs from the one assumed when constructing the recovery map.

## Mathematical Framework

### 1. Quantum Encoding

We consider the four-qubit LNCY code, which encodes one logical qubit into four physical qubits using the codewords

$$
|0_L\rangle = \frac{|0000\rangle+|1111\rangle}{\sqrt{2}},
$$

$$
|1_L\rangle = \frac{|0011\rangle+|1100\rangle}{\sqrt{2}}.
$$

The encoding is described by the isometry

$$
V:\mathbb{C}^2\longrightarrow\mathbb{C}^{16},
\qquad V^\dagger V=I_2,
$$

which induces the quantum channel

$$
\mathcal E(\rho)=V\rho V^\dagger.
$$

### 2. Amplitude-Damping Noise

Each physical qubit is independently subjected to amplitude-damping noise with parameter $\gamma\in[0,1]$.

The single-qubit channel has Kraus operators

$$
A_0=
\begin{pmatrix}
1&0\\
0&\sqrt{1-\gamma}
\end{pmatrix},
\qquad
A_1=
\begin{pmatrix}
0&\sqrt{\gamma}\\
0&0
\end{pmatrix}.
$$

The resulting noisy encoding channel is

$$
\Phi_\gamma=\mathcal N_\gamma^{\otimes 4}\circ\mathcal E.
$$

### 3. Quantum Recovery

We seek a completely positive trace-preserving recovery map

$$
\mathcal R:\mathcal B(\mathbb C^{16})
\longrightarrow\mathcal B(\mathbb C^2)
$$

that reconstructs the original logical state as accurately as possible.

Recovery performance is measured through entanglement fidelity:

$$
F_e(\mathcal R\circ\Phi_\gamma)=
\langle\Phi^+|
(\operatorname{id}\otimes\mathcal R\circ\Phi_\gamma)
(|\Phi^+\rangle\langle\Phi^+|)
|\Phi^+\rangle,
$$

where $|\Phi^+\rangle=(|00\rangle+|11\rangle)/\sqrt2$.

The principal optimization problem is

$$
\boxed{
\max_{\mathcal R\ \mathrm{CPTP}}
F_e(\mathcal R\circ\Phi_\gamma)
}
$$

For fixed encoding and noise, this optimization can be formulated as a semidefinite program using the Choi representation of the recovery channel.

## Methodology

The numerical study is organized into the following stages:

1. **Channel implementation:** Construct the LNCY encoding and amplitude-damping channel using Kraus representations.
2. **Validation:** Verify complete positivity, trace preservation, and analytical entanglement-fidelity identities.
3. **Petz recovery:** Implement the Petz transpose channel for the noisy encoding.
4. **Optimal recovery:** Formulate and solve the corresponding SDP using CVXPY.
5. **Numerical comparison:** Evaluate the entanglement fidelities of the different recovery strategies across a range of damping parameters.
6. **Further exploration:** Investigate the sensitivity of recovery performance to mismatches between assumed and actual noise parameters.

## Implementation

The project uses Python, primarily:

* NumPy and SciPy for numerical linear algebra.
* CVXPY for semidefinite optimization.
* Matplotlib for visualization.
* pytest for numerical and mathematical validation.

The planned repository structure is:

```text
qec-optimal-recovery/
├── README.md
├── requirements.txt
├── src/
│   ├── channels.py
│   ├── codes.py
│   ├── fidelity.py
│   ├── petz.py
│   ├── recovery.py
│   └── sdp.py
├── notebooks/
├── tests/
└── figures/
```

## Installation

Clone the repository and create a Python virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -e .
```

The tests can be run from the repository root using:

```bash
python -m pytest -v
```

## Project Status

In development.

The initial stage focuses on implementing the LNCY encoding and noise channels, establishing numerical validation tests, and constructing the Petz recovery map.

The SDP optimization, numerical comparisons, and exploratory extensions will be added progressively.

## References

1. D. W. Leung, M. A. Nielsen, I. L. Chuang, and Y. Yamamoto, *Approximate quantum error correction can lead to better codes*, Physical Review A 56, 2567 (1997). [DOI](https://doi.org/10.1103/PhysRevA.56.2567)
2. A. S. Fletcher, P. W. Shor, and M. Z. Win, *Optimum quantum error recovery using semidefinite programming*, Physical Review A 75, 012338 (2007). [DOI](https://doi.org/10.1103/PhysRevA.75.012338)
3. L. Burri, *Entanglement fidelity of the Petz decoder for one-shot entanglement transmission*, Journal of Mathematical Physics 67, 102203 (2026). [arXiv](https://arxiv.org/abs/2502.17411)
4. J. M. Renes, *Quantum Error Correction*, lecture notes, ETH Zürich. [Lecture notes](https://people.phys.ethz.ch/~renes/docs/QEC_20240521.pdf)

## Scope and Attribution

This is an independent educational and numerical research project reproducing and exploring established results in quantum error correction.

The LNCY code, Petz recovery, and SDP-based optimal recovery are existing constructions. The project does not claim originality for these methods. Any additional numerical investigations will be identified separately and compared with relevant prior literature.
