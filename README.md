# Multi-Mode Pinching-Antenna Systems: Mode Selection or Mode Combining? 

This repository provides the **MATLAB implementation** for the paper: **Multi-Mode Pinching-Antenna Systems: Mode Selection or Mode Combining?** 
[[arXiv Preprint](https://arxiv.org/abs/2603.08472)] 

The code reproduces the simulation results presented in the paper.



## ✨ Key Idea

### Concept of Multi-Mode PASS

>Pinching-Antenna Systems (PASS) have recently emerged as a promising architecture for flexible and efficient wireless communications. 
>This work proposes **new operating protocols** and optimization methods for **multi-mode PASS**. 
>- **Multi-mode PASS** allows a single waveguide to excite **multiple propagation modes** to simultaneously serve multiple users. 
>- This enables **mode-domain multiplexing** and offers extra **degrees of freedom (DoFs)**.
>- Compared to conventional PASS, where each waveguide can support only a *single* independent data stream, multi-mode PASS significantly improves spectral efficiency and resource utilization.



### 1. Operating Protocols

Two operating protocols are studied:

- **Mode Selection**: 
Each pinching antenna (PA) performs phase matching with only a single mode, thus maximizing the coupling strength and radiation power of this selected mode.
In this case, the propagation constant $\beta_{n}^{\mathrm{PA}}$ of a PA is selected from a discrete set $\\{\beta_{1},\beta_{2},\ldots,\beta_{M}\\}$ of the propagation constants of guided modes.

- **Mode Combining**: 
Each PA can flexbily radiate power of multiple modes without phase matching to a specific mode, thereby fully exploiting mode-domain multiplexing.
This is achieved by continuously tuning the propagation constant $\beta_{n}^{\mathrm{PA}}$ of each PA.

> **Uniform Mode Combining**: A practical operating protocol is uniform mode combining,
> where the propagation constant of each PA can be preconfigued at $\beta=(\beta_{1}+\beta_{2}+...\beta_{M})/M$. Our simulation results demonstrate the efficiency of this design.

### 2. Proposed Algorithm

The paper proposes the **Particle Swarm Optimization with KKT Parameterized Beamforming (PSO-KPBF) Algorithm** to jointly optimize the digital beamforming, pinching locations, and PA propagation constants. 
**KPBF** reconstructs KKT-conditioned beamforming solutions of WMMSE problem, which is parameterized by dual varaiables. 
Then, PSO jointly predicts KPBF dual parameters, pinching locations, and propagation constants of PAs. 
> **Benefit**:
>- Reconstructing stationary beamforming solutions without WMMSE iterations in a low-complexity way.
>- Guiding black-box swarm search by KKT solutions, significantly reducing the searching space. 



## 🚀 Reproducing the Results

The repository includes scripts to reproduce the main numerical results from the paper.

### Reproduce Fig. 2 - Achievable Rate vs. Maximum Transmit Power

Simply run [M_PASS_PSO_KPBF_Pmax.m](./M_PASS_PSO_KPBF_Pmax.m)

### Reproduce Fig. 3 - Achievable Rate vs. Number of PAs/Antennas

Simply run [M_PASS_PSO_KPBF_vs_N.m](./M_PASS_PSO_KPBF_vs_N.m)

---

## ❤️ Citation
If you use this code in your research, please cite the following paper: 

X. Xu, X. Mu, Y. Liu, and A. Nallanathan,
“Multi-Mode Pinching-Antenna Systems: Mode Selection or Mode Combining?”
arXiv:2603.08472, 2026.

## 📄 Related Work
This work is based on our previous paper ``Multi-Mode Pinching Antenna Systems Enabled Multi-User Communications'' [arXiv Preprint](https://arxiv.org/abs/2601.20780), 
which establishes the fundamental physic model of multi-mode PASS. In this previous work, the PAs are divided into two groups for fixed mode selection. 
