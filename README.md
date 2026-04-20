# Monte Carlo Simulation of Neutron Slowing Down

**Author:** Ching Kai Sing, Lucas 

## 1. Overview
This project implements a high-fidelity Monte Carlo simulation of neutron moderation in an infinite, capturing medium. It is based on the classic 1956 research by **Coveyou, Bate, and Osborn** at the Oak Ridge National Laboratory[cite: 1]. The simulation explores the effects of moderator temperature on neutron flux spectra, validating results against **Wigner-Wilkins (1944)** theory[cite: 2].

## 2. Physics Background
The simulation models neutron thermalization with the following parameters:
* **Target Motion:** Moderator atoms follow a Maxwell-Boltzmann distribution.
* **Cross-sections:** $1/v$ absorption cross-section is assumed.
* **Scattering:** Elastic scattering with temperature-dependent target velocity sampling.
* **Reproduction of Key Phenomena:**
    * **Spectral Hardening:** The thermal peak shifts to higher energies as absorption ($K$) increases.
    * **1/v Tail:** High-speed flux follows the theoretical $1/v$ behavior.
    * **Temperature Shift:** Deviation from the ideal Maxwell-Boltzmann distribution.

## 3. Implementation Details
* **Language:** C (for simulation) and Python (for data analysis).
* **Scale:** 1,000,000 neutron histories per case to ensure low statistical noise.
* **Sampling:** Implements rejection sampling for target velocities as described in the 1956 paper.
* **Normalization:** Tallies are correctly normalized to flux for comparison with analytical models.

## 4. Repository Structure
* `monte_carlo_neutron_slowing.c`: Main simulation engine in C.
* `plot.py`: Python script for visualization and statistical analysis.
* `data/`: CSV files containing simulation results for Hydrogen, Carbon, and High Absorption cases.
* `result/`: Generated plots showing flux comparisons and hardening factors.
* `docs/`: Original reference papers by Coveyou et al. and Wigner-Wilkins.

## 5. Usage
### Compilation
```bash
gcc -o monte_carlo_neutron_slowing monte_carlo_neutron_slowing.c -lm -O3
```

### Execution
```bash
./monte_carlo_neutron_slowing
python plot.py
```

## 6. Key Results & Validation
The simulation reproduces the theoretical relationship between moderator temperature ($T_m$) and effective neutron temperature ($T_e$):
$$T_m/T_e = 1 + 1.11 \times A \times K$$

| Case | Moderator | A | K | Peak Position ($y$) |
| :--- | :--- | :--- | :--- | :--- |
| 1 | Hydrogen | 1 | 0.18 | $\approx 1.094$ |
| 2 | Carbon | 12 | 0.10 | $\approx 1.406$ |
| 3 | High Absorption | 1 | 0.36 | $\approx 1.219$ |


## 7. References
1. Coveyou, R. R., Bate, R. R., & Osborn, R. K. (1956). "Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium". *Journal of Nuclear Energy*. 
2. Wigner, E. P., & Wilkins, J. E. (1944). "Effect of the Temperature of the Moderator on the Velocity Distribution of Neutrons". *AECD-2275*.
