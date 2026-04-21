# Monte Carlo Simulation of Neutron Slowing Down
## Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium

**Author:** Ching Kai Sing, Lucas

**Based on Coveyou et al. (1956) "Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium
" in Oak Ridge Laboratory**

![Repo Languages](https://github-readme-stats.vercel.app/api/top-langs/?username=lucas-cks&layout=compact&hide=html,css&repo=Neutron-Slowing-Monte-Carlo)(https://github.com/anuraghazra/github-readme-stats)

![C](https://img.shields.io/badge/Language-C-blue?logo=c)
![Python](https://img.shields.io/badge/Language-Python-yellow?logo=python)

![NumPy](https://img.shields.io/badge/Library-NumPy-013243?logo=numpy&logoColor=white)
![Pandas](https://img.shields.io/badge/Library-Pandas-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Library-Matplotlib-ffffff?logo=matplotlib&logoColor=black)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## 1. Overview
This project implements a high-fidelity Monte Carlo simulation of neutron moderation in an infinite, capturing medium. It is based on the classic 1956 research by **Coveyou, Bate, and Osborn** at the Oak Ridge National Laboratory[cite: 1]. The simulation explores the effects of moderator temperature on neutron flux spectra, validating results against **Wigner-Wilkins (1944)** theory[cite: 2]. For a high-level overview of the physics and methodology, please see my presentation slides.

[Back to Top](#readme-top)


## 2. Physics Background
The simulation models neutron thermalization with the following parameters:
* **Target Motion:** Moderator atoms follow a Maxwell-Boltzmann distribution.
* **Cross-sections:** $1/v$ absorption cross-section is assumed.
* **Scattering:** Elastic scattering with temperature-dependent target velocity sampling.
* **Reproduction of Key Phenomena:**
    * **Spectral Hardening:** The thermal peak shifts to higher energies as absorption ($K$) increases.
    * **1/v Tail:** High-speed flux follows the theoretical $1/v$ behavior.
    * **Temperature Shift:** Deviation from the ideal Maxwell-Boltzmann distribution.

[Back to Top](#readme-top)

## 3. Implementation Details
* **Language:** C (for simulation) and Python (for data analysis).
* **Scale:** 1,000,000 neutron histories per case to ensure low statistical noise.
* **Sampling:** Implements rejection sampling for target velocities as described in the 1956 paper.
* **Normalization:** Tallies are correctly normalized to flux for comparison with analytical models.

[Back to Top](#readme-top)

## 4. Repository Structure
* `monte_carlo_neutron_slowing.c`: Main simulation engine in C.
* `monte_carlo_neutron_slowing.exe`: Main simulation engine in exe.
* `plot.py`: Python script for visualization and statistical analysis.
* `data/`: CSV files containing simulation results for Hydrogen, Carbon, and High Absorption cases.
* `result/`: Generated plots showing flux comparisons and hardening factors.
* `docs/`: Original reference papers by Coveyou et al., Wigner-Wilkins, and presentation slides.
* `README.md`
* `LICENSE`: MIT License

[Back to Top](#readme-top)

## 5. Code Structure
```text
main()
├──init_simulation()        
├──run_simulation()        
│  └──run_history()     
|    ├──calculate_effective_sigma_s()
│    ├──sigma_a(v) = σa0 / v
│    ├──sample_target_velocity()
│    ├──elastic_collision()
│    └──tally speed into bins
└──calculate_flux()
```      
[Back to Top](#readme-top)

## 6. Usage
### Installation
```bash
git clone [https://github.com/your_username/neutron-slowing-mc.git](https://github.com/your_username/neutron-slowing-mc.git)
```

### Compilation
```bash
gcc -o monte_carlo_neutron_slowing monte_carlo_neutron_slowing.c -lm -O3
```

### Execution
```bash
./monte_carlo_neutron_slowing
python plot.py
```

Requirements

- **C compiler**
- **Python 3.x** with:
  - `numpy`, `pandas`, `matplotlib`

Install Python dependencies:
```bash
pip install numpy pandas matplotlib
```
[Back to Top](#readme-top)

## 7. Key Results & Validation
The simulation reproduces the theoretical relationship between moderator temperature ($T_m$) and effective neutron temperature ($T_e$):
$$T_m/T_e = 1 + 1.11 \times A \times K$$

| Case | Moderator | A | K | Peak Position ($y$) |
| :--- | :--- | :--- | :--- | :--- |
| 1 | Hydrogen | 1 | 0.18 | $\approx 1.094$ |
| 2 | Carbon | 12 | 0.10 | $\approx 1.406$ |
| 3 | High Absorption | 1 | 0.36 | $\approx 1.219$ |

See result/ for corresponding flux distribution plots.

[Back to Top](#readme-top)

## 8. License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## 9. References
1. Coveyou, R. R., Bate, R. R., & Osborn, R. K. (1956). "Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium". *Journal of Nuclear Energy*. (PDF available in `docs/`)
2. Wigner, E. P., & Wilkins, J. E. (1944). "Effect of the Temperature of the Moderator on the Velocity Distribution of Neutrons". *AECD-2275*. (PDF available in `docs/`)

## 10. Contact

For questions or suggestions, please open an issue on this repository or contact the author directly.

**Ching Kai Sing, Lucas**  
Department of Physics, The Chinese University of Hong Kong  
*Project Link:* [https://github.com/lucas-cks/Monte-Carlo-Simulation-of-Neutron-Slowing-Down](https://github.com/lucas-cks/Monte-Carlo-Simulation-of-Neutron-Slowing-Down)

[Back to Top](#readme-top)
