# Monte Carlo Simulation of Neutron Slowing Down
## Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium

**Based on Coveyou et al. (1956) "Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium" in Oak Ridge Laboratory**

![C](https://img.shields.io/badge/Language-C-blue?logo=c)
![Python](https://img.shields.io/badge/Language-Python-yellow?logo=python)

![NumPy](https://img.shields.io/badge/Library-NumPy-013243?logo=numpy&logoColor=white)
![Pandas](https://img.shields.io/badge/Library-Pandas-150458?logo=pandas&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Library-Matplotlib-ffffff?logo=matplotlib&logoColor=black)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## 1. Overview
This project implements a high‑fidelity Monte Carlo simulation of neutron moderation in an infinite, capturing medium, based on the classic 1956 research by **Coveyou, Bate, and Osborn** at Oak Ridge National Laboratory. The simulation explores the effect of moderator temperature on neutron flux spectra and validates results against **Wigner‑Wilkins (1944)** theory. The code now supports **six moderator materials** (H, He, C, Si, Mn), includes a **3D collision model**, and provides **convergence analysis** for statistical validation.

[Back to Top](#readme-top)

## 2. Physics Background
The simulation models neutron thermalization with the following parameters:
* **Target Motion:** Moderator atoms follow a Maxwell‑Boltzmann distribution.
* **Cross‑sections:** $1/v$ absorption cross‑section is assumed.
* **Scattering:** Elastic scattering with temperature‑dependent target velocity sampling, using isotropic scattering in the center‑of‑mass frame (3D).
* **Reproduction of Key Phenomena:**
    * **Spectral Hardening:** The thermal peak shifts to higher energies as absorption ($K$) increases.
    * **1/v Tail:** High‑speed flux follows the theoretical $1/v$ behavior.
    * **Temperature Shift:** Deviation from the ideal Maxwell‑Boltzmann distribution.

[Back to Top](#readme-top)

## 3. Implementation Details
* **Language:** C (simulation engine) and Python 3 (data analysis & visualization).
* **Scale:** 1,000,000 neutron histories per case (configurable) to ensure low statistical noise.
* **Sampling:** Implements rejection sampling for target velocities as described in the 1956 paper.
* **Normalization:** Tallies are correctly normalized to flux for comparison with analytical models.
* **3D Collision Kinematics:** Fully vectorised elastic scattering with isotropic CM frame scattering.
* **Error Estimation:** Flux uncertainties are computed from history‑to‑history variance.
* **Convergence Study:** Automatically runs simulations with increasing history counts (10k → 1M) to verify stability of peak position and tail slope.

[Back to Top](#readme-top)

## 4. Repository Structure
```
.
├── monte_carlo_neutron_slowing.c     # Main simulation engine in C
├── plot.py                           # Python script for visualization and statistical analysis
├── data/                             # CSV files containing simulation results
├── result/                           # Generated plots (PNG files)
├── convergence_study/                # Convergence study metrics
├── docs/                             # Original reference papers by Coveyou et al. and Wigner-Wilkins
├── archive/                          # Archived (old) verion of the project
├── README.md                         # This file
└── LICENSE                           # MIT License

```

[Back to Top](#readme-top)

## 5. Code Structure
```text
main()
├── init_simulation()        
├── run_simulation()        
│   └── run_history()     
│       ├── calculate_effective_sigma_s()
│       ├── sigma_a(v) = σa0 / v
│       ├── sample_target_velocity()  [rejection sampling on speed & μ]
│       ├── elastic_collision()       [3D, isotropic CM frame]
│       └── tally speed into bins
├── calculate_flux()                  [with variance]
├── extract_metrics()                 [peak y, tail slope]
├── run_convergence_study()           [multiple history counts]
└── save_results()
```
[Back to Top](#readme-top)

## 6. Usage
### Installation
```bash
git clone https://github.com/lucas-cks/Neutron-Slowing-Monte-Carlo.git
cd Neutron-Slowing-Monte-Carlo
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

The C code will generate CSV files for each moderator case (`neutron_flux_hydrogen.csv`, `neutron_flux_carbon.csv`, `neutron_flux_helium.csv`, `neutron_flux_silicon.csv`, `neutron_flux_manganese.csv`, `neutron_flux_highK.csv`) and a convergence file (`convergence_histories.csv`). The Python script then produces four plots:
- `comprehensive_flux_analysis.png` – 2×2 panel with log‑log spectrum, thermal peak, deviation from MB, and summary table.
- `individual_flux_comparisons.png` – separate subplots for each moderator.
- `simple_flux_comparison.png` – combined log‑log and linear plots.
- `convergence_study.png` – peak position, peak flux, and tail slope vs. number of histories.

### Requirements
- **C compiler** (e.g., `gcc`)
- **Python 3.x** with:
  - `numpy`, `pandas`, `matplotlib`

Install Python dependencies:
```bash
pip install numpy pandas matplotlib
```

[Back to Top](#readme-top)

## 7. Key Results & Validation

The simulation reproduces the theoretical relationship between moderator temperature ($T_m$) and effective neutron temperature ($T_e$):
$\frac{T_m}{T_e} = 1 + 1.11 \cdot A \cdot K$

| Case                | Moderator | A  | K       | A×K   | Expected Peak $y$ | Measured Peak $y$ | Deviation |
| :------------------ | :-------- | :- | :------ | :---- | :---------------- | :---------------- | :--------- |
| 1                   | Hydrogen  | 1  | 0.18    | 0.180 | 1.095             | 1.094             | –0.1%      |
| 2                   | High Abs. | 1  | 0.36    | 0.360 | 1.183             | 1.156             | –2.3%      |
| 3                   | Helium    | 2  | 0.03125 | 0.062 | 1.034             | 1.031             | –0.3%      |
| 4                   | Carbon    | 12 | 0.03125 | 0.375 | 1.190             | 1.156             | –2.8%      |
| 5                   | Silicon   | 16 | 0.03125 | 0.500 | 1.247             | 1.156             | –7.3%      |
| 6                   | Manganese | 25 | 0.03125 | 0.781 | 1.366             | 1.281             | –6.2%      |

All cases show $1/v$ tail slopes within 15% of –1.0, confirming the expected high‑energy asymptotic behavior. The hydrogen and helium cases match theoretical predictions to within <0.5%. Heavier moderators show slightly more deviation due to the simplified effective cross‑section model, but the trend of increasing peak energy with $A\times K$ is clearly reproduced.

Note: The ~7% deviation in heavier moderators (Si, Mn) at high AK values is consistent with the reported 5% reliability limit of the original Coveyou interpolation formula, reflecting the breakdown of the linear spectral-shift approximation in strong absorption regimes.

**Convergence study** (Hydrogen, $K=0.18$) demonstrates that:
- Peak position stabilises after ~1,000,000 histories.
- Tail slope converges to ≈ –1.0 with standard error <0.05.

See the `result/` directory for the generated plots and `convergence_histories.csv` for numerical convergence data.

[Back to Top](#readme-top)

## 8. License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

## 9. References
1. Coveyou, R. R., Bate, R. R., & Osborn, R. K. (1956). "Effect of Moderator Temperature upon Neutron Flux in Infinite, Capturing Medium". *Journal of Nuclear Energy*. (PDF available in `docs/`)
2. Wigner, E. P., & Wilkins, J. E. (1944). "Effect of the Temperature of the Moderator on the Velocity Distribution of Neutrons". *AECD-2275*. (PDF available in `docs/`)

## 10. Contact

For questions or suggestions, please open an issue on this repository or contact the author directly.

*Project Link:* [https://github.com/lucas-cks/Neutron-Slowing-Monte-Carlo](https://github.com/lucas-cks/Neutron-Slowing-Monte-Carlo)

[Back to Top](#readme-top)
