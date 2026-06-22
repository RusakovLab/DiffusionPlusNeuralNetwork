# Neurotransmitter Spillover Shapes Network Computation: Monte-Carlo Diffusion and Hopfield Attractor Simulations

Leonid Savtchenko & Dmitri Rusakov  
Institute of Neurology, University College London

---

## Overview

This repository contains simulation code for the accompanying manuscript. The work addresses two interconnected questions:

1. **How does glutamate escape the synaptic cleft and interact with perisynaptic astroglia?** — modelled via 3D Monte-Carlo particle diffusion through a realistic neuropil geometry.
2. **How does neurotransmitter spillover affect associative memory recall?** — modelled via Hopfield-type attractor networks augmented with a spatially graded "spillover" lateral field.

All simulations are implemented in MATLAB (R2018b or later).

---

## Repository Structure

```
├── Diffusion/                    # Monte-Carlo glutamate diffusion model
│   ├── Fig_SAN2000Psi01A01.m     # Main simulation script
│   ├── InputParametersSR.m       # Parameter loader (reads statisticSR.txt)
│   ├── NMDA_SpaceSR_New.m        # NMDA receptor occupancy analysis
│   └── statisticSR.txt           # User-editable configuration file
│
├── NetworkCode/                  # Hopfield attractor network models
│   ├── Fig_Hopfield_LIF_Sparse.m              # LIF spiking Hopfield network
│   ├── Fig_Sup_Hopfield_GeneralizedNetwork.m   # Generalised Hopfield network
│   ├── Fig_Sup_Hopfield_LIF.m                  # LIF validation model
│   ├── Fig_FG_Hopfield_H_V.m                   # Paper panels F & G (binary Hopfield + spillover kernel)
│   └── Fighopfieldlifspillover.m               # Spiking LIF Hopfield with biophysical spillover reformulation
│
├── data/                         # Output data and intermediate results
└── README.md
```

---

## Part 1 — Monte-Carlo Glutamate Diffusion (`Diffusion/`)

### Model Description

The simulation environment is a 4 µm³ cube containing:

- A central synaptic cleft (~220 nm diameter, ~20 nm gap).
- Surrounding neuronal and astrocytic processes represented as randomly sized, overlapping spheres (radii 50–350 nm), consistent with EM-derived neuropil ultrastructure.
- Sphere identity (astrocyte vs. neuron) is assigned stochastically according to the parameter `ProbabilityofAstrocytes`.

Glutamate molecules undergo Brownian motion with free-space diffusivity *D* = 0.5 µm²/ms. Binding to astroglial surfaces follows a time-dependent probability:

$$P(t) = 1 - e^{-t / \Psi}$$

where Ψ (`TimeINsideAdhesiveZone`) controls the build-up of binding probability during membrane contact.

### Files

| File | Purpose |
|------|---------|
| `Fig_SAN2000Psi01A01.m` | Main script: generates neuropil geometry, releases glutamate, simulates diffusion with astroglial binding, writes positional snapshots and diffusion statistics. |
| `InputParametersSR.m` | Reads key–value pairs from `statisticSR.txt` and returns structured parameter sets (geometry, diffusion constants, binding kinetics). |
| `NMDA_SpaceSR_New.m` | Post-processing script: converts radial glutamate distributions into spatiotemporal concentration maps and computes NMDA receptor open-state occupancy via ODE integration. |
| `statisticSR.txt` | Configuration file for particle count, astrocyte fraction, binding parameters, and sphere packing density. |

### Running the Diffusion Simulation

1. Edit `statisticSR.txt` to set desired parameters:
   ```
   np                      = 2000
   Trials                  = 10
   MaxProbAdhesive         = 1.0
   TimeINsideAdhesiveZone  = 0.1
   ProbabilityofAstrocytes = 0.10
   Numbersphere            = 1512
   UnboundProb             = 0.0
   ```

2. Run the main script in MATLAB:
   ```matlab
   Fig_SAN2000Psi01A01
   ```

3. Outputs include:
   - Molecular position snapshots at specified time points (`PD_*.txt`)
   - Radial distribution histograms for bound and free particles (`DistanceBound*.txt`, `DistanceFree*.txt`)
   - Diffusion coefficient traces (`Diffusion_Set_ms_*.txt`)

4. To compute NMDA receptor activation from diffusion output:
   ```matlab
   NMDA_SpaceSR_New
   ```
   This produces spatiotemporal maps of glutamate concentration and NMDA open-state probability as a function of distance and time.

### Key Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `np` | 2000 | Number of glutamate molecules released |
| `ProbabilityofAstrocytes` | 0.10 | Fraction of spheres representing astroglial processes |
| `TimeINsideAdhesiveZone` (Ψ) | 0.1 ms | Time constant for astroglial binding probability |
| `MaxProbAdhesive` | 1.0 | Maximum per-step binding probability |
| `Numbersphere` | 1512 | Total spheres defining the neuropil volume fraction |
| `TotalTime` | 8 ms | Simulation duration |
| `DifCoefINIT` | 5 × 10⁵ nm²/ms | Free-space glutamate diffusion coefficient |

---

## Part 2 — Hopfield Attractor Networks with Spillover (`NetworkCode/`)

### Model Description

We implement Hopfield associative memory networks augmented with a spatially graded lateral interaction term ("spillover field") that models the extrasynaptic influence of neurotransmitter diffusion on network computation. The spillover is realised as a signed ring convolution with a Gaussian kernel of width σ and range *M*.

Retrieval quality is quantified as:

$$Q(\eta) = \frac{1 + m(\eta)}{2}$$

where *m*(η) is the overlap between the retrieved state and the target pattern after presenting a cue corrupted by noise fraction η.

### Files and Figures

| File | Figure | Description |
|------|--------|-------------|
| `Fig_Hopfield_LIF_Sparse.m` | *Memory recall in a sparse LIF network.* Panel shows recall quality (%) vs. noise level (%) for given *M* and σ. | Sparse spiking Hopfield network with LIF neurons, local connectivity, and biophysical spillover mixing. Reports recall quality vs. noise for varying *M* and σ. |
| `Fig_Sup_Hopfield_GeneralizedNetwork.m` | *Recall stability in the generalised Hopfield network.* Panels show retrieval quality *Q*, overlap *m*, and precision vs. noise η for multiple σ values. | Generalised (non-spiking) Hopfield network with signed spillover convolution. Demonstrates equivalence of attractor retrieval dynamics across model formulations. |
| `Fig_Sup_Hopfield_LIF.m` | *Noise-dependent retrieval quality in the LIF Hopfield network.* Panel shows *Q*(η) curves for different σ, with AUC and critical noise analysis. | Updated LIF implementation confirming reproducibility. Includes noise robustness analysis (AUC, critical noise η_c). |
| `Fig_FG_Hopfield_H_V.m` | *Paper panels F and G.* Panel F: noise-free recall (stored pattern → probe → recalled). Panel G: 2×3 grid of recalled patterns at two noise levels (20%, 30%) × three spillover widths (σ = 0, 5, 10). | Self-contained binary (±1) Hopfield network with a 1D Gaussian spillover kernel applied via circular shift. Spillover amplitude is auto-calibrated to the Hebbian drive. Convergence detected by a sliding window; display trial selected by an isotropy filter that excludes spurious striped attractors. Exports `Panel_F.png` and `Panel_G.png`. |
| `Fighopfieldlifspillover.m` | *Quality of recall vs. noise for multiple spillover coefficients M.* Full sweep of noise levels (0–100%, 10% steps) × M values [1, 10, 25, 50, 100, 200] × 100 trials, plus paper panels at 0% and 20% noise. | Biophysical reformulation of the spiking LIF Hopfield network. M is reinterpreted as a **neurotransmitter spillover coefficient**: two fixed synaptic matrices coexist — a sparse wired Hebbian core (`W_core`, fan-in = 3 nearest neighbours) and a spatially gated full Hebbian spillover matrix (`W_spill`, Gaussian footprint of width σ = 5). Retrieval drive is `W_eff = (1 − α_M)·W_core + α_M·W_spill` where `α_M = α_max·(1 − exp(−M/M_tau))`. Pattern activity is decoded from the mean firing rate in the final 5-second window. |

### Running the Network Simulations

**LIF spiking network** (`Fig_Hopfield_LIF_Sparse.m`):
```matlab
Fig_Hopfield_LIF_Sparse
```
Outputs a CSV file (`M=*.csv`) and an errorbar plot of recall quality vs. noise.

**Generalised Hopfield network** (`Fig_Sup_Hopfield_GeneralizedNetwork.m`):
```matlab
params = struct();
params.visualize_patterns = true;
params.save_data_csv = true;
results = run_spiking_hopfield_v2(params);
```

**LIF validation model** (`Fig_Sup_Hopfield_LIF.m`):
```matlab
results = run_spiking_hopfield_LIF_v3();
```

**Paper panels F & G — binary Hopfield with spillover kernel** (`Fig_FG_Hopfield_H_V.m`):
```matlab
Fig_FG_Hopfield_H_V
```
Saves `Panel_F.png` and `Panel_G.png` to the working directory and prints Q mean ± s.d. to the console.

**Spiking LIF with biophysical spillover sweep** (`Fighopfieldlifspillover.m`):
```matlab
run Fighopfieldlifspillover
```
Runs automatically on execution. Set `GENERATE_PAPER_FIGURES = true` (default) to produce the 0% and 20% noise paper panels in addition to the full Q-vs-noise errorbar plot for all M values.

### Key Network Parameters

| Parameter | LIF Sparse | Generalised | Panels F & G | LIF Spillover Sweep | Description |
|-----------|------------|-------------|--------------|---------------------|-------------|
| *N* | 400 | 400 | 400 | 400 | Number of neurons |
| *P* | 3 | 3 | 3 | 3 | Number of stored patterns |
| *M* | 10 | 200 | 200 (fixed) | 1–200 (swept) | Spillover range / coefficient |
| σ | 5 | 0, 2, 5, 10 | 0, 5, 10 | 5 (fixed) | Spatial width of spillover kernel |
| η | 0–1 | 0–1 | 0, 0.2, 0.3 | 0–1 | Noise level (fraction of flipped bits) |
| Trials | 5 | 100 | 50 | 100 | Repetitions per condition |
| Neuron model | LIF spiking | Binary ±1 | Binary ±1 | LIF spiking | Update rule |

---

## System Requirements

### Hardware

- **Processor:** Any modern x86-64 CPU (Intel or AMD); multi-core recommended for parameter sweeps
- **RAM:** ≥ 8 GB (16 GB recommended for large diffusion runs with `np` > 5000)
- **Disk:** ~500 MB free space for code, outputs, and intermediate files
- **GPU:** Not required (all computations are CPU-based)

Typical runtimes on a standard desktop workstation (e.g., Intel i7, 16 GB RAM):
- Hopfield network scripts: 1–30 minutes depending on `trials` and network size
- Monte-Carlo diffusion (`Fig_SAN2000Psi01A01.m`): 1–8 hours per trial depending on `np` and `Numbersphere`
- NMDA post-processing (`NMDA_SpaceSR_New.m`): ~10–30 minutes

### Software Dependencies

| Software | Version | Required | Purpose |
|----------|---------|----------|---------|
| **MATLAB** | R2018b or later (tested on R2022b, R2024a) | Yes | All simulations |
| MATLAB Parallel Computing Toolbox | Any compatible version | No (optional) | Accelerated file processing in `NMDA_SpaceSR_New.m` (vectorised variant) |

No additional MATLAB toolboxes, third-party libraries, Python packages, or compiled binaries are required. All code uses core MATLAB functions only.

### Operating System

The code is platform-independent and runs on any OS supported by MATLAB:
- Windows 10 / 11
- macOS 12 (Monterey) or later
- Linux (Ubuntu 20.04+, CentOS 7+, or equivalent)

### Installation Instructions

1. **Install MATLAB** (R2018b or later) from [mathworks.com](https://www.mathworks.com/products/matlab.html). A standard MATLAB installation with no additional toolboxes is sufficient.

2. **Clone or download this repository:**
   ```bash
   git clone https://github.com/<username>/<repository>.git
   cd <repository>
   ```
   Alternatively, download the ZIP archive from the repository page and extract it.

3. **Verify the installation** by launching MATLAB and running:
   ```matlab
   ver   % Confirm MATLAB version is R2018b or later
   ```

4. **No additional setup is needed.** All scripts are self-contained and can be run directly from the MATLAB command window or editor. See the "Running" sections below for each simulation component.

> **Note:** No compilation, `mex` building, or path configuration is required. Simply navigate to the appropriate directory (`Diffusion/` or `NetworkCode/`) in MATLAB and run the scripts.

---

## Outputs

All scripts generate figures directly in MATLAB. Network simulation scripts additionally produce:

- CSV files with mean ± s.d. of recall quality across noise levels
- `.mat` result files for downstream analysis
- PNG figure exports (when `save_dir` is specified)

The diffusion simulation writes timestamped text files for molecular positions, radial distributions, and diffusion coefficients.

---

## Citation

If you use this code, please cite the accompanying manuscript:

> Savtchenko L, Rusakov D. *[Manuscript title]*. [Journal and DOI to be added upon publication]

---

## License

This code is provided for research reproducibility. Please contact the authors before any commercial or derivative redistribution.

---

## Contact

Prof Dmitri Rusakov — University College London  
Dr Leonid Savtchenko — University College London  
Institute of Neurology, UCL
