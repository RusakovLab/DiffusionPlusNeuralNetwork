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
│   └── Fig_Sup_Hopfield_LIF.m                  # LIF validation model
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

### Key Network Parameters

| Parameter | LIF Sparse | Generalised | Description |
|-----------|------------|-------------|-------------|
| *N* | 400 | 400 | Number of neurons |
| *P* | 3 | 3 | Number of stored patterns |
| *M* | 10 | 200 | Spillover range (connections) |
| σ | 5 | 0, 2, 5, 10 | Spatial width of spillover kernel |
| η | 0–1 | 0–1 | Noise level (fraction of flipped bits) |
| Trials | 5 | 100 | Repetitions per condition |

---

## Requirements

- **MATLAB** R2018b or later (no additional toolboxes required for core simulations)
- **Parallel Computing Toolbox** (optional, for accelerated diffusion file processing in `NMDA_SpaceSR_New.m`)

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
