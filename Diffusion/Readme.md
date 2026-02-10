# Monte-Carlo Simulation of Glutamate Diffusion and Astroglial Uptake (SpaceSR Model)

This repository contains MATLAB code implementing a 3D Monte-Carlo model of synaptic glutamate diffusion, extracellular tortuosity, and astrocyte-mediated glutamate uptake via surface binding. The model follows the overlapping-spherical neuropil representation commonly used in quantitative ultrastructural reconstructions of neocortical layer 1.

The simulation reproduces **transient glutamate escape** from synapses and its interactions with neurons and astrocytes under variable transporter densities.



---

## Repository Structure

| File | Description |
|------|-------------|
| **Fig4SAN2000Psi01A01.m** | Main simulation script. Generates stochastic neuropil geometry, releases glutamate, simulates diffusion and astrocyte interactions, and writes output snapshots and statistics. |
| **InputParametersSR.m** | Parameter loader. Reads the simulation settings from `statisticSR.txt` and constructs MATLAB parameter structures. |
| **NMDA_SpaceSR_New.m** | Optional analysis script for computing receptor activation (e.g., NMDA receptor occupancy) from diffusion output. |
| **statisticSR.txt** | User-editable configuration file specifying particle number, astrocyte fraction, binding probability parameters, and neuropil packing density. |

---

## Model Overview

The simulation environment is a cube (typically 4 µm across) containing:

- A central synaptic cleft (≈220 nm × 20 nm).
- Surrounding neuronal and astrocytic processes modeled as **randomly sized, overlapping spheres**.
- Sphere identities determined by `ProbabilityofAstrocytes`.
- Glutamate molecules propagate via Brownian motion with free-space diffusivity of 0.5 µm²/ms.
- Astrocyte surfaces allow **surface binding** according to:
  
\[
P(t) = 1 - e^{-t / \Psi}
\]

where **Ψ** (TimeINsideAdhesiveZone) controls the binding build-up during membrane proximity.

---

## How to Run the Simulation

1. Open MATLAB (R2018b or later recommended).
2. Ensure all `.m` files and `statisticSR.txt` are in the same directory.
3. Edit simulation parameters in:

np                      = 2000
Trials                  = 1
MaxProbAdhesive         = 1.0
TimeINsideAdhesiveZone  = 0.1
ProbabilityofAstrocytes = 0.10
Numbersphere            = 1500
UnboundProb             = 0.0


4. Run the main script:

```matlab
Fig4SAN2000Psi01A01

Output files are generated automatically, including:

Diffusion coefficient traces

Molecular position snapshots at specified time points

Radial distribution histograms of bound/unbound particles

License

This code is distributed for research use only.
Please contact the authors for permission before commercial or derivative redistribution.

Contact

Leonid  Savtchenko
Computation Neuroscience, University College London
Email: savtchenko@yahoo.com


