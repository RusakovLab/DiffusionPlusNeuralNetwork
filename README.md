# NatNeurosciSub
Network simulation with synaptic spillover (MATLAB)

This repository contains MATLAB code used to generate all simulations and figures for a network model including synaptic spillover / extrasynaptic signaling, accompanying a manuscript submitted to Nature Neuroscience.

The model implements:

Network dynamics with conventional synaptic coupling

A phenomenological spillover term representing diffusive neurotransmitter signaling

Fixed random seeds for full reproducibility

Usage

All simulations are run from MATLAB by executing the figure-specific scripts:

run simulations/figureX_spillover.m


Each script reproduces one figure or panel in the manuscript and saves outputs to the results/ directory.

Reproducibility

MATLAB version: R2021b or later

All parameters and random seeds are explicitly defined in the scripts

No external toolboxes are required

Code availability

The exact version of the code used in the manuscript corresponds to a tagged release of this repository.
