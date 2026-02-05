# Tunnel Karstic Inflow Estimation

## Overview

This repository presents a stochastic modelling framework to estimate groundwater inflows into tunnels excavated in karstified carbonate rocks.  
The codes provided here support the methodology developed in a scientific study. Raw input data are not publicly available.

## Related publication

Dallalba et al. (2022), *Journal of Hydrology*.  
https://www.sciencedirect.com/science/article/pii/S0013795222004355

## Methodological workflow

![Methodology workflow](images/methodology.png)

The approach is based on four main steps:

1. **3D karst network simulation**  
   Stochastic generation of karst conduit networks using the **Pykaso** package.

2. **Conduit aperture simulation**  
   Aperture fields derived from Gaussian Random Fields (GRF) using the geostatistical package **Geone**, combined with statistical distributions from the literature.

3. **3D flow simulation**  
   Hydraulic simulations of flow within the karst network and tunnel–karst interactions performed with the proprietary software **DISCO**.

4. **Stochastic analysis of tunnel inflows**  
   Monte Carlo simulations to explore the full range of possible breakthrough discharges.

An a posteriori analysis of conduit friction coefficients is performed and compared with analytical formulations and literature values.

## Example of results

![Simulation results](images/results.png)

## Repository structure

Tunnel_karstic_inflow_estimation/
├── 00_test_bernouilli_weisbach/
├── 01_test_karst_network_generator/
└── 02_final_results/


## Notes

- Raw geological and tunnel data are not included.
- Codes are provided to support the methodology and may require adaptation for other applications.

