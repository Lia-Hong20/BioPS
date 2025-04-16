# BioPS

**BioPS** (Biological Phase Separation Simulator) is a simulation tool designed to model the coupling between stochastic biochemical reactions and liquid-liquid phase separation (LLPS) in biological systems. In **BioPS**, biochemical reaction products can induce spatial reorganization by triggering phase separation, while the resulting phase-separated compartments, in turn, modulate the local kinetics of biochemical reactions—capturing the regulatory role of phase separation in gene expression within cells.

## Introduction

This repository provides a comprehensive suite of Julia functions for simulating coupled LLPS–reaction systems at the mesoscopic scale. The tools in this repository enable users to:

- Simulate the influence of stochastic biochemical reactions on phase behavior.
- Investigate the feedback effects of phase-separated compartments on the dynamics of reaction networks.
- Model non-equilibrium phase separation arising from the competition between biochemical reaction timescales and condensate relaxation dynamics.
- Capture in detail the interplay between protein production, degradation, and polymer chain dynamics.
- Flexibly configure reaction networks, polymer properties, and coupling parameters via parameterized input files.

These functions are designed to facilitate the study of intracellular phase behavior, spatial-temporal feedback, and nonequilibrium biophysical modeling.

## Environment Setup

**BioPS** can be run locally or on a server. Follow these steps to set up the environment:

1. **Download the repository**  
   Click the green "Code" button on the top-right of this page and choose "Download ZIP" to get the repository. Extract it to a working directory.

2. **Load the scripts**  
   Make sure all `.jl` scripts are in the same directory or properly referenced. The main script (`BioPS.jl`) will include them accordingly.

## Running Simulations, Analysis, and Visualization

To run the full simulation, execute `BioPS.jl`. This will demonstrate the evolution of a biochemical system coupled with LLPS, generating time-resolved outputs and visualizations.

### Required Input Files

- **polymer_topology.txt**: Defines the connectivity and sequence of polymer chains.
- **interaction_matrix.txt**: Specifies the interaction energies between bead types.
- **constants.jl**: Sets fixed simulation parameters such as system size and total steps.
- **parameters.jl**: Describes biochemical network structure, reaction rates, and timescales.

### Output Files

- **simulation_results.xlsx**: Contains system-level observables (e.g., total energy, largest cluster size), polymer configurations, condensate identities, and metadata (e.g., parameters).
- **simulation_results.eps**: Exports plots visualizing time-evolving system behavior and 3D bead positions.

## Script Structure and Functional Description

The repository is modular, with each script supporting specific tasks. Descriptions are as follows:

- **parameter_initialization.jl**  
  Manages loading and initialization of system parameters. Key responsibilities include:
  1. Configuring and validating simulation parameters.
  2. Initializing the spatial configuration of polymer chains.

- **data_initialization.jl**  
  Parses user-defined input files, validates data consistency, and adjusts the number and spatial distribution of polymers according to the specified initial conditions to match the lattice size and physical constraints.

- **system_initialization.jl**  
  Constructs the lattice model using validated parameters, places polymer chains in 3D space, initializes the interaction list, and calculates the initial energy state of the system.

- **simulation.jl**  
  The main driver script for the coupled dynamics simulation. It coordinates the execution of both biochemical reactions and phase separation in an iterative framework, ensuring their temporal and spatial integration.

- **biochemical_reaction.jl**  
  Implements stochastic biochemical reaction simulation based on the Gillespie algorithm. It dynamically inserts or removes polymer chains in the system to simulate synthesis or degradation events.

- **phase_separation.jl**  
  The main simulation function for phase separation processes.

- **move_behavior.jl**  
  Implements motion patterns for LLPS, defining rules for the dynamics of beads, chain segments, and clusters. This includes monomer diffusion, segmental rotation, and cluster-level translation, forming the move set for Monte Carlo sampling.

- **analysis.jl**  
  Analyzes simulation results by identifying the largest cluster and defining droplets (e.g., clusters containing more than 10 chains). It also computes statistical characteristics such as pair distribution functions, mean squared displacement, and radius of gyration to reveal spatial correlation features of the system.

- **visualization.jl**  
  Translates simulation outputs into visual graphics, supporting the generation of polymer spatial maps, energy evolution curves, and cluster size trajectories to help users intuitively interpret simulation results.

## Example

An example of the Birth-Death process is provided. The input files for this example can be found in the folder **"BirthDeath_example"**. These files follow the same structure as above but are tailored for simulating the Birth-Death process.

## Citation

If you find this code useful for your research, please cite it as follows:

XXXXXXX. (2025). XXXXXXXXXXXXXXXXX in Julia. GitHub repository, [https://github.com/ZhangLab/BioPS](https://github.com/ZhangLab/BioPS).
