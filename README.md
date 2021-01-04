# CVOFLS-DPD
Coupled Volume of Fluid, Level Set and Dissipative Particle Dynamics Solver

Lead developer: Raihan Tayeb, Doctoral Candidate, University of Missouri-Columbia

## Overview
A software for mesoscopic simulation developed using open-sourced softwares, [OpenFOAM](https://openfoam.org), [OpenSMOKE++](https://www.opensmokepp.polimi.it), [LIGGGHTS](https://www.cfdem.com/liggghtsr-open-source-discrete-element-method-particle-simulation-code) and [CFDEM](https://www.cfdem.com/cfdemrcoupling-open-source-cfd-dem-framework)

It is employed to investigate the evaporation induced self-assembly of charged polymeric nanoparticles in microdroplet solution

The software can handle
* Complex geometric configuration
* Moving mesh
* Spurious current reduction
* Contact line pinning
* Multi-component liquids (and gases) with mixing and phase change (evaporation and condensation)
* Chemkin style property tables
* Fluid-particle interaction (including fluid interface forces on particles)
* Particle-particle interaction (including DLVO forces and friction force)
* Heat transfers between particles and between particles and fluids
* Parallel processing
