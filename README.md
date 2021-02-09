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

If you use this solver in a project or scholarly work, please include the following citation, [Tayeb and Zhang (2019)](https://doi.org/10.1115/IMECE2019-11953). 

## Installation
The current version of the code uses the [OpenFOAM 2.3.1 libraries](http://www.openfoam.org/archive/2.3.1/download/source.php). It uses [isoadvector](https://github.com/isoAdvector/isoAdvector) library for interface tarcking and adection. It also calculates the thermodynamics and transport properties of gases (diffusion coefficients, thermal conductivity, heat capacities, and viscosity) based on the correlations available in the [OpenSMOKE++](https://www.opensmokepp.polimi.it) library. The coupling with the particle solver which is [LIGGGHTS](https://www.cfdem.com/liggghtsr-open-source-discrete-element-method-particle-simulation-code) is done using [CFDEM](https://www.cfdem.com/cfdemrcoupling-open-source-cfd-dem-framework) library. Once all the libraries are installed, you can proceed to install the CVOFLSDPD solver. 

**To Install:**
Navigate to a working folder in a shell terminal, clone the git code repository, and build.
```
$ git clone https://github.com/rtyme14/CVOFLS-DPD.git CVOFLS-DPD
$ cd CVOFLS-DPD/CVOFLS-DPD
$ source Allwmake.sh
```

The solver can be validated using the case in the [`validation_case`](validation_case). The nanoparticle-fluid case is in the [`colloid_case`](colloid_case).  

## Microgravity droplet evaporation case ([`validation_case`](validation_case))

This case demonstrates evaporation of a spherical micro-droplet,hanging from a thread, in ambient air. The gravity is ignored. Temperature of droplet continues to drop until it reaches wet bulb temperature. Except for a small initial phase, the dimensionless squared diameter follows the linear D2 law. More in the paper.

<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/validation_case/results/Figure3.jpg" width="450" height="400" />

To run the case do the following

```
$ comment out the particle coupling part in CVOFLSDPD.C file
$ wclean && wmake
$ cd ../../validation_case
$ blockMesh
$ setFields
$ decomposePar -force
$ mpirun -np ${no. of processors} CVOFLSDPD -parallel
```
[<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/validation_case/results/Figure4.jpg" width="1000" height="300" />](validation_case/results/Figure4.jpg)

[<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/validation_case/results/Figure5.jpg" width="1000" height="400" />](validation_case/results/Figure5.jpg)

## Nanoparticle self-assembly case ([`colloid_case`](colloid_case))
The self-assembly of nanoparticles inside a liquid micro-droplet is studied using the CVOFLSDPD solver. Be sure to uncomment the particle coupling part and recompile, if it is commented out in the validation case.
There are two folders for the case. The `CFD` folder has `0`, `constant`, `system` and `newdir` just like other case. The `DEM` folder has `in.colloid` file which initializes and solve the particle trajectories. The `Sphere2370FromDump.csv` file containes the initial coordinates for the 2370 nanoparticles. The case can be run in the same way as before using `decomposePar -force && mpirun -np ${no. of processors} CVOFLSDPD -parallel` commands. During simulation, a `liggghts.restart` file will appear which can be used to restart the simulation. This file should be tranferred to the `DEM` folder and in `in.colloid` file the `read_data` command should be commented out and `read_restart` should be uncommented.

The DPD simulation requires much smaller time steps than the CFD process. The time step for the DPD case is 1×10-9 s whereas the time step for the CFD case is 1×10-7 s which means that for each CFD iteration DPD runs for 100 iterations.

The pinned contact angle evaporation and the colloid deposition for three contact angles are shown in the [figure](colloid_case/results/Picture1.jpg) below.

[<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/colloid_case/results/Picture1.jpg" width="1000" height="600" />](colloid_case/results/Picture1.jpg)

The following videos show colloid deposition for 90 degree contact angle.
<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/colloid_case/results/Bottomviewbest.gif" />
<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/colloid_case/results/Sideviewbest.gif" />

## Algorithm
An overview of the algorithm used in CVOFLS-DPD solver is given below.

[<img src="https://github.com/rtymea14/CVOFLS-DPD/blob/main/CVOFLS_DPD/Figure2.jpg" width="400" height="600" />](CVOFLS_DPD/Figure2.jpg)
  
## Contribute
Open to collaboration with other investigators studying phase-change and multiscale particle-fluid flows. Please [contact us](mailto:rthvc@umsystem.edu) if you are interested in expanding the solver or find bugs to correct. Limited support (on a case-by-case basis) or consulting servies can also be provided.

## Acknowledgements
This research was generously supported by the U.S. National Science Foundation.

## References
* [Tayeb, R, Zhang, Y. "Controlling Evaporation Induced Self-Assembly of Polymeric Nanoparticles: A VOF-DPD Study." Proceedings of the ASME 2019 International Mechanical Engineering Congress and Exposition. Volume 8: Heat Transfer and Thermal Engineering. Salt Lake City, Utah, USA. November 11–14, 2019. V008T09A083. ASME.](https://doi.org/10.1115/IMECE2019-11953)
