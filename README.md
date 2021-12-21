# PHT OpenFOAM Lab53 IACP FEB RAS

OpenFOAM solvers for simulating heat transfer processes in porous media (porous heat transfer, PHT)

# Overview

In project **PHT OpenFOAM Lab53 IACP FEB RAS**, we develop the OpenFOAM solvers and libraries to simulate the transient fluid flow and heat transfer in porous media. The solvers are developed in the Laboratory of Fluid Mechanics of [Institute of Automation and Control Processes FEB RAS](www.iacp.dvo.ru) and [Far Eastern Federal University](www.dvfu.ru).

**phtEnthalpyFoam** &mdash; fluid flow in porous media without phase changes and chemical reactions. This solver allows one to simulate flows through the "inert" porous materials like rocks, metallic balls, etc. as in sensible heat thermal energy storage or flow through the self-heating porous media.

**porousHeatTransferFoam** &mdash; fluid flow in porous media. Optional modeling of phase transitions in the condenced phase like granular PCM. This solver allow to model flows in packed bed thermal energy storages with phase changes.

# Description

The solvers were developed on [OpenFOAM from OpenCFD](www.openfoam.com) v2012 and tested for versions v1906, v1912, v2006, v2012. Please do the following steps to compile:

1. Download ZIP this project from Github and move it to your $HOME/OpenFOAM Ubuntu directory. 

2. Unzip the archive by running "unzip pht-openfoam-lab53iacp-maser" and run "mv pht-openfoam-lab53iacp-maser OpenFOAM-PHT53IACP" in the prompt.

3. Run "echo 'source $HOME/OpenFOAM/OpenFOAM-PHT53IACP/etc/bashrc' >> ~/.bashrc" in the prompt.

4. Run ./Allwmake script.

# Usage
The project **PHT OpenFOAM Lab53 IACP FEB RAS** is released under [Creative Commons Attribution-NonCommercial 4.0 International Public License](https://creativecommons.org/licenses/by-nc/4.0/legalcode) (see file LICENSE in the main directory or www.creativecommons.org/licenses/by-nc/4.0/legalcode). The project **PHT OpenFOAM Lab53 IACP FEB RAS** is released in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please feel free to extend, modify and use our solvers and libraries according to your requirement, but only for non-commerical usage. When using this project, please cite the project as follows:

For solver **phtEnthalpyFoam**, cite work [1].

For solver **porousHeatTransferFoam**, cite work [1].

# References

1. S. S. Fetsov, A. S. Kim. On a novel finite-volume solver for OpenFOAM for simulating hydrodynamic and thermal processes in porous media // XVII All-Russian Seminar with International Participation "Dynamics of Multiphase Media". Abstracts. Novosibirsk: ITAM SB RAS. 2021. P. 55. (in Russian)
