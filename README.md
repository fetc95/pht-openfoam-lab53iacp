# PHT OpenFOAM Lab53 IACP FEB RAS

OpenFOAM solvers for simulating heat transfer processes in porous media (porous heat transfer, PHT)

# Overview

In project **PHT OpenFOAM Lab53 IACP FEB RAS**, we develop the OpenFOAM solvers and libraries to simulate the transient fluid flow and heat transfer in porous media. The solvers are developed in the Laboratory of Fluid Mechanics of [Institute of Automation and Control Processes FEB RAS](https://iacp.dvo.ru/) and (earlier) [Far Eastern Federal University](https://dvfu.ru/).

**phtEnthalpyFoam** &mdash; fluid flow in porous media without phase changes and chemical reactions. This solver allows one to simulate flows through the "inert" porous materials like rocks, metallic balls, etc. as in sensible heat thermal energy storage or flow through the self-heating porous media.

**porousHeatTransferFoam** &mdash; fluid flow in porous media. Optional modeling of phase transitions in the condenced phase like granular PCM. This solver allows one to model flows in packed bed thermal energy storages with phase changes.

**phtMetalExtractionFoam** &mdash; gas flow through a four-component reacting solid porous phase consisting of an inert, a combustible, a melting/evaporating component and an optional solid combustion product.

# Description

The solvers were developed on [OpenFOAM from OpenCFD](www.openfoam.com) and tested for version v2412. Please do the following steps to compile:

1. Download ZIP this project from Github and, through the prompt, move it to your $HOME/OpenFOAM Ubuntu directory. 
2. Unzip the archive by running "unzip pht-openfoam-lab53iacp-master.zip" and run "mv pht-openfoam-lab53iacp-master OpenFOAM-PHT53IACP" in the prompt.
3. Run "echo 'source $HOME/OpenFOAM/OpenFOAM-PHT53IACP/etc/bashrc' >> ~/.bashrc" in the prompt.
4. Run ". $HOME/.bashrc" in the prompt.
5. Change directory to OpenFOAM-PHT53IACP and run ./Allwmake script.

# Usage
The project **PHT OpenFOAM Lab53 IACP FEB RAS** is released under [Creative Commons Attribution-NonCommercial 4.0 International Public License](https://creativecommons.org/licenses/by-nc/4.0/legalcode) (see file LICENSE in the main directory or www.creativecommons.org/licenses/by-nc/4.0/legalcode). The project **PHT OpenFOAM Lab53 IACP FEB RAS** is released in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. Please feel free to extend, modify and use our solvers and libraries according to your requirement, but only for non-commerical usage. When using this project, please cite the project as follows:

For solver **phtEnthalpyFoam**, cite work [1].

For solver **porousHeatTransferFoam**, cite work [2].

# References

1. N. A. Lutsenko, S. S. Fetsov, K. G. Borovik, A. S. Kim. Gas flow and heat transfer in granular energy-releasing materials:
Novel computational model and important features // Int J Heat Mass Transf. 2022. V. 199, Article 123464. https://doi.org/10.1016/j.ijheatmasstransfer.2022.123464

2. S. S. Fetsov, N. A. Lutsenko. A novel computational model and OpenFOAM solver for simulating thermal energy storages based on granular phase change materials: Advantages and applicability // J Energy Storage. 2023. V. 65. Article 107294. https://doi.org/10.1016/j.est.2023.107294
