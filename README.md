# Aplysia Feeding Neuromechanical Model (Bennington et al. 2024)
Paper: "Incorporating buccal mass planar mechanics and anatomical features improves neuromechanical modeling of Aplysia feeding behavior"

System-level neuromechanical model of Aplysia feeding coupling a discrete-state Boolean neural model with a continuous rigid body biomechanical model.

Archived source code is available through Zenodo: [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.13773332.svg)](http://dx.doi.org/10.5281/zenodo.13773332)

Authors:

* Michael J. Bennington, CMU
* Ashlee S. Liao, CMU
* Ravesh Sukhnandan, CMU
* Bidisha Kundu, ULincoln, Bristol Medical School
* Stephen M. Rogers, ULincoln
* Jeffrey P. Gill, CWRU
* Jeffrey M. McManus, CWRU
* Gregory P. Sutton, ULincoln
* Hillel J. Chiel, CWRU
* Victoria A. Webster-Wood, CMU
Last Updated 9/17/2024

This repository contains all code and Simulink files required to implement the neuromechanical model presented in the above-referenced paper along with all of the computational experiments whose results are reported in the manuscript. Additionally, all animal datasets (mostly in the form of curves digitized from the literature) utilized in this paper are incorporated. All code is written for Matlab r2022b. Due to storage limitations of GitHub, the output data files for the computational experiments are not included in this repository but are archived in the Zenodo repository (link above). To recreate all data and figures presented in the paper, run the "ModelController.m" script. This script utilizes many of the other files in the repo, so all files should be kept in the same folder.

To perform comparison to the neuromechanical model presented in Webster-Wood et al. 2020 [1], the source code for that model in included here as well (see folder: WebsterWood2020_Model). The archived version of this code is available through Zenodo (doi:10.5281/zenodo.3978414).

[1] Webster-Wood, et al. "Control for multifunctionality : bioinspired control based on feeding in Aplysia californica." Biol. Cybern., 114(6):557-588, 2020.
