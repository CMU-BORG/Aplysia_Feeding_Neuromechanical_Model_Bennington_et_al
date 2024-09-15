# Aplysia Feeding Neuromechanical Model (Bennington et al. 2024)
Paper: "Incorporating buccal mass planar mechanics and anatomical features improves neuromechanical modeling of _Aplysia_ feeding behavior"

System-level neuromechanical model of Aplysia feeding coupling a discrete-state Boolean neural model with a continuous rigid body biomechanical model. 

Archived source code is available through Zenodo: 

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

Last Updated 9/13/2024

This repository contains all code and Simulink files required to implement the neuromechanical model presented in the above-referenced paper along with all of the computational experiments whose results are reported in the manuscript. Additionally, all animal datasets (mostly in the form of curves digitized from the literature) utilized in this paper are incorporated. All code is written for Matlab r2022b. To recreate all figures presented in the paper, the "ModelController.m" script can be run. This script utilizes many of the other files in the repo, so all files should be kept in the same folder.
