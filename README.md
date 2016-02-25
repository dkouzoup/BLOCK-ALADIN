# BLOCK-ALADIN
A prototype implementation of the block-based ALADIN scheme for optimal control

# INSTRUCTIONS

- Open script runOCP.m, insert the simulation details of your choice and click run.

- The integrator included in this repository is compiled for MAC OS X. To use on other operating systems, the user must install the latest version of ACADO (master branch is recommended) and set the EXPORT_INTEGRATOR flag to 1.

- To generate MEX interfaces for the algorithmic parts that are written in C, set the GENERATE_MEX flag to 1 (provided files are for MAC OS X and block sizes 5 and 10). 

- Users must also install the dense QP solver qpOASES (http://www.qpoases.org/)