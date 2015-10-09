# BLOCK-ALADIN
A prototype implementation of the block-based ALADIN scheme for optimal control

# INSTRUCTIONS

- Simply open script runOCP.m, insert the simulation details of your choice and click run.

- The integrator included in this repository is compiled for MAC OS X. To use on other operating systems, the user must install the latest version of ACADO (master branch) and run the generateSim.m script.

- For non MAC OS X users, the provided QP solver qpDUNES must be also compiled. Simply type "make" on the MATLAB command window inside the folder qpdunes_dev/interfaces/matlab.

- The algorithm at this stage is provided as a proof of concept for the block based ALADIN scheme. The complete code containing all reported implementations and timings will become public once the paper is accepted for publication.