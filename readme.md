## About GTasb3D

GTasb3D is a 3D thermophysical code using the General finite difference model to simulate the Thermal evolution of active small bodies. The features of this code includes:
* Capable of simulating various thermophysical processes involved in the transfer of heat fluxes, the sublimation/condensation of ice, and the diffusion of gas through a porous matrix within an icy body.
* Applicable to modeling a body with various kinds of shapes.
* Parallelized implimentaion in C++ with reasonable accuracy and efficiency.
More details on the mathematical formulate, code architecture, and validation tests can be found in <a href="https://doi.org/10.3847/PSJ/acc4c4">Zhang \& Hartzell (2023)</a>

## About this repository

Here includes the GTasb3D source code as well as an example with input files and MATLAB scripts. 

To compile the code, create a folder, e.g., "bin", in the top directory, type
```sh
cmake ../
```
This will generate a Makefile in the folder, and type
```sh
make
```
Voilà!  You should find the executable file named "GTasb3D" in this folder, and now you can start to perform a GTasb3D simulation.

Steps to run GTasb3D:
* Copy the obtained GTasb3D to the simulation directory, e.g., example/133P.
* Edit "gfdm.dat" file to adjust the simulation setup and physical parameters.
* Create a folder named "result" in this directory and launch a simulation by 
```sh
./GTasb3D -n $NUMBER_OF_PROCESSORS gfdm.dat
```

Please cite <a href="https://doi.org/10.3847/PSJ/acc4c4">Zhang \& Hartzell (2023)</a> if you use this code for your study. Please contact the author via the gmail address at yzhangastro if you would like to have more details about how to use and modify this code. 
