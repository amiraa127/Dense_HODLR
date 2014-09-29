#RIGID BODY
This folder contains examples on how to use the HODLR solver to effciciently solve the particale in a stokes flow problem. Currently, `lattice.cpp` contains a sample code corresponding to particles in a structured grid. 

####Build
To build the code simply run `make` in the `rigidBody` directory.

####Usage 
```
./latticeTest path/to/vertex/file d N_x N_y N_z
```
where `d` is the lattice parameter and `N_x`, `N_y` and `N_z` are the grid dimensions in the x, y and z coordinates respectively. 