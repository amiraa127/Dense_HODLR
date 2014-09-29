#RIGID BODY
This folder contains examples on how to use the HODLR solver to effciciently solve the particale in a stokes flow problem. Currently, `lattice.cpp` contains a sample code corresponding to particles in a structured grid. 

####Build
Follow the instruction on building the package as explained in the root directory. After running CMake, simply go to the `rigidBody` directory under `build` and type `make`.

####Usage 
```
./latticeTest path/to/vertex/file d N_x N_y N_z
```
where `d` is the lattice parameter and `N_x`, `N_y` and `N_z` are the grid dimensions in the x, y and z coordinates respectively. 

####Note
You might notice that the solver takes some time. This is because we calculate both the actual error and the relative l2 error. Hence, we take an arbitary `x_e` to be the solution. We then compute the right hand side `F` using a dumb (N^2) matrix vector multiplication module:
```
F = A * x_e
```
We then use the HODLR solver to solve for `x_HODLR`:

```
A * x_HODLR = F
```

The l2 and actual error are calculates as follows:

```
e_l2 = ||A * x_HODLR - F|| / ||F||
e_a  = ||x_HODLR - x_e|| / ||x_e||
```
The dumb matrix vector product used to calculate both `F` and `e_l2` is the most time consuming part of the calculation.