#HODLR PACKAGE : Hierarchically off-diagonal low-rank matrix manipulation package 

This software package is a highly flexible HODLR solver that can be used in a variety of applications. The package can solve matrices that are preallocated in memory or defined by a kernel function. 

A variety of low-rank approximation methods (SVD, partial pivoting ACA, full pivoting ACA, BDLR, etc.) are available. Furthermore, the package also accepts any custom HODLR partitioning from the user which is very useful in many applications. 

####Author :  

Amirhossein Aminfar: amirhossein.aminfar@gmail.com

####Citation:

If you use the implementation or any part of the implementation in your work, kindly cite as follows:

####Articles:

@article{aminfar2014fast,

author={{A}minfar, {A}mirhossein and {A}mbikasaran, {S}ivaram and {D}arve, {E}ric},

title={A Fast Block Low-Rank Dense Solver with Applications to Finite-Element Matrices},

journal={arXiv:1403.5337},

year={2014}

}


@article{ambikasaran2014fastdet,

title={Fast Direct Methods for {G}aussian Processes and the Analysis of {NASA} {K}epler Mission Data},

author={Ambikasaran, Sivaram and Foreman-Mackey, Daniel and Greengard, Leslie F. and Hogg, David W. and O'Neil, Michael},

journal={arXiv preprint arXiv:1403.6015},

year={2014}

}

@article{ambikasaran2013fast,

title={An $\mathcal{O}({N} \log {N})$ Fast Direct Solver for Partial Hierarchically Semi-Separable Matrices},

author={Ambikasaran, Sivaram and Darve, Eric F.},

journal={Journal of Scientific Computing},

volume={57},

number={3},

pages={477--501},

year={2013},

publisher={Springer}

}

@article{ambikasaran2014fastsym,

title={Fast symmetric factorization of hierarchical matrices with applications},

author={Ambikasaran, Sivaram and O'Neil, Michael},

journal={arXiv preprint arXiv:1405.0223},

year={2014}

}

####Code

@MISC{aminfar2014HODLR,

author = {{A}minfar, {A}mirhossein},

title = {A fast direct solver for HODLR matrices},

howpublished = {https://github.com/amiraa127/Dense_HODLR},

year = {2014}

}

####Version 1.00

Date: September 30th, 2014

Copyleft 2014: Amirhossein Aminfar 

Developed by Amirhossein Aminfar

####License


This program is free software; you can redistribute it and/or modify it under the terms of MPL2 license. The Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

####Build
This package has the following dependencies:

1. [Eigen] (http://eigen.tuxfamily.org/index.php?title=Main_Page) C++ library for all matrix manipulations.

2. [CPP Unit] (http://sourceforge.net/projects/cppunit/) for unit testing.

3. (Optional) [SOTCH] (http://www.labri.fr/perso/pelegrin/scotch/) for nested dissection and matrix reordering.

4. (Optional) [PaStiX] (http://pastix.gforge.inria.fr/files/README-txt.html) for the extended sparsification method. It is important to compile PaStiX with the -DFORCE_NOMPI flag and link PaStiX to SCOTCH not ptSCOTCH.

The easiest way to build the library is to use [CMake](http://www.cmake.org).
Go to the project directory and run:

```
mkdir build
cd build
cmake ../
make
```

The new version of the HODLR package can take advantage of PaStiX for factorizing HODLR matrices using the extended sparsification method. This requires the installation of PaStiX and SCOTCH for the sparse embedding functionality. To compile with PaStiX support, use:

```
mkdir build
cd build
cmake -DPASTIX=ON ../
make
```

####Documentation

The documentation for this package is currently located at http://amiraa127.github.io/Dense_HODLR/. It's currently incomplete but I will add material to it gradually.