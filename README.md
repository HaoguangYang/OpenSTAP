# OpenSTAP

### OpenSTAP is a Finite Element Method (FEM) Solver.

##### It is based on:
  - STAP (KJ Bath, FORTRAN IV)
  - STAP90 (Prof. Xiong Zhang, Tsinghua University, FORTRAN 90)

##### It includes:
  - Basic element types:
    - Truss
    - 3T (Triangular)
    - 4Q
    - 8Q
    - 9Q
    - 4T (Tetahedral)
    - 8H (Hexahedral)
    - 20H (B Matrix Unfinished)
    - Beam
    - Timoshenko Beam
    - Plate
    - Shell
    - Plate 8Q
    - Shell 8Q

  - Advanced element types:
    - Infinite 4Q
    - Plastic truss
 
  - Multiple solver choices
    - Classical Skyline storage and LDLT decomposition solver
    - Intel MKL Pardiso sparse matrix solver
    - LANCZOS eigen value solver (still not completely functional)
    - MKL FEAST eigen value solver (to be tested)
 
  - And am ABAQUS Python plugin for converting the input file of this program.

##### It works with:
 - Intel MKL
 - ParaVIEW for visualizing output data (actually it converts a text output file for printing and a vtk file in v3.0 format.)

##### It runs fast
Minimum RAM requirement: 400MB.
Benchmark of the program on a Intel i7-4790 Processor with 16G of RAM:
 - 4087 nodes with 2884 elements
     - about 3s using classical Skyline
     - about 0.2s using Intel MKL Pardiso
 - 37067 nodes with 30370 elements
     - about 3.5s using Intel MKL Pardiso
 - 259391 nodes with 232720 elements
     - about 35s using Intel MKL Pardiso

Our program also utilized MKL's out-of-core mode which allows external storage to be used when solving large problems. An example of which is running the third benchmark on a Surface Pro 4 with 8GB of RAM. Typically the program would create a ~950MB temporary file in the disk to store global stiffness matrix. This would also extend the solution time to about 2~3 minutes.

##### Plese refer to /doc for more information.
