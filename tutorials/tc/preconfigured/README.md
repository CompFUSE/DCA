This directory contains ready-to-use input files and jobs scripts for the [T<sub>c</sub> tutorial](https://github.com/CompFUSE/DCA/wiki/Tutorial:-Tc).
The files are preconfigured to run the example with a single process (no MPI) using the threaded Monte Carlo solver with 3 walkers and 5 accumulators.
The total number of measurements is 100 000.

The temperature steps are given by

    T = [1, 0.75, 0.5, 0.25, 0.125, 0.1, 0.09, 0.08, 0.07]

  and `T=0.1` is the first temperature for which we compute the leading d-wave eigenvalue.

Prerequisites
-------------
Configure and build the applications `main_dca` and `main_analysis` using the following CMake
options:

    CMAKE_BUILD_TYPE = Release
    DCA_BUILD_ANALYSIS = ON
    DCA_BUILD_DCA = ON
    DCA_CLUSTER_SOLVER = CT-AUX
    DCA_LATTICE = square
    DCA_MODEL = tight-binding
    DCA_POINT_GROUP = D4
    DCA_RNG = std::mt19937_64  # or std::ranlux48
    DCA_WITH_THREADED_SOLVER = ON

Usage
-----
We assume the *current working directory* to be the directory containing this file, i.e. `DCA/tutorials/tc/preconfigured`.

1. Copy the binaries `main_dca` and `main_analysis` into the current working directory.

2. Execute the [`dca` job script](https://github.com/CompFUSE/DCA/blob/master/tutorials/tc/preconfigured/job.dca_U=6_d=0.95_Nc=4.sh):

        $ ./job.dca_U=6_d=0.95_Nc=4.sh > out.dca.txt

    Runtime on an Intel Core i7-4980HQ @ 2.8GHz (4 cores) is ~110 min.  
    The file [`out.dca_reference.txt`](https://github.com/CompFUSE/DCA/blob/master/tutorials/tc/preconfigured/out.dca_reference.txt) provides reference output.
    Note that most of the numbers are subject to the stochastic Monte Carlo error.

3. When all `main_dca` runs are finished, execute the [`analysis` job script](https://github.com/CompFUSE/DCA/blob/master/tutorials/tc/preconfigured/job.analysis_U=6_d=0.95_Nc=4.sh):

        $ ./job.analysis_U=6_d=0.95_Nc=4.sh > out.analysis.txt

    The file [`out.analysis_reference.txt`](https://github.com/CompFUSE/DCA/blob/master/tutorials/tc/preconfigured/out.analysis_reference.txt) provides reference output.
    Most of the numbers again are subject to the stochastic Monte Carlo error of the `main_dca` runs.

4. When the `main_analysis` runs are finished, you can determine the superconducting transition temperature _T<sub>c</sub>_ with the Python script [`compute_tc.py`](https://github.com/CompFUSE/DCA/blob/master/tutorials/tc/compute_tc.py) provided in the parent directory:

        $ python ../compute_tc.py T=*

    The last command produces the plot `eigval_vs_temp.pdf` in the current working directory.
    You can compare the results with the reference plot [`eigval_vs_temp_reference.png`](https://github.com/CompFUSE/DCA/blob/master/tutorials/tc/eigval_vs_temp_reference.png) in the parent directory.
