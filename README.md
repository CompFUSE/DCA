# DCA++

With the DCA++ project, we aim to gain insight into the fascinating physics of strongly correlated electron systems by employing modern quantum cluster methods.
The DCA++ code provides a state of the art implementation of the dynamical cluster approximation (DCA) and its DCA<sup>+</sup> extension.
High scalability and portable performance allow to exploit today's leadership computing systems.


## Getting Started

Information about prerequisites and the CMake-based building procedure can be found in the [Wiki](https://github.com/CompFUSE/DCA/wiki).


## Contributing Back

Questions, bug reports and feature requests can be submitted via the
[issue tracking system](https://github.com/CompFUSE/DCA/issues).

Please read the [contribution guidelines](https://github.com/CompFUSE/DCA/blob/master/CONTRIBUTING.md) for details on how to contribute code.


## Citation Guidelines

Please follow the [citation guidelines](https://github.com/CompFUSE/DCA/blob/master/CITATION.md), if you use the DCA++ code for scientific publications.


## CDash Dashboard

[Dashboard](http://cdash.cscs.ch/index.php?project=DCA) displaying the latest build status

## SC19 reproducibility

Please find in the folders `SC_paper_summit_runs` and `SC_paper_titan_runs` the inputs, scripts, and
executables used to produce the data reported in the Super Computing 2019 submission.
The runs for the new code were made using commit `7652c6956eaaa36ac375522b9571058c470ce47a`, while
commit `d25b256db2d61090963cd7b2975158557fd5f098` was used for the old version. Mind that only the
`main_dca` executable is guaranteed to work with the old code.
