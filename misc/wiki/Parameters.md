All applications read simulation parameters, e.g. the DCA cluster, the temperature, or output filenames, from a JSON-formatted input file, in which parameters are thematically grouped and sub-grouped in JSON objects. This page provides descriptions of all input parameters including their type and default value using the format  

`"name":` type (default value)  
Description.  

In addition, we provide a complete (including mutually exclusive groups) [sample input file](https://github.com/CompFUSE/DCA/blob/master/tools/complete_input.json).

## Output parameters

Defined in <tt>output_parameters.hpp</tt>.  

**Group** `"output":`

`"directory":` string ("./")  
Directory to write the output to.

`"autoresume":` bool (false)
If true, looks for a file named `<filename-dca>.tmp` generated from an aborted run, and starts 
the DCA loop from the last completed iteration. If the read is successful the parameter 
`initial-self-energy` is ignored.

`"directory-config-read":` string ("")  
If not empty, the Monte Carlo configuration will be initialized with the configurations stored in this directory.

`"directory-config-write":` string ("")  
If not empty, after the last Monte Carlo iteration, the configurations are written in this directory.

`"filename-dca":` string ("dca.hdf5")  
Filename for the output of the application <tt>main_dca</tt>.

`"filename-analysis":` string ("analysis.hdf5")  
Filename for the output of the application <tt>main_analysis</tt>.

`"filename-ed":` string ("ed.hdf5")  
Filename for the ED output in the application <tt>cluster_solver_check</tt>.

`"filename-qmc":` string ("qmc.hdf5")  
Filename for the QMC output in the application <tt>cluster_solver_check</tt>.

`"filename-profiling":` string ("profiling.json")  
Filename for the profiling output. The file format is always JSON.

`"dump-lattice-self-energy":` boolean (false)  
Write out the lattice self-energy in DCA<sup>+</sup>.

`"dump-cluster-Greens-functions":` boolean (false)  
Write out the cluster Green's functions.

`"dump-Gamma-lattice":` boolean (false)  
Write out the &Gamma; function of the BSE lattice solver.

`"dump-chi-0-lattice":` boolean (false)  
Write out the &chi;<sub>0</sub> function of the BSE lattice solver.

#### Example

    {
        "output": {
            "directory": "./T=0.5",
            "autoresume" : true,
            "filename-dca": "dca.hdf5",
            "filename-analysis": "analysis.hdf5",
            "filename-ed": "ed.hdf5",
            "filename-qmc": "qmc.hdf5",
            "filename-profiling": "profiling.json",
            "dump-lattice-self-energy": false,
            "dump-cluster-Greens-functions": true,
            "dump-Gamma-lattice": false,
            "dump-chi-0-lattice": false
        }
    }


## Physics parameters

Defined in <tt>physics_parameters.hpp</tt>.  

**Group** `"physics":`

`"beta":` double (1.)  
Inverse temperature.

`"density":` double (1.)  
Target electron density.

`"chemical-potential":` double (0.)  
Initial value of the chemical potential.  
Note that this value is overwritten, if an output file with an initial self-energy is given (see section DCA parameters).

`"adjust-chemical-potential":` boolean (true)  
Adjust chemical potential to obtain specified density.

#### Example

    {
        "physics": {
            "beta": 0.5,
            "density": 0.9,
            "chemical-potential": 0.,
            "adjust-chemical-potential": true
        }
    }


## Model parameters

Defined in <tt>model_parameters.hpp</tt>.  
The following groups are *mutually exclusive*.

### Single-band Hubbard model

**Group** `"single-band-Hubbard-model":`  
Used if *square lattice* or *triangular lattice* are selected.

`"t":` double (0.)  
Nearest neighbor hopping parameter.

`"t-prime":` double (0.)  
Next nearest neighbor hopping parameter.

`"U":` double (0.)  
On-site Coulomb repulsion.

`"V":` double (0.)  
Nearest neighbor coulomb repulsion, opposite spins.

`"V-prime":` double (0.)  
Nearest neighbor coulomb repulsion, same spins.

#### Example

    {
        "single-band-Hubbard-model": {
            "t": 1.,
            "t-prime": 0.5,
            "U": 8.,
            "V": 2.,
            "V-prime": 2.
        }
    }


### Bilayer Hubbard model

**Group** `"bilayer-Hubbard-model":`  
Used if *bilayer lattice* is selected.

`"t":` double (0.)  
Nearest neighbor hopping parameter.

`"t-prime":` double (0.)  
Next nearest neighbor hopping parameter.
 
`"t-perp":` double (0.)  
Hopping parameter between layers.

`"U":` double (0.)  
On-site Coulomb repulsion.

`"V":` double (0.)  
Coulomb repulsion between layers, opposite spins.

`"V-prime":` double (0.)  
Coulomb repulsion between layers, same spins.


#### Example

    {
        "bilayer-Hubbard-model": {
            "t": 1.,
            "t-prime": 0.5,
            "t-perp": 0.2,
            "U": 8.,
            "V": 2.,
            "V-prime": 2.
        }
    }


### Material model

**Group** `"material-model":`  
Used if a material lattice is selected.

`"t_ij-filename":` string ("t_ij.txt")  
Name of the CSV file containing hopping parameters *t<sub>ij</sub>*.

`"U_ij-filename":` string ("U_ij.txt")  
Name of the CSV file containing values of the Coulomb repulsion *U<sub>ij</sub>*.


#### Example

    {
        "material-model": {
            "t_ij-filename": "NiO_t_ij.txt",
            "U_ij-filename": "NiO_U_ij.txt"
        }
    }


## DCA parameters

Defined in <tt>dca_parameters.hpp</tt>.  

**Group** `"DCA":`

`"initial-self-energy":` string ("zero")  
Either the name of the file with the initial self-energy (usually the <tt>main_dca</tt> output file of the previous temperature) or "zero" indicating that the initial self-energy should be zero.

`"iterations":` integer (1)  
Number of DCA<sup>(+)</sup> iterations.

`"accuracy"`: double (0.)  
Stop the DCA<sup>(+)</sup> loop if this accuracy has been reached.

`"self-energy-mixing-factor":` double (1.)  
&Sigma;<sup>(n)</sup> = &alpha; &Sigma;<sup>(n)</sup><sub>QMC</sub> + (1-&alpha;)&Sigma;<sup>(n-1)</sup><sub>coarsegrained</sub>, where &alpha; is the self-energy mixing factor.  

`"interacting-orbitals":` array of integers ([0])  
Indices of orbitals that are treated interacting.  
Note that this parameter must be consistent with the model that is used.

`"do-finite-size-QMC":` boolean (false)  
Do a finite-size QMC calculation (no mean-field).


<br></br>
**Subgroup** `"coarse-graining":`  
`"k-mesh-recursion":` integer (0)  
Number of recursion steps in the creation of the momentum space mesh.  
See the appendix of [1] for details.

`"periods":` integer (0)  
Number of "periods" of the interlaced coarse-graining patches, i.e. degree of their interleaving.    
Restricted by *k-mesh-recursion*. See the appendix of [1] for details.

`"quadrature-rule":`  integer (1)  
Determines the quadrature rule used in the coarse-graining.  
quadrature-rule < 0: use a flat mesh.  
quadrature-rule >= 0: use the Grundmann-Moeller rule of index = *quadrature-rule*. A rule of index *s* is exact for polynomials up to order 2*s*+1 [2].

`"threads":` integer (1)  
Number of threads used in the coarse-graining.

`"tail-frequencies":` integer (0)  
Number of tail frequencies used for updating the chemical potential.


<br></br>
**Subgroup** `"DCA+":`  
`"do-DCA+"`: boolean (false)  
Use the DCA<sup>+</sup> algorithm with a continuous lattice self-energy.

`"deconvolution-iterations":` integer (16)  
Maximum number of iterations in the deconvolution step.

`"deconvolution-tolerance":` double (1.e-3)  
Termination criteria for the deconvolution step.

`"HTS-approximation":` boolean (false)  
Use high-temperature series approximation for lattice mapping.

`"HTS-threads": ` integer (1)  
Number of threads used in the HTS-solver.

#### Example

    {
        "DCA": {
            "initial-self-energy": "./T=0.5/dca.hdf5",
            "iterations": 3,
            "accuracy": 1.e-3,
            "self-energy-mixing-factor": 0.5,
            "interacting-orbitals": [0],

            "do-finite-size-QMC": false,

            "coarse-graining": {
                "k-mesh-recursion": 3,
                "periods": 2,
                "quadrature-rule": 1,
                "threads": 8,
                "tail-frequencies": 10
            },

            "DCA+": {
                "do-DCA+": true,
                "deconvolution-iterations": 16,
                "deconvolution-tolerance": 1.e-3,
                "HTS-approximation": true,
                "HTS-threads": 8
            }
        }
    }
	

## Domains parameters

Defined in <tt>domains_parameters.hpp</tt>.

**Group** `"domains":`

<br></br>
**Subgroup** `"real-space-grids":`  
`"cluster":` array of arrays of integers (lattice basis, e.g. in 2D: [ [1, 0], [0, 1] ])  
Real space DCA cluster.  
Given as coordinates with respect to the lattice basis.  
To choose a cluster of a given size, consult this [exhaustive list](https://github.com/CompFUSE/DCA/blob/master/tools/cluster_definitions.txt) of 2D and 3D clusters.

`"sp-host":` array of arrays of integers (lattice basis, e.g. in 2D: [ [1, 0], [0, 1] ])  
Real space host grid for single-particle functions, also called *(sp-)lattice*.  
In addition, it is used in the coarse-graining step of **both DCA and DCA<sup>+</sup>** as the grid for the free dispersion relation &epsilon;<sub>**k**</sub>.  
Given as coordinates with respect to the lattice basis.

`"tp-host":` array of arrays of integers (lattice basis, e.g. in 2D: [ [1, 0], [0, 1] ])  
Real space host grid for two-particle functions, also called *tp-lattice*. Only used in **DCA<sup>+</sup>**.  
Given as coordinates with respect to the lattice basis.


<br></br>
**Subgroup** `"imaginary-time":`  
`"sp-time-intervals":` integer (128)  
Used to initialize <tt>time_domain</tt>, whose elements are  
&nbsp;&nbsp;&nbsp;&nbsp; -&beta;+&epsilon;, -&beta;+&Delta;, ..., -&Delta;, -&epsilon;, &epsilon;, &Delta;, ..., &beta;-&Delta;, &beta;-&epsilon;,  
where &Delta; = &beta;/*sp-time-intervals*.  
More precisely, it determines the *G*<sub>0</sub>-interpolation domain in CT-AUX and the *F*-interpolation domain in SS-CT-HYB.  
**TODO**: Make this clearer.

`"time-intervals-for-time-measurements":` integer (1)  
Used to initialize <tt>vertex_time_domain<TP_TIME_DOMAIN></tt> and <tt>vertex_time_domain<TP_TIME_DOMAIN_POSITIVE></tt>, where the latter is the domain for additional time measurements.


<br></br>
**Subgroup** `"imaginary-frequency"`  
`"sp-fermionic-frequencies":` integer (256)  
Number of positive and negative fermionic Matsubara frequencies, respectively.  
Used to initialize <tt>frequency_domain</tt>, whose elements are (*n* = sp-fermionic-frequencies)  
&nbsp;&nbsp;&nbsp;&nbsp; -(2*n*-1)&pi;/&beta;, ..., -3&pi;/&beta;, -&pi;/&beta;, &pi;/&beta;, 3&pi;/&beta;, ..., (2*n*-1)&pi;/&beta;.

`"HTS-bosonic-frequencies":` integer (0)  
Used to initialize <tt>vertex_frequency_domain<EXTENDED_BOSONIC></tt>, which is the bosonic Matsubara frequency mesh employed in the HTS-solver.
The elements of the domain are (*n* = HTS-bosonic-frequencies):  
&nbsp;&nbsp;&nbsp;&nbsp; -2*n*&pi;/&beta;, ..., -2&pi;/&beta;, 0, &pi;/&beta;, ..., 2*n*&pi;/&beta;

`"four-point-fermionic-frequencies":`  integer (1)  
Used to initialize <tt>vertex_frequency_domain<COMPACT/COMPACT_POSITIVE/EXTENDED/EXTENDED_POSITIVE></tt>, that is the fermionic Matsubara frequency mesh for measuring four-point quantities.

<br></br>
**Subgroup** `"real-frequency":`  
`"min":` double (-10.)  
`"max":` double (10.)  
`"frequencies":` integer (3)  
Parameters for the initialization of <tt>frequency_domain_real_axis</tt>:  
&nbsp;&nbsp;&nbsp;&nbsp; *min*,  *min*+&Delta;, *min* + 2&Delta;, ..., *max*-&Delta;, *max*, &nbsp;&nbsp; &Delta; = (*max*-*min*)/(*frequencies*-1).

`"imaginary-damping":` double (0.01)  
Small *positive* shift, used to shift the real frequencies &omega; away from the real axis, &omega; -> &omega;+i&eta;.

#### Example

    {	
        "domains": {
            "real-space-grids": {
                "cluster": [[2, 0],
                            [0, 2]],
                "sp-host": [[20, 20],
                            [20,-20]],
                "tp-host": [[8, 8],
                            [8,-8]]
            },

            "imaginary-time": {
                "sp-time-intervals": 128,
                "time-intervals-for-time-measurements": 128
            },

            "imaginary-frequency": {
                "sp-fermionic-frequencies": 256,
                "HTS-bosonic-frequencies": 32,
                "four-point-fermionic-frequencies": 16
            },

            "real-frequency": {
                "min": -10,
                "max": 10,
                "frequencies": 128,
                "imaginary-damping": 0.01
            }
        }
    }


## Monte Carlo integration parameters

Defined in <tt>mci_parameters.hpp</tt>.

**Group** `"Monte-Carlo-integration":`

`"seed":` integer | string (985456376)  
Seed for the random number generator(s) used in the Monte Carlo integration.  
Instead of an integer, the string `"random"` can be passed to generate a random seed.

`"warm-up-sweeps":` integer (20)  
Number of warm-up sweeps.

`"sweeps-per-measurement":` double (1.)  
Number of sweeps per measurement.

`"measurements":` integer (100)  
Number of independent measurements in each Monte Carlo iteration.

`"error-computation-type:"`  string("NONE")\
Determines the type of error computation that will be performed during the last Monte Carlo iteration. Possible options are:
- "NONE"
- "STANDARD_DEVIATION"
- "JACK_KNIFE"

`"time-correlation-window":` integer (0)  
 Maximum distance (in MC time) considered when computing the correlation between configurations.
 If 0, no auto-correlation is computed.
 
`"compute-G-correlation":` boolean (true)
If `time-correlation-window` is larger than 0, G(r = 0, t = 0) is included in the observables
whose autocorrelation is computed. This measurements requires some device memory. 
 
`"stamping-period"` integer (0)  
If larger than 0, the master MPI rank logs the walker configuration every `stamping-period` sweeps. 
 
`"store-configuration":` : boolean (false)
If true, the vertex configuration is stored between DCA iterations to initialize the walkers of 
the following iteration.

<br></br>
**subgroup** `"threaded-solver":`  
Additional parameters for threaded Monte Carlo solvers.  
`"walkers":` integer (1)  
Number of Monte Carlo walkers.  
`"accumulators":` integer (1)  
Number of Monte Carlo accumulators.  
`"shared-walk-and-accumulation-thread":` boolean (false)\
When this mode is activated, each thread will run an instance of a walker and an accumulator. This parameter is ignored if the numbers of walkers and accumulators are different.\
`"fix-meas-per-walker":` boolean(false)\
The number of sweeps performed by each walker is fixed a priori, avoiding possible bias in the sampling of faster walkers.

#### Example
	
    {
        "Monte-Carlo-integration": {
            "seed": 42,
            "warm-up-sweeps": 40,
            "sweeps-per-measurement": 4.,
            "measurements": 1000000,
            "error-computation-type": "JACK_KNIFE",
            "time-correlation-window" : 100,
            "compute-G-correlation" : true,
            "stamping-period" : 0,
            "store-configuration" : true,

            "threaded-solver": {
                "walkers": 3,
                "accumulators": 5,
                "shared-walk-and-accumulation-thread": false
            }
        }
    }
	
	
## Monte Carlo solver parameters

Defined in <tt>mc_solver_parameters.hpp</tt>.  
The following parameter groups are *mutually exclusive*.

### CT-AUX

**Group** `"CT-AUX":`  
Used if *CT-AUX* is selected as the cluster solver.

`"expansion-parameter-K":` double (1.)  
The perturbation order in the CT-AUX algorithm increases linearly with the expansion parameter *K*.  
While *K* is only subject to the restriction of being positive, values of the order of 1 have proven to be a good choice [3].

`"initial-configuration-size":` integer (10)  
The CT-AUX solver is initialized with `"initial-configuration-size"` random **interacting** vertices.

`"initial-matrix-size":` integer (128)  
Initial size of the CT-AUX matrices.

`"max-submatrix-size":` integer (128)  
Maximum number of single spin updates per submatrix update (see [4]).

`"neglect-Bennett-updates":` boolean (false)  
Neglect the two types of Bennett updates:
1. removal of an *interacting* spin that has already been proposed for removal,  
2. removal of a *non-interacting* spin that has been proposed for insertion.  

Turning on this option leads to a larger number of delayed spins (~`"max-submatrix-size"`) at the cost of a systematic error due to the violation of detailed balance.


`"additional-time-measurements":` boolean (false)  
Do additional time measurements.

#### Example

    {
        "CT-AUX": {
            "expansion-parameter-K": 1.,
            "initial-configuration-size": 100,
            "initial-matrix-size": 128,

            "max-submatrix-size": 64,
            "neglect-Bennett-updates": true,

            "additional-time-measurements": true
        }
    }

### CT-AUX

**Group** `"CT-INT":`  
Used if *CT-INT* is selected as the cluster solver.

`"initial-configuration-size"` integer(0)  
Size of the initial random configuration (if a previous configuration is not read).   
`"alpha-dd-pos"` double(0.501)  
Strength of the alpha field. A random variable in {`alpha-dd-pos`, -`alpha-dd-pos`} is added to the 
diagonal of the D matrix to mitigate the sign problem. This parameter applies only to 
density-density interactions with positive coupling.   
`"alpha-dd-neg"` double(0.)  
Strength of the alpha field for interactions with negative coupling.  
`"adjust-alpha-dd"` double(1e-4)  
If true, the value of G0(0+) is added to `alpha-dd-pos`.  
`adjust_alpha-dd"` boolean(false)  
Strength of the alpha field for non-density-density interactions.  
`"double-update-probability"` double(0.)  
Probability to propose a double vertex insertion/removal.  
`"all-sites-partnership"` bool(false)  
Determines if a double update can be proposed with vertices on different sites.    
`"max-submatrix-size"` integer(1)  
Maximum number of vertices inserted/removed in a single submatrix update.  

#### Example

    {
        "CT-INT": {
            "initial-configuration-size": 1000,
            "alpha-dd-pos": 0.1,
            "alpha-dd-neg": 0,
            "alpha-ndd" : 0.01,
            "adjust-alpha-dd" : true,
            "double-update-probability" : 0.5,
            "all-sites-partnership" : false,
            "max-submatrix-size": 128,
        }
    }

### CT-HYB

**Group** `"SS-CT-HYB":`  
Used if *SS-CT-HYB* is selected as the cluster solver.

`"self-energy-tail-cutoff":` integer (0)  
Cut-off parameter for imaginary frequency tail of self-energy.

`"steps-per-sweep":` double (0.5)  
`"shifts-per-sweep":` double (0.5)  
The fraction of insert/removal steps is given by  
&nbsp;&nbsp;&nbsp;&nbsp; *steps-per-sweep/(steps-per-sweep + shifts-per-sweep)*,  
while the fraction of segments shifts is  
&nbsp;&nbsp;&nbsp;&nbsp; *shifts-per-sweep/(steps-per-sweep + shifts-per-sweep)*.

#### Example

    {
        "SS-CT-HYB": {
            "self-energy-tail-cutoff": 10,
            "steps-per-sweep": 0.6,
            "shifts-per-sweep": 0.4
        }
    }


## Four-point parameters
 
Defined in <tt>four_point_parameters.hpp</tt>.

**Group** `"four-point":`

`"channels":` array of strings ([])  
List of four point channels to be accumulated in the main dca application. Select only one for the analysis application.  
The possible list elements are:
* "PARTICLE_PARTICLE_UP_DOWN"
* "PARTICLE_HOLE_CHARGE"
* "PARTICLE_HOLE_MAGNETIC"
* "PARTICLE_HOLE_LONGITUDINAL_UP_UP"
* "PARTICLE_HOLE_LONGITUDINAL_UP_DOWN"
* "PARTICLE_HOLE_TRANSVERSE"

`"type":` string ("NONE") Deprecated.  
Add a single channel to the accumulated channels list.

`"momentum-transfer":` array of doubles (null vector with the dimension of the lattice, e.g. for a 2D lattice: [0., 0.])  
Transferred momentum **q**.  
Must be an element of the reciprocal lattice of the DCA cluster (|| **q** - **K**<sub>DCA</sub> ||<sub>2</sub> < 10<sup>-3</sup>).  

`"frequency-transfer":` integer (0)  
Transferred frequency &omega;.   
Given as the index *n* of the bosonic Matsubara frequency &omega;<sub>*n*</sub> = 2*n*&pi;/&beta;.

`"compute-all-transfers":` boolean (false)\
When this mode is activated all possible positive frequency transfers, up to and including the transfer with index "frequency-transfer". All possible momentum transfers will be computed as well, ignoring the parameter "momentum-transfer".

#### Example
	
    {
        "four-point": {
            "channels": ["PARTICLE_PARTICLE_UP_DOWN"],
            "momentum-transfer": [0., 0.],
            "frequency-transfer": 0,
            "compute-all-transfers": false
        }
    }
		
		
## Analysis parameters

Defined in <tt>analysis_parameters.hpp</tt>.

**Group** `"analysis":`

`"symmetrize-Gamma":` boolean (true)  
Symmetrize &Gamma;<sub>cluster</sub> and &Gamma;<sub>lattice</sub> according to the symmetry group and diagrammatic symmetries.

`"Gamma-deconvolution-cut-off":` double (0.5)  
Cut-off parameter for the deconvolution of the vertex &Gamma;(**K**<sub>1</sub>, **K**<sub>2</sub>).  
The deconvolution is done by Fourier transforming to real space, dividing by the coarse-graining patch functions in real space &phi;(**r**) and then transforming back to reciprocal space. The real space patch functions are only inverted down to the value of this parameter,  
&nbsp;&nbsp;&nbsp;&nbsp; &phi;<sup>-1</sup><sub>cut-off</sub>(**r**) = &phi;<sup>-1</sup>(**r**),&nbsp;&nbsp; if &phi;(**r**) > *Gamma-deconvolution-cut-off*,  
&nbsp;&nbsp;&nbsp;&nbsp; &phi;<sup>-1</sup><sub>cut-off</sub>(**r**) = 0,&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; else.

`"project-onto-chrystal-harmonics":` boolean (false)  
Project &Gamma; &chi;<sub>0</sub> onto the crystal harmonic functions and diagonalize in this subspace.

`"projection-cut-off-radius":` double (1.5)  
For the projection use crystal-harmonic functions of lattice vectors **r** with ||**r**||<sub>2</sub> < *projection-cut-off-radius*.


#### Example

    {
        "analysis": {
            "symmetrize-Gamma": true,
            "Gamma-deconvolution-cut-off": 0.5,
            "project-onto-crystal-harmonics": true,
            "projection-cut-off-radius": 1.5
        }
    }


## ED-solver parameters
Defined in <tt>ed_solver_parameters.hpp</tt>.  
Only required if the exact diagonalization solver is used, e.g. in the application <tt>cluster_solver_check</tt>.

**Group** `"ED":`

`"eigenvalue-cut-off":` double (1.e-6)  
Only keep energy eigenvalues *E* for which e<sup>-*&beta;E*</sup> > *eigenvalue-cut-off*.

`"threads":` int (1)  
Number of threads participating in the ED solver.


#### Example

    {
        "ED": {
            "eigenvalue-cut-off": 1.e-6,
            "threads" : 8
        }
    }


## Double-counting parameters

Defined in <tt>double_counting_parameters.hpp</tt>.  
Only required for LDA+DMFT.

**Group** `"double-couting":`

`"method":` string ("none")  
The double-counting method to use, options are:

* "none": no double-counting correction
* "constant-correction-without-U-correction"
* "constant-correction-with-U-correction"

`"correction"`: double (0.)  
The value of the double-counting correction.

#### Example

    {
        "double-counting": {
            "method" : "constant-correction-without-U-correction",
            "correction" : 4.2
        }
    }

			
### References
[1]	P. Staar, M. Jiang, U. R. Hähner, T. C. Schulthess, and T. A. Maier, Physical Review B 93, 165144 (2016).  
[2] https://people.sc.fsu.edu/~jburkardt/cpp_src/simplex_gm_rule/simplex_gm_rule.html  
[3] E. Gull, P. Werner, O. Parcollet, and M. Troyer, Europhys. Lett. 82, 57003 (2008).  
[4] E. Gull, P. Staar, S. Fuchs, P. Nukala, M. S. Summers, T. Pruschke, T. C. Schulthess, and T. Maier, Physical Review B 83, 075122 (2011).
