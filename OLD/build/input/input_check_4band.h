{
    "filename-parameters" : 
    {
	"output-format"         : "HDF5",

	"output-ED"  : "output__ED.hdf5",
	"output-CPE" : "output_CPE.hdf5",
	"output-QMC" : "output_QMC.hdf5"	
    },

    "physics-parameters" :
    { 
        "beta"                      :  10,

	"adjust-chemical-potential" : "false", 
        "density"                   :  1.0,
        "chemical-potential"        :  0.2
    },

    "fourband-model" : 
    {
	"ei0" : 0.6,
	"eb0" : 0.1,
	"t0"  : 0.5,
	"ei1" : 0.1,
	"eb1" : 0.6,
	"t1"  : 1,
        "U0"  : 7, 
	"U1"  : 8,
	"V"  :  6,
        "V_prime"  :  5
       },

    "DCA" : 
    {
    	"do-DCA+" : "false",
	"interacting-bands" : [0,1],
	"DCA-iterations"    : 2,
	"DCA-mixing-factor" : 1.0,
	
	"cluster" : [[ 1, 0],
		     [ 0, 1]],

	"cluster-mapping" : 
	{
	   "k-mesh-refinement" : 3,
	   "quadrature-rule"   : 1,
	   "number-of-periods" : 0,

	   "precompute-Hamiltonian"      : "true",
	   "phi(k) integration accuracy" : 0.01,

	   "print-phi"                   : "false"
	},

	"lattice-mapping" : 
	{
	   "interpolation-method"         : "wannier-interpolation",

	   "deconvolution-tolerance"      : 0.01,
	   "max-deconvolution-iterations" : 16
	}
    },

    "Monte-Carlo-Integration" : 
    {
	"Sigma-file" : "zero",

	"warm-up-sweeps"         : 200,
	"sweeps-per-measurement" : 10,
	"measurements"           : 100000,

	"adaptive-double-counting" : "false",
	
	"RNG-seed" : 985456376,

	"MC-posix-parameters" : 
	{
	    "nr-walkers"       : 2,
	    "nr-accumulators"  : 2,
	    "additional-steps" : 0
	}
    },

    "CT-AUX-solver" : 
    {
	"submatrix-size"      : 16,
	"initial-matrix-size" : 256,

	"K-parameter"         : 100
    },

    "SS-CT-HYB-solver" : 
    {
	"steps-per-sweep"  : 1,
	"swaps-per-sweep"  : 0,
	"shifts-per-sweep" : 1,

	"Sigma-tail-cutoff" : 50
    },

    "ED-solver-parameters" : 
    {   
    	"eigenvalue-cut-off" : 1e-8
    },

    "function-parameters" : 
    {
	"single-particle-functions" : 
	{
	    "H(k) grid-size"        : [1, 1],	

	    "time-intervals"        : 128,

	    "fermionic-frequencies" : 512,
	    "bosonic-frequencies"   : 128,
	    
	    "sp-cluster" : [[8, 8],
		      	    [8,-8]]
	},

	"two-particle-functions" : 
	{
	    "time-intervals"        : 32,

	    "fermionic-frequencies" : 16,
	    "bosonic-frequencies"   : 16,
	
	    "tp-cluster" : [[8, 8],
		      	    [8,-8]]
	},

	"real-axis-functions" :
	{
	    "lower-bound" : -10,
	    "upper-bound" :  10,

	    "nr-intervals" : 256,
	    "real-axis-off-set" : 0.01
	}
     },

    "equal-time-observables" : 
    {
        "do-equal-time-measurements" : "false"
    },

    "CPE-parameters" : 
    {
	"do-CPE" : "false",

	"max-CPE-iterations" : 1000,
	"max-CPE-error"      : 0.01,

	"number-of-matsubara-freqencies" : 32,

	"smoothing-factor" : 1,	
	
	"simulate-Gaussian-noise" : "true",
	"nr-of-samples"           : 100,
	"simulated-stddev"        : 0.01,

	"compute-free-spectrum"    : "false",
	"compute-lattice-spectrum" : "false",
	"compute-cluster-spectrum" : "true"
    }
}
