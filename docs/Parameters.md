# Parameters

## Monte Carlo Integration (MCI) parameters

Monte Carlo Integration parameters are defined in the group `"Monte-Carlo-Integration"`.
	
### Seed

Seed for the random number generator used in the Monte Carlo integration.

Field name: `"RNG-seed"`    
Value type: `int | "random"`  
Default value: `985456376`

If the seed option `"random"` is given, a random seed is generated.

### Example
	
	{
    	"Monte-Carlo-Integration": {
    		"Sigma-file": "sigma.hdf5",

			"warm-up-sweeps": 20,
			"sweeps-per-measurement": 4.,
			"measurements": 100,

			"adaptive-double-counting": "false",

			"RNG-seed": 42,
			
			"MC-posix-parameters": {
				"nr-walkers": 3,
				"nr-accumulators": 5,
				"additional-steps": 1,
				
				"HTS-threads": 1
			}
		}
	}