# PizDaint-gpu cooldown and running instructions

**Required modifications** in <tt>cooldown.py</tt>:

	batch_tmpl = "job_daint_gpu.slm"  # Add path if necessary.
	run_command = "srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK"

Adjust **resources** (number of nodes and walltime) in <tt>cooldown.py</tt> if necessary:

	nodes = ...
	walltime = ...

On Daint-gpu we recommend running **1 MPI task per node**, setting **OMP_NUM_THREADS = 12** and using **22 threads** in the Monte Carlo solver to leave space for the OS. Based on these recommendations we put these default values in <tt>job_daint_gpu.slm</tt>:

	ntasks-per-node=1
	cpus-per-task=24
	
	export OMP_NUM_THREADS=12

In the input files then set *walkers* and *accumulators*  such that ***walkers* + *accumulators* = 22**, e.g.

	"walkers": 8
	"accumulators": 14