# PizDaint-mc cooldown and running instructions

**Required modifications** in <tt>cooldown.py</tt>:

	batch_tmpl = "job_daint_mc.slm"  # Add path if necessary.
	run_command = "srun -n $SLURM_NTASKS --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK"

Adjust **resources** (number of nodes and walltime) in <tt>cooldown.py</tt> if necessary:

	nodes = ...
	walltime = ...

On Daint-mc we recommend running **2 MPI tasks per node**, setting **OMP_NUM_THREADS = 18** and using **34 threads** in the Monte Carlo solver to leave space for the OS. Based on these recommendations we put these default values in <tt>job_daint_mc.slm</tt>:

	ntasks-per-node=2
	cpus-per-task=36
	
	export OMP_NUM_THREADS=18

In the input files then set *walkers* and *accumulators*  such that ***walkers* + *accumulators* = 34**, e.g.

	"walkers": 14
	"accumulators": 20