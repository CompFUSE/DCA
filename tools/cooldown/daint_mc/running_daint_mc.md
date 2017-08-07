# Piz Daint (multi-core) cooldown and running instructions

**Required modifications** in <tt>cooldown.py</tt>:

	batch_tmpl = "job_daint_mc.slm"  # Add path if necessary.
	run_command = "srun -n $SLURM_NTASKS --ntasks-per-core=$SLURM_NTASKS_PER_CORE --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK"

Specify **resources** (number of nodes and walltime) in <tt>cooldown.py</tt>:

	nodes = ...
	walltime = ...

On Piz Daint (multi-core) we recommend using **2 MPI tasks per node**, setting **OMP\_NUM\_THREADS = 18**, and running **34 threads** in the Monte Carlo solver to leave space for the OS.  Accordingly, `"walkers"` and `"accumulators"` need to be set in `input_sp.json.in` and `input_tp.json.in` such that **walkers + accumulators = 34**, e.g.

	"walkers": 14
	"accumulators": 20