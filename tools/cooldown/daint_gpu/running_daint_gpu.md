# Piz Daint (GPU) cooldown and running instructions

**Required modifications** in <tt>cooldown.py</tt>:

	batch_tmpl = "job_daint_gpu.slm"  # Add path if necessary.
	run_command = "srun -n $SLURM_NTASKS --ntasks-per-core=$SLURM_NTASKS_PER_CORE --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK"

Specify **resources** (number of nodes and walltime) in <tt>cooldown.py</tt>:

	nodes = ...
	walltime = ...

On Piz Daint (GPU) we recommend using **1 MPI task per node**, setting **OMP\_NUM\_THREADS = 12**, and running **22 threads** in the Monte Carlo solver to leave space for the OS. Accordingly, `"walkers"` and `"accumulators"` need to be set in `input_sp.json.in` and `input_tp.json.in` such that **walkers + accumulators = 22**, e.g.

	"walkers": 8
	"accumulators": 14