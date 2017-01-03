# Titan cooldown and running instructions

On Titan we recommend running **2 MPI tasks per node** each running **8 threads**. The following instructions are based on this recommendation.

**Required modifications** in <tt>cooldown.py</tt>:

	tasks_per_node = 2
	tasks = nodes*tasks_per_node
	threads = 8
	batch_tmpl = "job_titan.pbs"  # Add path if necessary.
	run_command = "aprun -n " + str(tasks) + " -N " + str(tasks_per_node) + " -d " + str(threads)

Adjust **resources** (number of nodes and walltime) in <tt>cooldown.py</tt> if necessary:

	nodes = ...
	walltime = ...

Preset in <tt>job_titan.slm</tt>:

	export OMP_NUM_THREADS=8

Set *walkers* and *accumulators* in the input files such that ***walkers* + *accumulators* = threads = 8**, e.g.

	"walkers": 3
	"accumulators": 5