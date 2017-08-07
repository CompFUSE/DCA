# Titan cooldown and running instructions

On Titan we recommend using **2 MPI tasks per node**, each task running **8 threads**. The following instructions are based on this recommendation.

**Required modifications** in <tt>cooldown.py</tt>:

	tasks_per_node = 2
	tasks = nodes*tasks_per_node
	threads = 8
	batch_tmpl = "job_titan.pbs"  # Add path if necessary.
	run_command = "aprun -n " + str(tasks) + " -N " + str(tasks_per_node) + " -d " + str(threads)

Specify **resources** (number of nodes and walltime) in <tt>cooldown.py</tt>:

	nodes = ...
	walltime = ...

Preset in <tt>job_titan.slm</tt>:

	export OMP_NUM_THREADS=8

`"walkers"` and `"accumulators"` need to be set in `input_sp.json.in` and `input_tp.json.in` such that **walkers + accumulators = 8**, e.g.

	"walkers": 3
	"accumulators": 5