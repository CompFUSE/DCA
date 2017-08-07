# Local machine cooldown and running instructions

**Required modifications** in <tt>job_local.sh</tt>:

	mpirun=...
	export OMP_NUM_THREADS=...  # optional

**Required modifications** in <tt>cooldown.py</tt>:

	batch_tmpl = "job_local.sh"  # Add path if necessary.
	run_command = "$mpirun -n " + str(tasks)

Specify **resources** (number of MPI tasks) in <tt>cooldown.py</tt>:

	tasks = ...
