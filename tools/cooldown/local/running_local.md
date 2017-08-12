# Local machine cooldown and running instructions

**Required modifications** in <tt>cooldown.py</tt>:

* `batch_tmpl = "job_local.sh"` (Add path if necessary.)
* `run_command_dca`: e.g. `"mpirun -n 8"`
* `run_command_analysis`: e.g. `"mpirun -n 1"`
