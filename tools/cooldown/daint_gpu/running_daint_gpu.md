# Piz Daint (GPU) cooldown and running instructions

1. Set `"walkers"` and `"accumulators"` in `input_sp.json.in` and `input_tp.json.in` such that **walkers + accumulators = 22**, e.g.:

        "walkers": 8
        "accumulators": 14

2. **Required modifications** in <tt>cooldown.py</tt>:

    * `batch_tmpl = "job_daint_gpu.slm"` (Add path if necessary.)
    * `run_command_dca = "srun -n $SLURM_NTASKS --ntasks-per-core=$SLURM_NTASKS_PER_CORE --ntasks-per-node=$SLURM_NTASKS_PER_NODE -c $SLURM_CPUS_PER_TASK"`

3. Specify **resources**, i.e. number of nodes and walltime, in the generated batch scripts before submitting them:
    
        #SBATCH --nodes= ...
        #SBATCH --time= ...
