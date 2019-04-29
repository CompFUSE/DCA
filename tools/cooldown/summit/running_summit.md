# Summit cooldown and running instructions

1. Set `"walkers"` to the  number of physical cores available per process, usually 7, in 
   `input_sp.json.in` and `input_tp.json.in`. Usually 1 accumulator is fine.

        "walkers": 7
        "accumulators": 1

2. **Required modifications** in <tt>cooldown.py</tt>:

    * `batch_tmpl = "job_summit.pbs"` (Add path if necessary.)
    * `run_command_dca = "srun -n <num-tasks-dca>  -a 1 -g 1 -c 7 -b rs"`
    * `run_command_analysis = "srun -n <num-tasks-analysis>"`  

    Replace `<num-tasks-dca>` with the number of GPUs used for the DCA run, so 6 times the number of 
    compute nodes, and `<num-tasks-analysis>` with the number of processes for the analysis 
    application.


3. Specify **resources**, i.e. number of nodes and walltime, in the generated batch scripts before 
   submitting them:
   ``` 
   #BSUB -nnodes nodes= ...
   #BSUB -W walltime= ...
   ```
