# Titan cooldown and running instructions

1. Set `"walkers"` and `"accumulators"` in `input_sp.json.in` and `input_tp.json.in` such that **walkers + accumulators = 8**, e.g.:

        "walkers": 3
        "accumulators": 5

2. **Required modifications** in <tt>cooldown.py</tt>:

    * `batch_tmpl = "job_titan.pbs"` (Add path if necessary.)
    * `run_command_dca = "aprun -n <num-tasks-dca> -N 2 -d 8"`
    * `run_command_analysis = "aprun -n <num-tasks-analysis> -N 2 -d 8"`  

    Replace `<num-tasks-dca>` and `<num-tasks-analysis>` with twice the number of nodes you intend to use for the **dca** and **analysis** runs, respectively.


3. Specify **resources**, i.e. number of nodes and walltime, in the generated batch scripts before submitting them:
    
        #PBS -l nodes= ...
        #PBS -l walltime= ...
