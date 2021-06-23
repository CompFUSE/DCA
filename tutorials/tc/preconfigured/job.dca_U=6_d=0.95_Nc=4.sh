#!/bin/bash

#This will run all the tc calculations

#uncomment the runner apropriate to your system or add one.
#CADES
#RUNNER="srun -n1 --cpus-per-task=18 --gpus-per-task=1 --accel-bind=gv"

RUN_DCA="${RUNNER} ../../applications/dca/main_dca"
date

$RUN_DCA ./T=1/input_sp.json
$RUN_DCA ./T=0.75/input_sp.json
$RUN_DCA ./T=0.5/input_sp.json
$RUN_DCA ./T=0.25/input_sp.json
$RUN_DCA ./T=0.125/input_sp.json
$RUN_DCA ./T=0.1/input_sp.json
$RUN_DCA ./T=0.1/input_tp.json
$RUN_DCA ./T=0.09/input_sp.json
$RUN_DCA ./T=0.09/input_tp.json
$RUN_DCA ./T=0.08/input_sp.json
$RUN_DCA ./T=0.08/input_tp.json
$RUN_DCA ./T=0.07/input_sp.json
$RUN_DCA ./T=0.07/input_tp.json

date
