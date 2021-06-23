#!/bin/bash

#uncomment the runner apropriate to your system or add one.
#CADES
#RUNNER="srun -n1 --cpus-per-task=18 --gpus-per-task=1 --accel-bind=gv"

RUN_DCA="${RUNNER} ../../applications/analysis/main_analysis"

date

$RUN_DCA ./T=0.1/input_tp.json
$RUN_DCA ./T=0.09/input_tp.json
$RUN_DCA ./T=0.08/input_tp.json
$RUN_DCA ./T=0.07/input_tp.json

date
