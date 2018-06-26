#!/bin/bash

date

./main_dca ./T=1/input_sp.json
./main_dca ./T=0.75/input_sp.json
./main_dca ./T=0.5/input_sp.json
./main_dca ./T=0.25/input_sp.json
./main_dca ./T=0.125/input_sp.json
./main_dca ./T=0.1/input_sp.json
./main_dca ./T=0.1/input_tp.json
./main_dca ./T=0.09/input_sp.json
./main_dca ./T=0.09/input_tp.json
./main_dca ./T=0.08/input_sp.json
./main_dca ./T=0.08/input_tp.json
./main_dca ./T=0.07/input_sp.json
./main_dca ./T=0.07/input_tp.json

date
