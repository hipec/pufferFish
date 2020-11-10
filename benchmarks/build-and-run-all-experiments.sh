#!/bin/bash

# THIS IS A SHORTCUT SCRIPT TO BUILD AND RUN
# ALL EXPERIMENTS IN A SINGLE GO!
# PLOTTING THE RESULTS IS NOT INCLUDED HERE.

cd $PUFFERFISH
./clean.sh
cd $PUFFERFISH/runtime
./patch-hclib.sh
./build-hclib.sh
cd $PUFFERFISH/benchmarks
./experiment-seq.sh
./experiment.sh
cd $PUFFERFISH/runtime
./build-likwid-hclib.sh
cd $PUFFERFISH/benchmarks
./experiment-pcm.sh
