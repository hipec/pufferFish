#!/bin/bash

LOG_DIRECTORY="$PUFFERFISH/logs"
#Make sure this path is same in the scrtip test/teams/experiment.sh
BUILDS=( "sequential" )
BENCHMARKS=( "CilkSort" "SOR" "SORir" "LULESH" "LULESHir" "SRAD" "SRADir" )
ITERATIONS=10
#############################################
######### NO MODIFICATIONS BELOW ############
#############################################

mkdir -p $LOG_DIRECTORY 2>/dev/null
launch() {
    config=$1
    echo "=============Launching experiment: $config==============="
        current_iteration=0
        while [ $current_iteration -lt $ITERATIONS ]; do
            for exe in "${BENCHMARKS[@]}"; do
                FILE="$LOG_DIRECTORY/$exe.3096.1024.$config.threads-0.log"
		date
		who
                echo "Currently Running: $FILE "
                ./$exe 2>&1 | tee out
                if [ `cat out | grep "TEST PASSED" | wc -l` -eq 0 ]; then
                    echo "ERROR: $FILE did not give Success. Not appending result..."
                    #echo "ERROR: $FILE did not give Success. Terminating... Recheck and launch again"
                    #exit
                else
                    echo "Test Success: $FILE"
                    cat out >> $FILE
                fi
            done
            current_iteration=`expr $current_iteration + 1`
        done
}

run_experiment() {
    for config in "${BUILDS[@]}"; do
	./clean.sh
        make -f makefile.sequential 
	name="$config$NAME_PREFIX"
	launch $name
    done
}

run_experiment

