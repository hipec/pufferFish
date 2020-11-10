#!/bin/bash

#Make sure this path is same in the scrtip test/teams/experiment.sh
LOG_DIRECTORY="$PUFFERFISH/logs"
HPT_FILE_DIRECTORY="$PUFFERFISH/runtime/hclib/hpt"
NAME_PREFIX=
COMMAND=""
BUILDS=( "defaultRand" "hptDA" "cilk" "pufferFish" )
BENCHMARKS=( "CilkSort" "SOR" "SORir" "LULESH" "LULESHir" "SRAD" "SRADir" )
TIMEOUT=( 40 40 40 40 40 40 40 )
TIMED_ITERATIONS=10
DISCARD_ITERATIONS=2
THREADS=( 16 32 )
#############################################
export HCLIB_PERFCOUNTERS="RAPL_PKG_ENERGY:PWR1,L3_MISS:CPMC0,L3_ACCESS:CPMC1,RETIRED_INSTRUCTIONS:PMC0,DATA_CACHE_ACCESSES:PMC2,DATA_CACHE_REFILLS_ALL:PMC3,CPU_CLOCKS_UNHALTED:PMC1"
export HCLIB_STATS=1
export HCLIB_BIND_THREADS=1
#############################################
######### NO MODIFICATIONS BELOW ############
#############################################

mkdir -p $LOG_DIRECTORY 2>/dev/null

launch() {
    config=$1
    echo "=============Launching experiment: $config==============="
    NET_ITERATIONS=`expr $DISCARD_ITERATIONS + $TIMED_ITERATIONS`
    hpt_index=0
    for thread in "${THREADS[@]}"; do
        if [ `echo $config | grep "hptDA" | wc -l` -eq 1 ]; then
	    N=`expr $thread / 8`
            export HCLIB_HPT_FILE="$HPT_FILE_DIRECTORY/hpt-hippo-"$N"numa.xml"
	    unset HCLIB_WORKERS
	    hpt_index=`expr $hpt_index + 1`
        elif [ `echo $config | grep "pufferFish" | wc -l` -eq 1 ]; then
	    N=`expr $thread / 8`
            export HCLIB_HPT_FILE="$HPT_FILE_DIRECTORY/hpt-hippo-"$N"numa.xml"
	    unset HCLIB_WORKERS
	    hpt_index=`expr $hpt_index + 1`
        elif [ `echo $config | grep "default" | wc -l` -eq 1 ]; then
	    unset HCLIB_HPT_FILE
	    export HCLIB_WORKERS=$thread
        elif [ `echo $config | grep "cilk" | wc -l` -eq 1 ]; then
	    export CILK_NWORKERS=$thread
	else
	    echo "ERROR: UNSUPPORTED CONFIGURATION"
	    exit
	fi
	timeout_index=0
        for exe in "${BENCHMARKS[@]}"; do
            current_iteration=0
            while [ $current_iteration -lt $NET_ITERATIONS ]; do
                FILE="$LOG_DIRECTORY/$exe.3096.1024.$config.threads-$thread.log"
		timeout=${TIMEOUT[$timeout_index]}
		date
		who
                echo "Currently Running: $FILE as follows:"
		echo "./timedrun -t $TIMEOUT $COMMAND ./$exe"
                ./timedrun -t $TIMEOUT $COMMAND ./$exe 2>&1 | tee out
		if [ $current_iteration -ge $DISCARD_ITERATIONS ]; then
                    if [ `cat out | grep "TEST PASSED" | wc -l` -eq 0 ]; then
                        echo "ERROR: $FILE did not give Success. Not appending result..."
                        #echo "ERROR: $FILE did not give Success. Terminating... Recheck and launch again"
                        #exit
                    else
                        echo "Test Success: $FILE"
                        cat out >> $FILE
                    fi
		fi
                current_iteration=`expr $current_iteration + 1`
            done
        done
	timeout_index=`expr $timeout_index + 1`
    done
}

run_experiment() {
    for config in "${BUILDS[@]}"; do
        source $INSTALL_DIRECTORY/$config/bin/hclib_setup_env.sh
	./clean.sh
        if [ `echo $config | grep "default" | wc -l` -eq 1 ]; then
            make -f makefile.randws
        elif [ `echo $config | grep "cilk" | wc -l` -eq 1 ]; then
            make -f makefile.cilk
        elif [ `echo $config | grep "pufferFish" | wc -l` -eq 1 ]; then
            make -f makefile.pufferFish
	else
            make -f makefile.hptDA
	fi
	name="$config$NAME_PREFIX"
	launch $name
    done
}

#First touch except for hpt configs
run_experiment

#Interleaved experiments except for hpt configs
BUILDS=( "defaultRand" "cilk" )
NAME_PREFIX="iAll"

THREADS=( 16 )
COMMAND="numactl -i 0,1"
run_experiment

THREADS=( 32 )
COMMAND="numactl -iall"
run_experiment

