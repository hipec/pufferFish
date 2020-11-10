#!/bin/bash

rm -rf $INSTALL_DIRECTORY 2>/dev/null
cd $PUFFERFISH/runtime/hclib
./clean.sh
rm out 2>/dev/null
cd $PUFFERFISH/benchmarks
./clean.sh
rm out *.log 2>/dev/null
rm -rf $PUFFERFISH/logs* 2>/dev/null

