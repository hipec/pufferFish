
#!/bin/bash

#USE_PCM="--enable-likwid"

#############################################
######### DONT CHANGE ANYTHING BELOW ########
#############################################

verify_installation() {
  DIRNAME=$1
  if [ `cat out | grep "Libraries have been installed in" | wc -l` -eq 0 ]; then
    echo "INSTALLATION ERROR: $DIRNAME"
    exit
  fi
  rm -rf $INSTALL_DIRECTORY/$DIRNAME
  cp -rf hclib-install $INSTALL_DIRECTORY/$DIRNAME
  echo "export HCLIB_ROOT=$INSTALL_DIRECTORY/$DIRNAME" > $INSTALL_DIRECTORY/$DIRNAME/bin/hclib_setup_env.sh
  chmod +X $INSTALL_DIRECTORY/$DIRNAME/bin/hclib_setup_env.sh
}

rm -rf $INSTALL_DIRECTORY 2>/dev/null
mkdir -p $INSTALL_DIRECTORY 2>/dev/null

cd hclib

#Default implementation with random work-stealing
./clean.sh
HCLIB_FLAGS=" --enable-randws ${USE_PCM}" ./install.sh  2>&1 | tee out
verify_installation "defaultRand"

#Default implementation with HPT_DA
./clean.sh
HCLIB_FLAGS=" --enable-hptDA ${USE_PCM}" ./install.sh  2>&1 | tee out
verify_installation "hptDA"

#ELASTIC TASKING IN PUFFERFISH
./clean.sh
HCLIB_FLAGS=" --enable-pufferFish ${USE_PCM}" ./install.sh  2>&1 | tee out
verify_installation "pufferFish"
