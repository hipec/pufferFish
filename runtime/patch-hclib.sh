#!/bin/bash

rm -rf hclib 2>/dev/null
git clone https://github.com/habanero-rice/hclib.git
cd hclib
git checkout ab310a0
patch -p1 <../hclib.patch
