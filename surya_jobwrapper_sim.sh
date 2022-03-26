#!/bin/bash

source env.sh

mkdir temp

# ./ino_digi_read1_sim corsika76300_FLUKA_SIBYLL_3dFlux_20211214at_20220107aa_iter4.root $1 $2 $3
./ino_digi_read1_sim $1 $2 $3 $4

