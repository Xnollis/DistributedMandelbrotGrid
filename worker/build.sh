#!/bin/sh
THE_OUTPUT_DIR=output
rm -Rf $THE_OUTPUT_DIR
mkdir $THE_OUTPUT_DIR
cd $THE_OUTPUT_DIR
cmake ..
make
cd ..
