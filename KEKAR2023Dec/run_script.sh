#!/bin/bash

# Set the path to your ROOT installation
#ROOTSYS=/Users/haein/software/root/6.22.08


# Set the path to your compiled C++ program
CPP_PROGRAM=./analysis_basic

# Loop through specific numbers
for number in {269..531}
do
    echo "Running for specific number $number"
    # Set up ROOT environment
    #source $ROOTSYS/bin/thisroot.sh
    
    # Run your compiled C++ program with the specific number
    
    $CPP_PROGRAM $number
    
    echo "---------------------------------------"
done
