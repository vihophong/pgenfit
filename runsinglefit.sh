#!/bin/bash
./genFunctionCopy.sh $1
make
./main $1 $2 $3 $4 $5 $6
