#!/bin/bash
./genFunctionCopy.sh $1
make
nn=$(expr $6 / $8)
hh=""
for i in $(seq 1 $8)
do
	echo "seed = $(expr $i + 4356), number of MC samples = $nn"
        ./main $1 $2 $3_$i $4 $5 $nn $7 $(expr $i + 4356) &
        hh+=" $3_$i"
done
wait
echo $hh
hadd -f $3 $hh

for i in $(seq 2 $8)
do
    /bin/rm -rf $3_$i
    /bin/rm -rf $3_$i.txt
done
/bin/rm -rf $3_1
