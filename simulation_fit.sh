#!/bin/bash
>simoutput_fit.root.txt
for ((i=1; i<=$1; i++)); do
	./simulation parmsex.txt simparmsex.txt simoutput.root
	./runsinglefitT12.sh parmsex.txt simoutput.root simoutput_fit.root 0.005 500 0 parmsex.txt_effparms 4357
done


