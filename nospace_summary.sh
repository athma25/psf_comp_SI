#!/bin/bash

# Input batch ID

if [ $# -ne 3 ]
then
	echo Give ID, runs and rep
	exit 1
fi

dir="./output/$1"

if [ ! -d ${dir} ]
then
	echo Directory does not exists. Copy output from seawulf.
	exit 2
fi

Rscript -e "source('KmeansGapLib.R'); for(i in 0:$2){ kmeans_gap('$1',i,$3)};" > ${dir}/kmeans.log 2>&1
llgal -s -f --tx 300 --ty 300 -w 6 -d ${dir} > ${dir}/gallery.log 2>&1

exit 0
