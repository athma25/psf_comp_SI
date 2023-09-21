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
	echo Directory does not exists. Generate parameter files first.
	exit 2
fi

run=$3

cp parGen.sh nospace.m ${dir}
matlab.exe -batch "for i=0:$2; for j=0:$((run-1)); nospace('$1',i,j,'off'); end; end; exit" < /dev/null > ${dir}/matlab.log 2>&1
Rscript -e "source('KmeansGapLib.R'); for(i in 0:$2){ kmeans_gap('$1',i,$3)};" > ${dir}/kmeans.log 2>&1
llgal -s -f --tx 300 --ty 300 -w 6 -d ${dir} > ${dir}/gallery.log 2>&1

exit 0
