#!/bin/bash

# Input batch ID
if [ $# -ne 1 ]
then
	echo Give all inputs
	exit 1
fi

dir="./output/$1"

if [ ! -d ${dir} ]
then
	mkdir ${dir}
else
	echo Directory exists. Backup and clear
	exit 2
fi

# Constant (across the sweep) parameters
l=3
m=13

et0=0

nu=1

t0=0
tf=10

# Create full factorial list of parameter files based on the following ranges
s0=0
et=0
sg_c=(0.2 0.25 0.3 0.35 0.4 0.45 0.5)
sg_s=(1 2)
imm=0.001

n0=${#s0[@]}
n1=${#et[@]}
n2=${#sg_c[@]}
n3=${#sg_s[@]}
n4=${#imm[@]}

for i0 in $(seq 0 $((n0-1)))
do
	for i1 in $(seq 0 $((n1-1)))
	do
		for i2 in $(seq 0 $((n2-1)))
		do
			for i3 in $(seq 0 $((n3-1)))
			do
				for i4 in $(seq 0 $((n4-1)))
				do
					fname=$(printf "run%03d.txt" $((i0+n0*i1+n0*n1*i2+n0*n1*n2*i3+n0*n1*n1*n2*n3*i4)))
					echo -e "\"Model parameters\" \t \"\" \t \"\"" > $fname
					echo -e "l \t $l \t \"Soil/resource niche size\"" >> $fname
					echo -e "m \t $m \t \"Lattice size in niche space\"" >> $fname
					echo -e "s0 \t ${s0[i0]} \t \"Initial soil condition\"" >> $fname
					echo -e "nu \t ${nu} \t \"Sensitivity to competition\"" >> $fname
					echo -e "et \t ${et[i1]} \t \"Soil conditioning strength\"" >> $fname
					echo -e "sg_c \t ${sg_c[i2]} \t \"Resource niche width\"" >> $fname
					echo -e "sg_s \t ${sg_s[i3]} \t \"Soil niche width\"" >> $fname
					echo -e "et0 \t ${et0} \t \"Exogenous driver strength\"" >> $fname
					echo -e "imm \t ${imm[i4]} \t \"Immigration rate\"" >> $fname
					echo -e "\n \"Numerics - parameters\" \t \"\" \t \"\"" >> $fname
					echo -e "t0\t ${t0} \t \"Initial time\"" >> $fname
					echo -e "tf \t ${tf} \t \"Final time\"" >> $fname

					echo -e "${s0[i0]} \t ${et[i1]} \t ${sg_c[i2]} \t ${sg_s[i3]} \t ${imm[i4]}" >> parSweep.txt
					
				done
			done
		done
	done
done

mv run[0-9]*.txt ${dir}

ls ${dir} > sweep.txt
mv sweep.txt parSweep.txt ${dir}

echo Generated parameter files at ${dir}
echo Last file is $fname
exit 0
