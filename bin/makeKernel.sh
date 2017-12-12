#!/bin/sh

# .sif file <source> <int> <target>
input_network=$1
time_param=$2
output=$3
cut -f 1,3 $input_network \
    | `dirname $0`/order_keys.py \
    | sort -k 1,2 -t '	' -u \
    > /tmp/in.tab

sed -e "s/%TIME%/$time_param/g" \
	-e "s/%OUTPUT%/$output/g" \
	`dirname $0`/network_diffusion_kernel.m > startup.m

matlab -nojvm -nodesktop -nosplash;
rm -f startup.m 
