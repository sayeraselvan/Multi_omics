#!/bin/bash

# Set the queue and resource constraints
queue="multicore40"
mem_per_task=25000
total_memory=50000

for i in {1..19}; do
    bio_number="bio_$i"

    # Command to execute Rscript with bsub
    cmd="bsub -q $queue -R 'rusage[mem=$mem_per_task]' -M $total_memory \"Rscript /netscratch/dep_coupland/grp_fulgione/siva/scripts/bufferfrom143_1.R /netscratch/dep_coupland/grp_fulgione/siva/worldclim/wc2.1_30s_bio/tif/wc2.1_30s_${bio_number}.tif ${bio_number}.asc\""

    # Execute the command
    eval "$cmd"
done


