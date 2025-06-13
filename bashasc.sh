#!/bin/bash

for i in {1..19}
do
  echo "Submit $i"
  bsub -q multicore40 -R "rusage[mem=75000]" -M 80000 "Rscript asc.R $i"
done

echo "all jobs submitted"

