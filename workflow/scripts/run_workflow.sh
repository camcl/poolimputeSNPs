#!/bin/bash

apptainer run --bind /usr/bin/sbatch,/usr/lib64/slurm,/etc/slurm,/run/munge,/usr/lib64/libmunge.so.2 container.sif -c 16 --latency-wait=5
