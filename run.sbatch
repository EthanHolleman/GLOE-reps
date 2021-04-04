#!/bin/bash -login
#SBATCH -p high                      # partition, or queue, to assign to
#SBATCH -J GLOE_reps         # name for job
#SBATCH -o j%j.out # file to write STDOUT to
#SBATCH -e j%j.err # file to write error output to
#SBATCH -N 1                        # one "node", or computer
#SBATCH -n 1                        # one task for this node
#SBATCH -c 1                        # one cores per task
#SBATCH -t 48:00:00                  
#SBATCH --mem=10gb                  # ask for no more than 10 GB of memory
#SBATCH --mail-type=ALL
#SBATCH --mail-user=etholleman@ucdavis.edu

mkdir -p Logs
conda activate GLOE-reps
snakemake -j 14 -s Snakefile --cluster-config cluster.yml \
--cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} \
-n {cluster.cpus} --mem {cluster.mem} -J {cluster.name} -o {cluster.output} \
-e {cluster.output} --mail-type ALL --mail-user {cluster.email}" \
--conda-frontend=mamba \
--latency-wait 120 --verbose --use-conda -p -k $@