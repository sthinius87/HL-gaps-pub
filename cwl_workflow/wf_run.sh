#!/bin/bash
export TOIL_SLURM_ARGS="-t 100-00:00:00 --export=ALL,PATH"
#"-t 100-00:00:00 --export=ALL,PATH"

#source /home/sat/anaconda3/etc/profile.d/conda.sh
#conda activate py36toil
#source activate py36toil


#toil-cwl-runner --workDir ./tmp --jobStore ./job_store --batchSystem slurm --disableCaching true --disableProgress cat_grep.cwl params.yml > wf.log 2>&1 &
# --preserve-entire-environment

#toil-cwl-runner --workDir ./tmp --jobStore ./job_store --batchSystem slurm --disableCaching true --disableProgress --restart --clean never scatter_wf.cwl scatter_inp.yml > wf.log 2>&1 &


#cwltool --preserve-entire-environment scatter_wf.cwl scatter_inp.yml

#NORMAL
toil-cwl-runner --workDir ./tmp \
--jobStore ./job_store \
--debug \
--batchSystem slurm \
--caching FALSE \
--defaultDisk 1G \
--defaultMemory 1G \
--disableProgress TRUE \
--clean never \
--cleanWorkDir always \
--maxNodes 12 \
--maxJobs 500 \
--defaultCores 4 \
--stats \
scatter_wf.cwl \
scatter_inp.yml \

#RESTART
#toil-cwl-runner --workDir ./tmp \
#--jobStore ./job_store \
#--debug \
#--batchSystem slurm \
#--caching FALSE \
#--defaultDisk 1G \
#--defaultMemory 1G \
#--disableProgress TRUE \
#--clean never \
#--cleanWorkDir always \
#--maxNodes 12 \
#--maxJobs 500 \
#--defaultCores 4 \
#--stats \
#--restart \
#scatter_wf.cwl \
#scatter_inp.yml \




#--disableCaching TRUE \
#--stats \
# always,onError,never,onSuccess
#--preserve-entire-environment \
# --runCwlInternalJobsOnWorkers true \
#--cleanWorkDir always,onError,never,onSuccess

#toil-cwl-runner --workDir ./tmp \
#--jobStore ./job_store \
#--debug \
#--batchSystem slurm \
#--disableCaching true \
#--defaultDisk 1M \
#--defaultMemory 1M \
#--disableProgress \
#--clean never \
#--cleanWorkDir always \
#--maxNodes 12 \
#--defaultCores 6 \
#--stats \
#--restart \
#scatter_wf.cwl \
#scatter_inp.yml \
