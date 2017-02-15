#!/bin/bash

#########################################################################
# Job submission requested parameters, see 'man sbatch'.                #
#########################################################################
#SBATCH --job-name "vibrations"
#SBATCH --nodes 1 
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 3000             ## MB
##SBATCH --mem-per-cpu 3000    ## MB
#SBATCH --time 24:00:00        ## min,  min:sec,  h:min:sec, d-h:min:sec
#SBATCH --partition short      ## short, medium, long, verylong
#SBATCH --mail-type ALL        ## BEGIN, END and FAIL
#SBATCH --mail-user brandimarte@gmail.com
#SBATCH --output Vibrations.log
#SBATCH --error Vibrations.err
##SBATCH --dependency=afterany:job_id[:jobid...]

#########################################################################
# Setup environment.                                                    #
#########################################################################

ulimit -s unlimited
ulimit -l unlimited

#export OMP_NUM_THREADS=1

module load vibrations/master-intel-13.1.1

##VIBRATIONS=/home/pedro/local/bin/vibrations
VIBRATIONS="vibrations"

LSCRATCH="/home/scratch/${USER}/vibrations_${SLURM_JOB_ID}"
RESULT_DIR=${SLURM_SUBMIT_DIR}

#########################################################################
# Setup and copy files to computing node, '/home/scratch'.              #
#########################################################################

# Generate information files in original folder.
echo `scontrol show hostnames ${SLURM_NODELIST}` > \
    NODELIST_${SLURM_JOB_ID}
echo ${LSCRATCH} > LSCRATCH_${SLURM_JOB_ID}

# Make lscratch.
mkdir -p ${LSCRATCH}
echo ${SLURM_SUBMIT_DIR} > ${LSCRATCH}/orig.workdir

# Copy data to lcratch in *every* node.
rsync -av --exclude 'slurm-*' ${SLURM_SUBMIT_DIR}/* ${LSCRATCH}

# Switch to working directory.
cd ${LSCRATCH}

echo "------------------------------------------------------------------"
echo "Beginning of the job:" `date`
echo "Hostname: " `hostname`
echo "PWD: " ${PWD}

${VIBRATIONS} . bdtVibra.fdf full splitFC > vibrations.out

echo "End of the job:" `date`
echo "------------------------------------------------------------------"

# Cleanup, sleep, then copy back to original folder.
cp -a ./*JMOL* ./*.Meph ./vibrations.out ${RESULT_DIR}
cd ${SLURM_SUBMIT_DIR}
rm -rf ${LSCRATCH}

# Remove files with no more use.
rm ${SLURM_SUBMIT_DIR}/LSCRATCH_${SLURM_JOB_ID}
rm ${SLURM_SUBMIT_DIR}/NODELIST_${SLURM_JOB_ID}

############################## END OF JOB ###############################
