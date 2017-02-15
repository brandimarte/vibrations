#!/bin/bash

#########################################################################
# Job submission requested parameters, see 'man sbatch'.                #
#########################################################################
#SBATCH --job-name "i-smeagol"
#SBATCH --nodes 1 
#SBATCH --ntasks 32
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node 32
#SBATCH --mem 48000              ## MB
##SBATCH --mem-per-cpu 4096      ## MB
#SBATCH --time 160:00:00         ## min, min:sec, h:min:sec, d-h:min:sec
#SBATCH --partition medium       ## short, medium, long, verylong
#SBATCH --mail-type ALL          ## BEGIN, END and FAIL
#SBATCH --mail-user brandimarte@gmail.com
#SBATCH --output ISmeagol.log
#SBATCH --error ISmeagol.err
##SBATCH --dependency=afterany:job_id[:jobid...]

#########################################################################
# Setup environment.                                                    #
#########################################################################

ulimit -s unlimited
ulimit -l unlimited

#export OMP_NUM_THREADS=1

module load i-smeagol/master-impi-13.1.1

##ISMEAGOL=/home/pedro/local/bin/i-smeagol.2.0.dev.CDF.IMPI.i2013
##ISMEAGOL=/home/pedro/local/bin/i-smeagol.2.0.optim.CDF.IMPI.i2013
##ISMEAGOL=i-smeagol.2.0.dev.CDF.IMPI.i2013
ISMEAGOL=i-smeagol.2.0.optim.CDF.IMPI.i2013

#########################################################################
# Setup and copy files to computing node, '/home/scratch'.              #
#########################################################################

# Remove annoying or confusing files.
##rm ${SLURM_SUBMIT_DIR}/INPUT_TMP.*

LSCRATCH="/home/scratch/${USER}/ismeagol_${SLURM_JOB_ID}"
RESULT_DIR=${SLURM_SUBMIT_DIR}

# Generate information files in original folder.
echo `scontrol show hostnames ${SLURM_NODELIST}` > \
    NODELIST_${SLURM_JOB_ID}
echo ${LSCRATCH} > LSCRATCH_${SLURM_JOB_ID}

# Make lscratch.
srun --ntasks ${SLURM_JOB_NUM_NODES} --ntasks-per-node=1 \
    mkdir -p ${LSCRATCH}
echo ${SLURM_SUBMIT_DIR} > ${LSCRATCH}/orig.workdir

# Copy data to lcratch in *every* node.
srun --ntasks ${SLURM_JOB_NUM_NODES} --ntasks-per-node=1 \
    rsync -av --exclude 'slurm-*' ${SLURM_SUBMIT_DIR}/* ${LSCRATCH}
##cp ${SLURM_SUBMIT_DIR}/{INCAR,KPOINTS,POSCAR,POTCAR} ${LSCRATCH}

# Switch to working directory.
cd ${LSCRATCH}

echo "------------------------------------------------------------------"
echo "Beginning of the job:" `date`
echo "Hostname: " `hostname`
echo "PWD: " ${PWD}

srun -n ${SLURM_NTASKS} ${ISMEAGOL} < bdtVibra.fdf > bdtVibra.out

echo "End of the job:" `date`
echo "------------------------------------------------------------------"

# Cleanup, sleep, then copy back to original folder.
cp -a ./* ${RESULT_DIR}
cd ${SLURM_SUBMIT_DIR}
srun --ntasks ${SLURM_JOB_NUM_NODES} --ntasks-per-node=1 \
    rm -rf ${LSCRATCH}

# Remove files with no more use.
rm ${SLURM_SUBMIT_DIR}/LSCRATCH_${SLURM_JOB_ID}
rm ${SLURM_SUBMIT_DIR}/NODELIST_${SLURM_JOB_ID}

############################## END OF JOB ###############################

