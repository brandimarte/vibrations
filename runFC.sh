#!/bin/bash

#  *****************************************************  #
#             ** Phonon Vibration Analysis **             #
#                                                         #
#                    **  Version 2  **                    #
#                                                         #
#                         IF/USP                          #
#                                                         #
#   Advisor: Prof. Dr. Alexandre Reily Rocha              #
#                                                         #
#   Author: Pedro Brandimarte Mendonca                    #
#                                                         #
#   File: runFC.sh                                        #
#                                                         #
#   Versions: 1 - 10/10/2012                              #
#             2 - 08/01/2013                              #
#                                                         #
#  *****************************************************  #
#  This script receives an input 'fdf' file for a force   #
#  constants (FC) run, includes the flag 'PB.FCwriteHS'   #
#  and submits the calculation. Then, it changes the      #
#  'fdf' file, removing the flag 'PB.FCwriteHS',          #
#  including 'PB.OnlyS' and changing the geometry and     #
#  number of atoms.                                       #
#  *****************************************************  #
#  Input:                                                 #
#    ${1} :  mpi compiler (e.g. 'mpirun')                 #
#    ${2} :  file mapping procs (e.g. '${PBS_NODEFILE}')  #
#    ${3} :  SIESTA executable with path                  #
#    ${4} :  FC calculation directory                     #
#  *****************************************************  #

# Prints the header.
echo -e ""
echo -e "**  *****************************************************  **"
echo -e "**             ** Phonon Vibration Analysis **             **"
echo -e "**                                                         **"
echo -e "**                    **  Version 1  **                    **"
echo -e "**                                                         **"
echo -e "**                         IF/USP                          **"
echo -e "**                                                         **"
echo -e "**   Advisor: Prof. Dr. Alexandre Reily Rocha              **"
echo -e "**                                                         **"
echo -e "**   Author: Pedro Brandimarte Mendonca                    **"
echo -e "**                                                         **"
echo -e "**  *****************************************************  **"

# Checks if the number of arguments is correct.
if [ ${#} != 4 ]
then
    echo -e "\nvibranal: ERROR: wrong number of arguments!\n"
    echo -e "vibranal: Use: ./runFC.sh [mpi compiler] [file mapping procs] \ "
    echo -e "                 [SIESTA executable] [FC calculation directory]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "vibranal: Start of run: "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy

# Number of cores.
cores=$[ `cat ${2} | wc -l` ]

# FC and work directories.
check=`echo "${4}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    FCdir=${4}
else
    FCdir=${4}/
fi

# Checks if the files and FC folder exists and are accessible.
echo -e ""
echo -n "vibranal: Checking input... "
if [ ! -r ${2} ]
then
    echo -e "\nvibranal: ERROR: the file \"${2}\" doesn't exist or is not accessible!\n"
    exit -1
elif [ ! -r ${3} ]
then
    echo -e "\nvibranal: ERROR: the file \"${3}\" doesn't exist or is not accessible!\n"
    exit -1
elif [ ! -x ${3} ]
then
    echo -e "\nvibranal: ERROR: you don't have permission to execute \"${3}\"!\n"
    exit -1
elif [ ! -r ${FCdir} ]
then
    echo -e "\nvibranal: ERROR: the directory \"${FCdir}\" doesn't exist or is not accessible!\n"
    exit -1
else
    echo -e "ok!\n"
fi

# Gets '.fdf' files and creates 'pbFCrun.fdf'.
if [ -r ${FCdir}pbFCrun.fdf ]
then
    # Removes old 'pbFCrun.fdf' file.
    rm ${FCdir}pbFCrun.fdf
fi
if [ -r ${FCdir}pbCoord.fdf ]
then
    # Removes old 'pbCoord.fdf' file.
    rm ${FCdir}pbCoord.fdf
fi
if ls ${FCdir}*.fdf > /dev/null 2>&1
then
    > ${FCdir}FCfdf.tmp
    FCfdf=`ls ${FCdir}*.fdf`
    echo -e "vibranal: FC fdf input files found:"
    for fdf in ${FCfdf}
    do
	echo -e "\t\t\t\t\t${fdf}"
        # Copies to 'FCfdf.tmp' without comented lines.
	sed '/^[#]\+/d' ${fdf} >> ${FCdir}FCfdf.tmp
    done
else
    echo -e "\nvibranal: ERROR: couldn't find any FC fdf input file at \"${FCdir}\"!\n"
    exit -1
fi
cat ${FCdir}FCfdf.tmp > ${FCdir}pbFCrun.fdf
rm ${FCdir}FCfdf.tmp
FCfdf=${FCdir}pbFCrun.fdf

# Sets FC 'fdf' input file with the flag 'PB.FCwriteHS .true.'.
echo -e "\nvibranal: Including the option 'PB.FCwriteHS .true.' at 'pbFCrun.fdf' file.\n"
check=`grep -i "PB.FCWRITEHS" ${FCfdf}`
if [ "${check}" == "" ]
then
    echo -e "\nPB.FCwriteHS            .true." >> ${FCfdf}
else
    sed "s/${check}/PB.FCwriteHS            .true.\n/I" -i ${FCfdf}
fi

# Sets FC 'fdf' input file with the flag 'PB.FConlyS .false.'.
echo -e "vibranal: Including the option 'PB.FConlyS .false.' at 'pbFCrun.fdf' file.\n"
check=`grep -i "PB.FCONLYS" ${FCfdf}`
if [ "${check}" == "" ]
then
    echo -e "\nPB.FConlyS              .false." >> ${FCfdf}
else
    sed "s/${check}/PB.FConlyS              .false.\n/I" -i ${FCfdf}
fi

# Submits the job.
echo -e "vibranal: Running FC calculation...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *.alloc *.BONDS* *.STRUCT_* INPUT_TMP.*

# Sets FC 'fdf' input file with the flag 'PB.FCwriteHS .false.'.
echo -e "\nvibranal: Including the option 'PB.FCwriteHS .false.' at 'pbFCrun.fdf' file.\n"
check=`grep -i "PB.FCWRITEHS" ${FCfdf}`
sed "s/${check}/PB.FCwriteHS            .false.\n/I" -i ${FCfdf}

# Sets FC 'fdf' input file with the flag 'PB.FConlyS .true.'.
echo -e "vibranal: Including the option 'PB.FConlyS .true.' at 'pbFCrun.fdf' file.\n"
check=`grep -i "PB.FCONLYS" ${FCfdf}`
sed "s/${check}/PB.FConlyS              .true.\n/I" -i ${FCfdf}

# Gets the atoms displacement.
echo -n "vibranal: Atoms displacement: "
displ=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $2}'`
if [ "${displ}" == "" ]
then
    echo -e "\nvibranal: ERROR: can't find the 'MD.FCdispl' value at FC fdf input file!\n"
    exit -1
fi

# Gets the displacement unit (Bohr or Ang).
dunit=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $3}' | tr '[:upper:]' '[:lower:]'`
if [ "${dunit}" == "" ]
then
    echo -e "\nvibranal: ERROR: can't find the 'MD.FCdispl' unit at FC fdf input file!\n"
    exit -1
fi
echo -e "${displ} ${dunit}\n"

# Gets the atomic coordinates unit (Bohr or Ang).
aunit=`grep -i "ATOMICCOORDINATESFORMAT" ${FCfdf} | awk '{print $2}' \
    | tr '[:upper:]' '[:lower:]'`
if [ "${aunit}" == "" ]
then
    echo -e "\nvibranal: ERROR: can't find the 'AtomicCoordinatesFormat' at FC fdf input file!\n"
    exit -1
elif [ "${aunit}" == "notscaledcartesianang" ]
then
    aunit="ang"
elif [ "${aunit}" != "bohr" ]
then
    if [ "${aunit}" != "ang" ]
    then
	echo -e "\nvibranal: ERROR: unrecognized atomic coordinates format \"${aunit}\" at"
	echo -e "                 FC fdf input file! It must be \"Ang\" or \"Bohr\"!\n"
	exit -1
    fi
fi

# Converts the displacement unit if necessary.
echo -e "vibranal: Atomic coordinates unit: ${aunit}\n"
if [ "${dunit}" != "${aunit}" ]
then
    echo -e "vibranal: Converting the displacement to ${aunit}.\n"
    if [ "${dunit}" == "bohr" ]
    then
	displ=`echo "scale = 15; ${displ} * 0.52917721092" | bc`
    elif [ "${dunit}" == "ang" ]
    then
	displ=`echo "scale = 15; ${displ} / 0.52917721092" | bc`
    else
	echo -e "\nvibranal: ERROR: unrecognized displacement unit \"${dunit}\" at FC"
	echo -e "                 fdf input file! It must be \"Ang\" or \"Bohr\"!\n"
	exit -1
    fi
fi

# Gets the number of atoms.
natoms=`grep -i "NUMBEROFATOMS" ${FCfdf} | awk '{print $2}'`
if [ "${natoms}" == "" ]
then
    echo -e "\nvibranal: ERROR: can't find the 'NumberOfAtoms' at FC fdf input file!\n"
    exit -1
fi
echo -e "vibranal: Actual number of atoms: ${natoms}\n"

# Doubles the number of atoms.
echo -n "vibranal: Doubling the number of atoms from ${natoms} to "
check=`grep -i -m 1 "NUMBEROFATOMS" ${FCfdf}`
n2atoms=`echo "${natoms} * 2" | bc`
sed "s/${check}/NumberOfAtoms           ${n2atoms}\n/I" -i ${FCfdf}
echo -e "${n2atoms}\n"

# Gets and removes the 'AtomicCoordinatesAndAtomicSpecies' block.
echo -n "vibranal: Reading the block 'AtomicCoordinatesAndAtomicSpecies'... "
acoord=`grep -i -m 1 -A ${natoms} "ATOMICCOORDINATESANDATOMICSPECIES" ${FCfdf} \
    | tail -n ${natoms}`
if [ "${acoord}" == "" ]
then
    echo -e "\nvibranal: ERROR: can't find the block 'AtomicCoordinatesAndAtomicSpecies'"
    echo -e "                 at FC fdf input file!\n"
    exit -1
else
    echo -e "ok!\n"

    # Copies 'AtomicCoordinatesAndAtomicSpecies' block to 'pbCoord.tmp' file.
    echo -e "vibranal: Copying 'AtomicCoordinatesAndAtomicSpecies' to 'pbCoord.tmp' file.\n"
    grep -i -m 1 -A ${natoms} "ATOMICCOORDINATESANDATOMICSPECIES" ${FCfdf} | tail -n ${natoms} \
	> ${FCdir}pbCoord.tmp

    # Removes 'AtomicCoordinatesAndAtomicSpecies' block.
    echo -e "vibranal: Removing 'AtomicCoordinatesAndAtomicSpecies' from 'pbFCrum.fdf' file.\n"
    first=`nl -b a ${FCfdf} | grep -i -m 1 "ATOMICCOORDINATESANDATOMICSPECIES" \
	| awk '{print $1}'`
    last=`echo "${first} + ${natoms} + 1" | bc`
    sed "${first},${last}d" -i ${FCfdf}

    # Removes '%include' lines.
    echo -e "vibranal: Removing '%include' lines from 'pbFCrum.fdf' file.\n"
    check=`grep -i -m 1 "%INCLUDE" ${FCfdf}`
    while [ "${check}" != "" ]
    do
	sed "s/${check}/ /" -i ${FCfdf}	
	check=`grep -i -m 1 "%INCLUDE" ${FCfdf}`
    done
fi

# Add the external coordinate flag.
echo -e "vibranal: Adding '%include pbCoord.fdf' at 'pbFCrun.fdf' file.\n"
echo -e "%include                pbCoord.fdf" >> ${FCfdf}

# Gets the system label.
slabel=`grep -i "SYSTEMLABEL" ${FCfdf} | awk '{print $2}'`
if [ "${slabel}" == "" ]
then
    echo -e "\nvibranal: ERROR: can't find the 'SystemLabel' at FC fdf input file!\n"
    exit -1
fi

# Only S runs:

echo -e "vibranal: Including atomic coordinates with at 'pbCoord.fdf' file.\n"
echo -e "%block AtomicCoordinatesAndAtomicSpecies" > ${FCdir}pbCoord.fdf
awk '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2, $3, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf

# 1. '-x' move.
echo -e "vibranal: Changing the system label to '${slabel}_1'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_1/I" -i ${FCfdf}

echo -e "vibranal: Shifting half of the atoms to '-x' direction.\n"
awk -v var=${displ} '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1 - var, $2, $3, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies" >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibranal: Running 'onlyS' calculation for '-x' direction...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *_1.alloc *_1.BONDS* *_1.KP *_1.STRUCT_* INPUT_TMP.*

# 2. '+x' move.
echo -e "vibranal: Changing the system label to '${slabel}_2'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_2/I" -i ${FCfdf}

echo -e "\nvibranal: Shifting half of the atoms to '+x' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ} '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1 + var, $2, $3, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies" >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibranal: Running 'onlyS' calculation for '+x' direction...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *_2.alloc *_2.BONDS* *_2.KP *_2.STRUCT_* INPUT_TMP.*

# 3. '-y' move.
echo -e "vibranal: Changing the system label to '${slabel}_3'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_3/I" -i ${FCfdf}

echo -e "\nvibranal: Shifting half of the atoms to '-y' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ} '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2 - var, $3, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies" >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibranal: Running 'onlyS' calculation for '-y' direction...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *_3.alloc *_3.BONDS* *_3.KP *_3.STRUCT_* INPUT_TMP.*

# 4. '+y' move.
echo -e "vibranal: Changing the system label to '${slabel}_4'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_4/I" -i ${FCfdf}

echo -e "\nvibranal: Shifting half of the atoms to '+y' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ} '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2 + var, $3, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies" >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibranal: Running 'onlyS' calculation for '+y' direction...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *_4.alloc *_4.BONDS* *_4.KP *_4.STRUCT_* INPUT_TMP.*

# 5. '-z' move.
echo -e "vibranal: Changing the system label to '${slabel}_5'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_5/I" -i ${FCfdf}

echo -e "\nvibranal: Shifting half of the atoms to '-z' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ} '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2, $3 - var, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies" >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibranal: Running 'onlyS' calculation for '-z' direction...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *_5.alloc *_5.BONDS* *_5.KP *_5.STRUCT_* INPUT_TMP.*

# 6. '+z' move.
echo -e "vibranal: Changing the system label to '${slabel}_6'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_6/I" -i ${FCfdf}

echo -e "\nvibranal: Shifting half of the atoms to '+z' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ} '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2, $3 + var, $4}' \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies" >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibranal: Running 'onlyS' calculation for '+z' direction...\n\n"
${1} -n ${cores} -machinefile ${2} ${3} < ${FCfdf}
wait
rm *_6.alloc *_6.BONDS* *_6.KP *_6.STRUCT_* INPUT_TMP.*

# Finishing.
rm ${FCdir}pbCoord.tmp
rm ${FCdir}pbCoord.fdf
rm ${FCfdf}
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "vibranal: End of run: "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nvibranal: Run time: ${tempo} min\n"

exit 0
