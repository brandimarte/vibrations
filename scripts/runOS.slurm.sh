#!/bin/bash

#  *******************************************************************  #
#                  ** PhOnonS ITeratIVE VIBRATIONS **                   #
#                                                                       #
#                           **  Version 2  **                           #
#                                                                       #
#  Written by Pedro Brandimarte (brandimarte@gmail.com)                 #
#                                                                       #
#  Copyright (c), All Rights Reserved                                   #
#                                                                       #
#  This program is free software. You can redistribute it and/or        #
#  modify it under the terms of the GNU General Public License          #
#  (version 3 or later) as published by the Free Software Foundation    #
#  <http://fsf.org/>.                                                   #
#                                                                       #
#  This program is distributed in the hope that it will be useful, but  #
#  WITHOUT ANY WARRANTY, without even the implied warranty of           #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     #
#  General Public License for more details (file 'LICENSE_GPL'          #
#  distributed along with this program or at                            #
#  <http://www.gnu.org/licenses/gpl.html>).                             #
#  *******************************************************************  #
#                            runOS.slurm.sh                             #
#  *******************************************************************  #
#  This script receives an input 'fdf' file for a force constants (FC)  #
#  run, changes it by removing the flag 'FCwriteHS' and including the   #
#  flag 'FConlyS' and changing the geometry and number of atoms.        #
#                                                                       #
#  Input:  ${1} :  FC calculation main directory                        #
#          ${2} :  FC input file                                        #
#          ${3} :  command to run parallel jobs (e.g. 'srun')           #
#          ${4} :  total number of cores                                #
#          ${5} :  I-SMEAGOL executable with path                       #
#                                                                       #
#  Use: ./runOS.slurm.sh [FC calculation directory] [FC input file]     #
#       [parallel run command] [number of cores]                        #
#       [I-SMEAGOL executable]                                          #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    October 2012                                    #
#                       February 2014                                   #
#                       June 2015                                       #
#  *******************************************************************  #

# Prints the header.
linh="************************************"
echo ""
echo "   ${linh}${linh}"
echo ""
echo "                     ** PhOnonS ITeratIVE VIBRATIONS **"
echo ""
echo "                             **  Version 2  **"
echo ""
echo -n "                       "
date
echo ""
echo "      Written by Pedro Brandimarte (brandimarte@gmail.com)"
echo ""
echo "      Copyright (c), All Rights Reserved"
echo ""
echo "      This program is free software. You can redistribute it"     \
    "and/or"
echo "      modify it under the terms of the GNU General Public"        \
    "License"
echo "      (version 3 or later) as published by the Free Software"     \
    "Foundation"
echo "      <http://fsf.org/>. See the GNU General Public License for"  \
    "details."
echo ""
echo "   ${linh}${linh}"

# Checks if the number of arguments is correct.
if [ ${#} != 5 ]
then
    echo -e "\nvibrations: ERROR: wrong number of arguments!\n"
    echo -e "vibrations: Use: ./runOS.slurm.sh"                         \
	"[FC calculation directory] [FC input file] \ "
    echo -e "               [parallel run command] [number of cores]"   \
	"[I-SMEAGOL executable]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "vibrations: Start of runs: "
date
begin=$(date +%s%N) # initial time with nanoseconds accuracy

# FC and work directories.
check=`echo "${1}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    FCdir=${1}
else
    FCdir=${1}/
fi

# Checks if the files and FC folder exists and are accessible.
echo -e ""
echo -n "vibrations: Checking input... "
if [ ! -r ${FCdir} ]
then
    echo -e "\nvibrations: ERROR: the directory \"${FCdir}\" doesn't"   \
	"exist or is not accessible!\n"
    exit -1
elif [ ! -r ${2} ]
then
    echo -e "\nvibrations: ERROR: the file \"${2}\" doesn't exist or"   \
	"is not accessible!\n"
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
> ${FCdir}FCfdf.tmp
for fdf in ${2}
do
    # Copies to 'FCfdf.tmp' without comented lines.
    sed '/^[#]/d' ${fdf} >> ${FCdir}FCfdf.tmp
done
cat ${FCdir}FCfdf.tmp > ${FCdir}pbFCrun.fdf
rm ${FCdir}FCfdf.tmp
FCfdf=${FCdir}pbFCrun.fdf

# Sets FC 'fdf' input file with the flag 'FCwriteHS .false.'.
echo -e "\nvibrations: Including the option 'FCwriteHS .false.' at"     \
    "'pbFCrun.fdf' file.\n"
check=`grep -i "FCWRITEHS" ${FCfdf}`
if [ "${check}" == "" ]
then
    echo -e "\nFCwriteHS            .false." >> ${FCfdf}
else
    sed "s/${check}/FCwriteHS            .false.\n/I" -i ${FCfdf}
fi

# Sets FC 'fdf' input file with the flag 'FConlyS .true.'.
echo -e "vibrations: Including the option 'FConlyS .true.' at"          \
    "'pbFCrun.fdf' file.\n"
check=`grep -i "FCONLYS" ${FCfdf}`
if [ "${check}" == "" ]
then
    echo -e "\nFConlyS              .true." >> ${FCfdf}
else
    sed "s/${check}/FConlyS              .true.\n/I" -i ${FCfdf}
fi

# Gets the atoms displacement.
echo -n "vibrations: Atoms displacement: "
displ=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $2}'`
if [ "${displ}" == "" ]
then
    echo -e "\nvibrations: ERROR: can't find the 'MD.FCdispl' value at" \
	"FC input file!\n"
    exit -1
fi

# Gets the displacement unit (Bohr or Ang).
dunit=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $3}'                 \
       | tr '[:upper:]' '[:lower:]'`
if [ "${dunit}" == "" ]
then
    echo -e "\nvibrations: ERROR: can't find the 'MD.FCdispl' unit at"  \
	"FC input file!\n"
    exit -1
fi
echo -e "${displ} ${dunit}\n"

# Gets the atomic coordinates unit (Bohr or Ang).
aunit=`grep -i "ATOMICCOORDINATESFORMAT" ${FCfdf} | awk '{print $2}'    \
       | tr '[:upper:]' '[:lower:]'`
if [ "${aunit}" == "" ]
then
    echo -e "\nvibrations: ERROR: can't find the"                       \
	"'AtomicCoordinatesFormat' at FC fdf input file!\n"
    exit -1
elif [ "${aunit}" == "notscaledcartesianang" ]
then
    aunit="ang"
elif [ "${aunit}" == "notscaledcartesianbohr" ]
then
    aunit="bohr"
elif [ "${aunit}" != "bohr" ]
then
    if [ "${aunit}" != "ang" ]
    then
	echo -e "\nvibrations: ERROR: unrecognized atomic coordinates"  \
	    "format \"${aunit}\" at"
	echo -e "                 FC input file! It must be \"Ang\" or" \
	    " \"Bohr\"!\n"
	exit -1
    fi
fi

# Converts the displacement unit if necessary.
echo -e "vibrations: Atomic coordinates unit: ${aunit}\n"
if [ "${dunit}" != "${aunit}" ]
then
    echo -e "vibrations: Converting the displacement to ${aunit}.\n"
    if [ "${dunit}" == "bohr" ]
    then
	displ=`echo "scale = 15; ${displ} * 0.52917721092" | bc`
    elif [ "${dunit}" == "ang" ]
    then
	displ=`echo "scale = 15; ${displ} / 0.52917721092" | bc`
    else
	echo -e "\nvibrations: ERROR: unrecognized displacement unit"   \
	    "\"${dunit}\" at"
	echo -e "                 FC input file! It must be \"Ang\" or" \
	    " \"Bohr\"!\n"
	exit -1
    fi
fi

# Gets the number of atoms.
natoms=`grep -i "NUMBEROFATOMS" ${FCfdf} | awk '{print $2}'`
if [ "${natoms}" == "" ]
then
    echo -e "\nvibrations: ERROR: can't find the 'NumberOfAtoms' at FC" \
	"input file!\n"
    exit -1
fi
echo -e "vibrations: Actual number of atoms: ${natoms}\n"

# Doubles the number of atoms.
echo -n "vibrations: Doubling the number of atoms from ${natoms} to "
check=`grep -i -m 1 "NUMBEROFATOMS" ${FCfdf}`
n2atoms=`echo "${natoms} * 2" | bc`
sed "s/${check}/NumberOfAtoms           ${n2atoms}\n/I" -i ${FCfdf}
echo -e "${n2atoms}\n"

# Gets and removes the 'AtomicCoordinatesAndAtomicSpecies' block.
echo -n "vibrations: Reading the block"                                 \
    "'AtomicCoordinatesAndAtomicSpecies'... "
atcs="ATOMICCOORDINATESANDATOMICSPECIES"
acoord=`grep -i -m 1 -A ${natoms} ${atcs} ${FCfdf}| tail -n ${natoms}`
if [ "${acoord}" == "" ]
then
    echo -e "\nvibrations: ERROR: can't find the block"                 \
	"'AtomicCoordinatesAndAtomicSpecies'"
    echo -e "                 at FC input file!\n"
    exit -1
else
    echo -e "ok!\n"

    # Copies 'AtomicCoordinatesAndAtomicSpecies' block to 'pbCoord.tmp'.
    echo -e "vibrations: Copying 'AtomicCoordinatesAndAtomicSpecies'"   \
	"to 'pbCoord.tmp' file.\n"
    grep -i -m 1 -A ${natoms} "ATOMICCOORDINATESANDATOMICSPECIES"       \
	${FCfdf} | tail -n ${natoms} > ${FCdir}pbCoord.tmp

    # Removes 'AtomicCoordinatesAndAtomicSpecies' block.
    echo -e "vibrations: Removing 'AtomicCoordinatesAndAtomicSpecies'"  \
	"from 'pbFCrum.fdf' file.\n"
    first=`nl -b a ${FCfdf} | grep -i -m 1 ${atcs} | awk '{print $1}'`
    last=`echo "${first} + ${natoms} + 1" | bc`
    sed "${first},${last}d" -i ${FCfdf}

    # Removes '%include' lines.
    echo -e "vibrations: Removing '%include' lines from 'pbFCrum.fdf'.\n"
    check=`grep -i -m 1 "%INCLUDE" ${FCfdf}`
    while [ "${check}" != "" ]
    do
	sed "s/${check}/ /" -i ${FCfdf}	
	check=`grep -i -m 1 "%INCLUDE" ${FCfdf}`
    done
fi

# Add the external coordinate flag.
echo -e "vibrations: Adding '%include pbCoord.fdf' at 'pbFCrun.fdf'.\n"
echo -e "%include                pbCoord.fdf" >> ${FCfdf}

slabel=`grep -i "SYSTEMLABEL" ${FCfdf} | awk '{print $2}'`
if [ "${slabel}" == "" ]
then
    echo -e "ERROR: can't find the 'SystemLabel' at FC input file!\n"
    exit -1
fi

# --- Only S runs: ---

echo -e "vibrations: Including atomic coordinates with at"              \
    "'pbCoord.fdf'.\n"
echo -e "%block AtomicCoordinatesAndAtomicSpecies" > ${FCdir}pbCoord.fdf
awk '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2, $3, $4}'          \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
submit="${3} -n ${4} ${5}"

# 1. '-x' move.
echo -e "vibrations: Changing the system label to '${slabel}_1'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_1/I" -i ${FCfdf}

echo -e "vibrations: Shifting half of the atoms to '-x' direction.\n"
awk -v var=${displ}                                                     \
    '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1 - var, $2, $3, $4}'    \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies"                   \
    >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibrations: Running 'onlyS' calculation for '-x'"              \
    "direction...\n\n"
${submit} < ${FCfdf}
wait
rm *_1.alloc *_1.BONDS* *_1.KP *_1.STRUCT_* INPUT_TMP.*

# 2. '+x' move.
echo -e "vibrations: Changing the system label to '${slabel}_2'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_2/I" -i ${FCfdf}

echo -e "\nvibrations: Shifting half of the atoms to '+x' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ} \
    '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1 + var, $2, $3, $4}'    \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies"                   \
    >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibrations: Running 'onlyS' calculation for '+x'"              \
    "direction...\n\n"
${submit} < ${FCfdf}
wait
rm *_2.alloc *_2.BONDS* *_2.KP *_2.STRUCT_* INPUT_TMP.*

# 3. '-y' move.
echo -e "vibrations: Changing the system label to '${slabel}_3'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_3/I" -i ${FCfdf}

echo -e "\nvibrations: Shifting half of the atoms to '-y' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ}                                                     \
    '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2 - var, $3, $4}'    \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies"                   \
    >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibrations: Running 'onlyS' calculation for '-y'"              \
    "direction...\n\n"
${submit} < ${FCfdf}
wait
rm *_3.alloc *_3.BONDS* *_3.KP *_3.STRUCT_* INPUT_TMP.*

# 4. '+y' move.
echo -e "vibrations: Changing the system label to '${slabel}_4'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_4/I" -i ${FCfdf}

echo -e "\nvibrations: Shifting half of the atoms to '+y' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ}                                                     \
    '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2 + var, $3, $4}'    \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies"                   \
    >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibrations: Running 'onlyS' calculation for '+y'"              \
    "direction...\n\n"
${submit} < ${FCfdf}
wait
rm *_4.alloc *_4.BONDS* *_4.KP *_4.STRUCT_* INPUT_TMP.*

# 5. '-z' move.
echo -e "vibrations: Changing the system label to '${slabel}_5'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_5/I" -i ${FCfdf}

echo -e "\nvibrations: Shifting half of the atoms to '-z' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ}                                                     \
    '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2, $3 - var, $4}'    \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies"                   \
    >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibrations: Running 'onlyS' calculation for '-z'"              \
    "direction...\n\n"
${submit} < ${FCfdf}
wait
rm *_5.alloc *_5.BONDS* *_5.KP *_5.STRUCT_* INPUT_TMP.*

# 6. '+z' move.
echo -e "vibrations: Changing the system label to '${slabel}_6'.\n"
check=`grep -i "SYSTEMLABEL" ${FCfdf}`
sed "s/${check}/SystemLabel             ${slabel}_6/I" -i ${FCfdf}

echo -e "\nvibrations: Shifting half of the atoms to '+z' direction.\n"
first=$[ ${natoms} + 2 ]
last=$[ 2 * ${natoms} + 2 ]
sed "${first},${last}d" -i ${FCdir}pbCoord.fdf
awk -v var=${displ}                                                     \
    '{printf "\t% .10f\t% .10f\t% .10f%3d\n", $1, $2, $3 + var, $4}'    \
    ${FCdir}pbCoord.tmp >> ${FCdir}pbCoord.fdf
echo -e "%endblock AtomicCoordinatesAndAtomicSpecies"                   \
    >> ${FCdir}pbCoord.fdf

# Submits the job.
echo -e "vibrations: Running 'onlyS' calculation for '+z'"              \
    "direction...\n\n"
${submit} < ${FCfdf}
wait
rm *_6.alloc *_6.BONDS* *_6.KP *_6.STRUCT_* INPUT_TMP.*

# Finishing.
rm ${FCdir}pbCoord.tmp
rm ${FCdir}pbCoord.fdf
rm ${FCfdf}
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "vibrations: End of run: "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nvibrations: Run time: ${tempo} min\n"

exit 0
