#!/bin/bash

#  *******************************************************************  #
#                  ** PhOnonS ITeratIVE VIBRATIONS **                   #
#                                                                       #
#                           **  Version 2  **                           #
#                                                                       #
#   By: Pedro Brandimarte (brandimarte@gmail.com) and                   #
#       Alexandre Reily Rocha (reilya@ift.unesp.br)                     #
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
#                              setInput.sh                              #
#  *******************************************************************  #
#  This script collects required informations from the input 'fdf'      #
#  file of a previous force constants (FC) run and copies them at       #
#  'inputFC.in' file. Concatenates the '.FC' files at 'FC*' folders     #
#  and copies '.gHS' files at 'FC*' folders, renaming them              #
#  accordingly.                                                         #
#                                                                       #
#  Input:  ${1} :  FC calculation main directory                        #
#          ${2} :  FC input file                                        #
#          ${3} :  calculation type (full or onlyPh)                    #
#          ${4} :  splitFC (optional)                                   #
#                                                                       #
#  Use: ./setInputs.sh [FC calculation directory] [FC input file]       #
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

# FC directory
check=`echo "${1}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    FCdir=${1}
else
    FCdir=${1}/
fi

# Checks if the files and FC folder exists and are accessible.
echo -e ""
echo -n " Checking input... "
if [ ! -r ${FCdir} ]
then
    echo -e "\nvibrations: ERROR: the directory \"${FCdir}\" doesn't"   \
	"exist or is not accessible!\n"
    exit -1
elif [ ! -r ${FCdir}${2} ]
then
    echo -e "\nvibrations: ERROR: the file \"${2}\" doesn't exist or"   \
	"is not accessible!\n"
    exit -1
else
    echo -e "ok!\n"
fi

# Gets '.fdf' files.
> ${FCdir}FCfdf.tmp
for fdf in ${FCdir}${2}
do
    # Copies to 'FCfdf.tmp' without comented lines.
    sed '/^[#]/d' ${fdf} >> ${FCdir}FCfdf.tmp
done
FCfdf=${FCdir}FCfdf.tmp

# Puts all required data at 'inputFC.in'.
echo -e " Writing required data at '${FCdir}inputFC.in'."
> ${FCdir}inputFC.in

slabel=`grep -i "SYSTEMLABEL" ${FCfdf} | awk '{print $2}'`
if [ "${slabel}" == "" ]
then
    echo -e "ERROR: can't find the 'SystemLabel' at FC input file!\n"
    exit -1
else
    echo -e ${slabel} >> ${FCdir}inputFC.in
fi

natoms=`grep -i "NUMBEROFATOMS" ${FCfdf} | awk '{print $2}'`
if [ "${natoms}" == "" ]
then
    echo -e "ERROR: can't find the 'NumberOfAtoms' at FC input file!\n"
    exit -1
else
    echo -e ${natoms} >> ${FCdir}inputFC.in
fi

nSpecies=`grep -i "NUMBEROFSPECIES" ${FCfdf} | awk '{print $2}'`
if [ "${nSpecies}" == "" ]
then
    echo -e "ERROR: can't find the 'NumberOfSpecies' at FC input file!\n"
    exit -1
else
    echo -e ${nSpecies} >> ${FCdir}inputFC.in
fi

check=`grep -i -A${nSpecies} "CHEMICALSPECIESLABEL" ${FCfdf}            \
       | head -n $[ ${nSpecies} + 1 ] | tail -n ${nSpecies}             \
       | awk '{print $1, $2, $3}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the block 'ChemicalSpeciesLabel' at FC"  \
	"input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

check=`grep -i "SPINPOLARIZED" ${FCfdf} | awk '{print $2}'`
if [ "${check}" == "" ]
then
    echo -e "1" >> ${FCdir}inputFC.in
else
    TorF=`echo ${check} | grep -i "TRUE"`
    if [ "${TorF}" == "" ]
    then
	echo -e "1" >> ${FCdir}inputFC.in
    else
	echo -e "2" >> ${FCdir}inputFC.in
    fi
fi

FCfirst=`grep -i "MD.FCFIRST" ${FCfdf} | awk '{print $2}'`
if [ "${FCfirst}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FCfirst' at FC input file!\n"
    exit -1
else
    echo -e ${FCfirst} >> ${FCdir}inputFC.in
fi

FClast=`grep -i "MD.FCLAST" ${FCfdf} | awk '{print $2}'`
if [ "${FClast}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FClast' at FC input file!\n"
    exit -1
else
    echo -e ${FClast} >> ${FCdir}inputFC.in
fi

check=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $2}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FCdispl' value at FC input"      \
	"file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

check=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $3}'                 \
       | tr '[:upper:]' '[:lower:]'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FCdispl' unit at FC input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

dynAtoms=$[ ${FClast} - ${FCfirst} + 1 ]
check=`grep -i -A${FClast} "ATOMICCOORDINATESANDATOMICSPECIES" ${FCfdf} \
       | head -n $[ ${FClast} + 1 ] | tail -n ${dynAtoms}               \
       | awk '{print $4}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the block"                               \
	"'AtomicCoordinatesAndAtomicSpecies' at FC input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

rm ${FCdir}FCfdf.tmp

if [ "${4}" != "" ]
then
    # Concatenate calculated constant forces matrices and Fermi
    # energies, and copy '.gHS' files, renaming accordingly.
    echo "Force constants matrix" > ${FCdir}${slabel}.FC

    if [ "${3}" == "1" ]
    then
	head -n1 ${FCdir}FC${FCfirst}/${slabel}.ef > ${FCdir}${slabel}.ef
        nlines=$(( ${natoms} * 6 ))
	j=0
	for i in `seq ${FCfirst} ${FClast}`
	do
            # Concatenated FC matrix.
	    tail -n${nlines} ${FCdir}FC${i}/${slabel}.FC                \
		>> ${FCdir}${slabel}.FC
	
            # Get Fermi energies.
	    awk 'NR>=2 && NR<=7' ${FCdir}FC${i}/${slabel}.ef            \
		| awk -v j=$j '{printf ("%12d  % .14f\n", NR+j, $2);}'  \
		>> ${FCdir}${slabel}.ef

	    for k in `seq 1 6`
	    do
		intIn=`printf "%.3d" ${k}`
		intOut=`printf "%.3d" $(( ${j} + ${k} ))`
		cp ${FCdir}FC${i}/${slabel}_${intIn}.gHS                \
		    ${FCdir}${slabel}_${intOut}.gHS
	    done

	    j=$(( ${j} + 6 ))
	done

        # Copy one of the '000.gHS' and '.orb' to main FC folder.
	cp ${FCdir}FC${FCfirst}/${slabel}_000.gHS ${FCdir}
	cp ${FCdir}FC${FCfirst}/${slabel}.orb ${FCdir}

    else # onlyPh calculation
	for i in `seq ${FCfirst} ${FClast}`
	do
            # Concatenated FC matrix.
	    tail -n${nlines} ${FCdir}FC${i}/${slabel}.FC                \
		>> ${FCdir}${slabel}.FC
	done

    fi
fi

exit 0
