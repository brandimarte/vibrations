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
#   File: buildInput.sh                                   #
#                                                         #
#   Versions: 1 - 10/10/2012                              #
#             2 - 08/01/2013                              #
#                                                         #
#  *****************************************************  #
#  This script collects required informations from the    #
#  input 'fdf' file of a previous force constants (FC)    #
#  run and puts at 'inputFC.in' file.                     #
#  *****************************************************  #

# FC directory
check=`echo "${1}" | sed 's:.*\(.$\):\1:'`
if [ "${check}" == "/" ]
then
    FCdir=${1}
else
    FCdir=${1}/
fi

# Checks if FC directory exists and is readable.
echo -n " Checking FC data path... "
if [ ! -r ${FCdir} ]
then
    echo -e "ERROR: the directory \"${FCdir}\" doesn't exist or is not accessible!\n"
    exit -1
else
    echo -e "ok!\n"
fi

# Gets '.fdf' files.
if ls ${FCdir}*.fdf > /dev/null 2>&1
then
    > ${FCdir}FCfdf.tmp
    FCfdf=`ls ${FCdir}*.fdf`
    echo -e " FC fdf input files found:"
    for fdf in ${FCfdf}
    do
	echo -e "\t\t\t\t${fdf}"
        # Copies to 'FCfdf.tmp' without comented lines.
	sed '/^[#]\+/d' ${fdf} >> ${FCdir}FCfdf.tmp
    done
else
    echo -e " ERROR: couldn't find any FC fdf input file at \"${FCdir}\"!\n"
    exit -1
fi
FCfdf=${FCdir}FCfdf.tmp

# Puts all required data at 'inputFC.in'.
echo -e ""
echo -n " Writing required data at '${FCdir}inputFC.in'... "
> ${FCdir}inputFC.in

check=`grep -i "SYSTEMLABEL" ${FCfdf} | awk '{print $2}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the 'SystemLabel' at FC fdf input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

check=`grep -i "NUMBEROFATOMS" ${FCfdf} | awk '{print $2}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the 'NumberOfAtoms' at FC fdf input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

nSpecies=`grep -i "NUMBEROFSPECIES" ${FCfdf} | awk '{print $2}'`
if [ "${nSpecies}" == "" ]
then
    echo -e "ERROR: can't find the 'NumberOfSpecies' at FC fdf input file!\n"
    exit -1
else
    echo -e ${nSpecies} >> ${FCdir}inputFC.in
fi

check=`grep -i -A${nSpecies} "CHEMICALSPECIESLABEL" ${FCfdf} \
       | head -n $[ ${nSpecies} + 1 ] | tail -n ${nSpecies} \
       | awk '{print $1, $2, $3}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the block 'ChemicalSpeciesLabel' at FC fdf input file!\n"
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
    echo -e ${TorF}
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
    echo -e "ERROR: can't find the 'MD.FCfirst' at FC fdf input file!\n"
    exit -1
else
    echo -e ${FCfirst} >> ${FCdir}inputFC.in
fi

FClast=`grep -i "MD.FCLAST" ${FCfdf} | awk '{print $2}'`
if [ "${FClast}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FClast' at FC fdf input file!\n"
    exit -1
else
    echo -e ${FClast} >> ${FCdir}inputFC.in
fi

check=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $2}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FCdispl' value at FC fdf input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

check=`grep -i "MD.FCDISPL" ${FCfdf} | awk '{print $3}' | tr '[:upper:]' '[:lower:]'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FCdispl' unit at FC fdf input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

dynAtoms=$[ ${FClast} - ${FCfirst} + 1 ]
check=`grep -i -A${FClast} "ATOMICCOORDINATESANDATOMICSPECIES" ${FCfdf} \
       | head -n $[ ${FClast} + 1 ] | tail -n ${dynAtoms} | awk '{print $4}'`
if [ "${check}" == "" ]
then
    echo -e "ERROR: can't find the block 'AtomicCoordinatesAndAtomicSpecies' at FC fdf input file!\n"
    exit -1
else
    echo -e ${check} >> ${FCdir}inputFC.in
fi

echo -e "ok!\n"

rm ${FCdir}FCfdf.tmp

exit 0
