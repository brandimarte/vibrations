#!/bin/bash

#  *******************************************************************  #
#                       ** Vibration Analysis **                        #
#                                                                       #
#                           **  Version 2  **                           #
#                                                                       #
#                                IF/USP                                 #
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
#                               runFCs.sh                               #
#  *******************************************************************  #
#  This script receives an input 'fdf' file for a force constants (FC)  #
#  run, includes the flag 'FCwriteHS' and submits the calculation.      #
#                                                                       #
#  Input:  ${1} :  FC calculation main directory                        #
#          ${2} :  FC input file                                        #
#          ${3} :  PBS submit script                                    #
#                                                                       #
#  Use: ./runFCs.sh [FC calculation directory] [FC input file]          #
#       [PBS submit script]                                             #
#                                                                       #
#  Written by Pedro Brandimarte, Feb 2014.                              #
#  Instituto de Fisica                                                  #
#  Universidade de Sao Paulo                                            #
#  e-mail: brandimarte@gmail.com                                        #
#  ***************************** HISTORY *****************************  #
#  Original version:    October 2012                                    #
#                       February 2014                                   #
#  *******************************************************************  #

# Prints the header.
linh="************************************"
echo ""
echo "   ${linh}${linh}"
echo ""
echo "                       **  VIBRATIONAL ANALYSIS  **"
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
if [ ${#} != 3 ]
then
    echo -e "\nvibranal: ERROR: wrong number of arguments!\n"
    echo -e "vibranal: Use: ./runFCs.sh [FC calculation directory]"     \
	"[FC input file] \ "
    echo -e "               [PBS submit script]\n"
    exit -1
fi

# Time is running.
echo -e ""
echo -n "vibranal: Start of runs: "
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
echo -n "vibranal: Checking input... "
if [ ! -r ${FCdir} ]
then
    echo -e "\nvibranal: ERROR: the directory \"${FCdir}\" doesn't"     \
	"exist or is not accessible!\n"
    exit -1
elif [ ! -r ${2} ]
then
    echo -e "\nvibranal: ERROR: the file \"${2}\" doesn't exist or is"  \
	"not accessible!\n"
    exit -1
elif [ ! -r ${3} ]
then
    echo -e "\nvibranal: ERROR: the file \"${3}\" doesn't exist or is"  \
	"not accessible!\n"
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
> ${FCdir}FCfdf.tmp
for fdf in ${2}
do
    # Copies to 'FCfdf.tmp' without comented lines.
    sed '/^[#]/d' ${fdf} >> ${FCdir}FCfdf.tmp
done
cat ${FCdir}FCfdf.tmp > ${FCdir}pbFCrun.fdf
rm ${FCdir}FCfdf.tmp
FCfdf=${FCdir}pbFCrun.fdf

# Sets FC 'fdf' input file with the flag 'FCwriteHS .true.'.
echo -e "vibranal: Including the option 'FCwriteHS .true.' at"          \
    "'pbFCrun.fdf' file.\n"
check=`grep -i "FCWRITEHS" ${FCfdf}`
if [ "${check}" == "" ]
then
    echo -e "\nFCwriteHS            .true." >> ${FCfdf}
else
    sed "s/${check}/FCwriteHS            .true.\n/I" -i ${FCfdf}
fi

# Sets FC 'fdf' input file with the flag 'FConlyS .false.'.
echo -e "vibranal: Including the option 'FConlyS .false.' at"           \
    "'pbFCrun.fdf' file.\n"
check=`grep -i "FCONLYS" ${FCfdf}`
if [ "${check}" == "" ]
then
    echo -e "\nFConlyS              .false." >> ${FCfdf}
else
    sed "s/${check}/FConlyS              .false.\n/I" -i ${FCfdf}
fi

PBSlabel=`grep "#PBS[[:blank:]]*-N" ${3} | awk '{print $3}'`
PBSname="#PBS[[:blank:]]*-N[[:blank:]]*${PBSlabel}"
if [ "${PBSlabel}" == "" ]
then
    echo -e "ERROR: can't find the '#PBS -N' at PBS submit script!\n"
    exit -1
fi

slabel=`grep -i "SYSTEMLABEL" ${FCfdf} | awk '{print $2}'`
if [ "${slabel}" == "" ]
then
    echo -e "ERROR: can't find the 'SystemLabel' at FC input file!\n"
    exit -1
fi

MDFCfirst=`grep -i "MD.FCFIRST" ${FCfdf} | awk '{print $1}'`
FCfirst=`grep -i "MD.FCFIRST" ${FCfdf} | awk '{print $2}'`
if [ "${FCfirst}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FCfirst' at FC input file!\n"
    exit -1
fi

MDFClast=`grep -i "MD.FCLAST" ${FCfdf} | awk '{print $1}'`
FClast=`grep -i "MD.FCLAST" ${FCfdf} | awk '{print $2}'`
if [ "${FClast}" == "" ]
then
    echo -e "ERROR: can't find the 'MD.FClast' at FC input file!\n"
    exit -1
fi

# Get the pseudo-potential file extention.
if ls ${FCdir}*.psf > /dev/null 2>&1
then
    pseudo=".psf"
elif ls ${FCdir}*.vps > /dev/null 2>&1
    pseudo=".vps"
else
    echo -e " ERROR: couldn't find any pseudo-potential file (psf or"   \
	"vps) at \"${FCdir}\"!\n"
    exit -1
fi

for i in `seq ${FCfirst} ${FClast}`
do
    # Create a FC folder for each atom.
    mkdir ${FCdir}FC${i}

    # Copy needed files to these folders.
    cp ${FCdir}*${pseudo} ${FCdir}FC${i}/
    cp ${FCdir}${slabel}.DM ${FCdir}FC${i}/

    # Copy input file with only one dynamic atom.
    sed "s/${MDFCfirst}[[:blank:]]*${FCfirst}/${MDFCfirst}      ${i}/;
         s/${MDFClast}[[:blank:]]*${FClast}/${MDFClast}      ${i}/"     \
	${FCfdf} > ${FCdir}FC${i}/${2}
    
    # Copy PBS submite file and change name.
    sed "s/$PBSname/#PBS -N ${PBSlabel}-${i}/" ${3} > ${FCdir}FC${i}/${3}

    # Submite job.
    cd ${FCdir}FC${i}
    qsub ${3}
    cd ../
    
done

# Finishing.
rm ${FCfdf}
end=$(date +%s%N) # final time with nanoseconds accuracy
echo -e ""
echo -n "vibranal: End of run: "
date
tempo=`echo "scale = 10; (${end} - ${begin}) / 60000000000" | bc`
echo -e "\nvibranal: Run time: ${tempo} min\n"

exit 0
