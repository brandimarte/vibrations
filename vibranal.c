/**  *****************************************************  **/
/**             ** Phonon Vibration Analysis **             **/
/**                                                         **/
/**                    **  Version 2  **                    **/
/**                                                         **/
/**                         IF/USP                          **/
/**                                                         **/
/**   Advisor: Prof. Dr. Alexandre Reily Rocha              **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**   File: vibranal.c                                      **/
/**                                                         **/
/**   Versions: 1 - 10/10/2012                              **/
/**             2 - 08/01/2013                              **/
/**                                                         **/
/**  *****************************************************  **/
/**  This is a (client) program that reads the output data  **/
/**  from a siesta force constants (FC) calculation runned  **/
/**  with the flag 'PB.FCwriteHS .true.' at the fdf input   **/
/**  file. Computes the phonon energies and modes, the      **/
/**  hamiltonian derivatives with respect to the dynamic    **/
/**  atoms displacements and computes the electron-phonon   **/
/**  coupling matrix.                                       **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Utils.h"
#include "Phonon.h"

int main (int nargs, char *arg[])
{
   int nDynTot, nDynOrb;
   double *EigVec, *EigVal, *Meph;
   double time;
   clock_t inicial, final;

   /* Checks if the input were typed correctly. */
   if (nargs != 2) {
      fprintf (stderr, "\n\n Wrong number of arguments!\n");
      fprintf (stderr, "\n Use: vibranal"); /* arg[0] */
      fprintf (stderr, " [FC directory]\n\n"); /* arg[1] */
      fprintf (stderr, " Example: vibranal ~/MySystem/FCcalc\n\n");
      exit (EXIT_FAILURE);
   }

   /* Time is running. */
   inicial = clock();

   /* Writes the header on the screen. */
   PHONheader ();

   /* Reads info from FC fdf input file. */
   PHONreadFCfdf (arg[0], arg[1], &nDynTot, &nDynOrb);

   /* Computes phonon frequencies. */
   EigVec = UTILdoubleVector (nDynTot * nDynTot);
   EigVal = UTILdoubleVector (nDynTot);
   PHONfreq (EigVec, EigVal);

   /* Computes electron-phonon coupling matrices. */
   Meph = UTILdoubleVector (nDynOrb * nDynOrb * nDynTot);
   PHONephCoupling (EigVec, EigVal, Meph);

   /* Frees memory. */
   free (EigVec);
   free (EigVal);
   free (Meph);

   /* Calculates the execution time. */
   final = clock();
   time = (double)(final - inicial) / CLOCKS_PER_SEC;
   printf ("\n This calculation took %.2f seconds.\n\n", time);

   return 0;
   
} /* main */

