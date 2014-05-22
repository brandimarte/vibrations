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
/**   File: main.c                                          **/
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
   int nDynTot, nDynOrb, spinPol, calcType;
   double *EigVec, *EigVal, *Meph;
   double time;
   clock_t inicial, final;

   /* Checks if the input were typed correctly. */
   if (nargs < 4 || nargs > 5) {
      fprintf (stderr, "\n\n Wrong number of arguments!\n");
      fprintf (stderr, "\n Use: vibranal"); /* arg[0] */
      fprintf (stderr, " [FC directory]"); /* arg[1] */
      fprintf (stderr, " [FC input file]"); /* arg[2] */
      fprintf (stderr, " [calculation type]"); /* arg[3] */
      fprintf (stderr, " [splitFC]\n\n"); /* arg[4] */
      fprintf (stderr,
	       " Examples : vibranal ~/MySystem/FCdir runFC.in full\n");
      fprintf (stderr,
	       "            vibranal ~/MySystem/FCdir runFC.in full splitFC\n");
      fprintf (stderr,
	       "            vibranal ~/MySystem/FCdir runFC.in onlyPh\n");
      fprintf (stderr,
	       "            vibranal ~/MySystem/FCdir runFC.in onlyPh splitFC\n\n");
      exit (EXIT_FAILURE);
   }

   /* Time is running. */
   inicial = clock();

   /* Writes the header on the screen. */
   PHONheader ();

   /* Checks calculation option (full or onlyPh). */
   if (strcmp (arg[3], "full") == 0)
      calcType = 1;
   else {
      if (strcmp (arg[3], "onlyPh") == 0)
	 calcType = 2;
      else {
	 fprintf (stderr, "\n\n Wrong calculation type option!\n");
	 fprintf (stderr, "\n It must be 'full' or 'onlyPh'");
	 exit (EXIT_FAILURE);
      }
   }

   /* Reads info from FC fdf input file. */
   if (nargs == 4)
      PHONreadFCfdf (arg[0], arg[1], arg[2], calcType, " ",
		     &nDynTot, &nDynOrb, &spinPol);
   else
      PHONreadFCfdf (arg[0], arg[1], arg[2], calcType, arg[4],
		     &nDynTot, &nDynOrb, &spinPol);

   /* Computes phonon frequencies. */
   EigVec = UTILdoubleVector (nDynTot * nDynTot);
   EigVal = UTILdoubleVector (nDynTot);
   PHONfreq (EigVec, EigVal);

   /* Writes a 'xyz' file for each computed phonon mode. */
   PHONjmolVib (EigVec);

   /* At "onlyPh" calculations it is finished. */
   if (calcType == 2) {

      /* Frees memory. */
      free (EigVec);
      free (EigVal);

      /* Calculates the execution time. */
      final = clock();
      time = (double)(final - inicial) / CLOCKS_PER_SEC;
      printf ("\n This calculation took %.2f seconds.\n\n", time);

      return 0;

   }

   /* Computes electron-phonon coupling matrices. */
   Meph = UTILdoubleVector (nDynOrb * nDynOrb * nDynTot * spinPol);
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

