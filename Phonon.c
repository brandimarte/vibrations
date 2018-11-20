/**  *****************************************************  **/
/**             ** Phonon Vibration Analysis **             **/
/**                                                         **/
/**                    **  Version 2  **                    **/
/**                                                         **/
/**   By: Pedro Brandimarte (brandimarte@gmail.com) and     **/
/**       Alexandre Reily Rocha (reilya@ift.unesp.br)       **/
/**                                                         **/
/**  *****************************************************  **/
/**  Implementation of phonon analysis:                     **/
/**   - reads the output data from a siesta force           **/
/**   constants (FC) calculation runned with the script     **/
/**   'runFC.sh';                                           **/
/**   - builds the force constants matrix and computes the  **/
/**   phonon frequencies and modes by finite differences;   **/
/**   - computes the Hamiltonian and overlap matrices       **/
/**   derivatives with respect to the dynamic atoms         **/
/**   displacements and apply suitable corrections;         **/
/**   - calculates the electron-phonon coupling matrix.     **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Extern.h"
#include "Check.h"
#include "Utils.h"
#include "Phonon.h"

/* Structure for atomic number and mass. */
typedef struct ATNMASS nmass;
struct ATNMASS {
   int Z; /* atomic number */
   double A; /* atomic mass */
};

/* Structure for chemical element info. */
typedef struct CHEMIC element;
struct CHEMIC {
   int id; /* number of the atom at system structure */
   nmass atom; /* atomic number and mass */
};

/* Structure for 'xyz' coordinates. */
typedef struct ATCOORD xyzcoord;
struct ATCOORD {
   char name[2]; /* element specie */
   double x; /* x coordinate */
   double y; /* y coordinate */
   double z; /* z coordinate */
};

static char *workDir; /* work directory */
static char *FCdir; /* FC directory */
static char sysLabel[30]; /* system label */
static int nAtoms; /* number of atoms */
static int nspin; /* spin polarization */
static int FCfirst; /* first dynamic atom */
static int FClast; /* last dynamic atom */
static int nDyn; /* number of dynamic atoms */
static int no_u; /* number of basis orbitals from unit cell */
static int *orbIdx; /* first orbital index of each atom. */
static double FCdispl; /* atoms displacement */
static double *ef; /* Fermi energy values */
static element *dynAtoms; /* dynamic atoms chemical info */
static nmass periodicTable[95] = { /* {Z,A} - from SIESTA 3.1 */
   { 0,  0.00},{ 1,  1.01},{ 2,  4.00},{ 3,  6.94},{ 4,  9.01},
   { 5, 10.81},{ 6, 12.01},{ 7, 14.01},{ 8, 16.00},{ 9, 19.00},
   {10, 20.18},{11, 22.99},{12, 24.31},{13, 26.98},{14, 28.09},
   {15, 30.97},{16, 32.07},{17, 35.45},{18, 39.95},{19, 39.10},
   {20, 40.08},{21, 44.96},{22, 47.88},{23, 50.94},{24, 52.00},
   {25, 54.94},{26, 55.85},{27, 58.93},{28, 58.69},{29, 63.55},
   {30, 65.39},{31, 69.72},{32, 72.61},{33, 74.92},{34, 78.96},
   {35, 79.90},{36, 83.80},{37, 85.47},{38, 87.62},{39, 88.91},
   {40, 91.22},{41, 92.91},{42, 95.94},{43, 98.91},{44,101.07},
   {45,102.91},{46,106.42},{47,107.87},{48,112.41},{49,114.82},
   {50,118.71},{51,121.75},{52,127.60},{53,126.90},{54,131.29},
   {55,132.91},{56,137.33},{57,138.91},{58,140.12},{59,140.91},
   {60,144.24},{61,146.92},{62,150.36},{63,151.97},{64,157.25},
   {65,158.93},{66,162.50},{67,164.93},{68,167.26},{69,168.93},
   {70,173.04},{71,174.97},{72,178.49},{73,180.95},{74,183.85},
   {75,186.21},{76,190.20},{77,192.22},{78,195.08},{79,196.97},
   {80,200.59},{81,204.38},{82,207.20},{83,208.98},{84,208.98},
   {85,209.99},{86,222.02},{87,223.02},{88,226.03},{89,227.03},
   {90,232.04},{91,231.04},{92,238.03},{93,237.05},{94,244.06}
};


/* ********************************************************* */
/* Reads the file 'inputFC.in' and assigns static global     */
/* variables.                                                */
static void assignGlobalVar ()
{
   register int i, j, len;
   int nSpecies;
   element *species;
   char FCdisplUnit[5]; /* unit of the displacement */
   char **name; /* element name */
   char *inputFile;
   FILE *FCinfo;

   /* Sets the file name 'inputFC.in' with work directory path. */
   len = strlen (workDir);
   inputFile = CHECKmalloc ((len + 11) * sizeof (char));
   sprintf (inputFile, "%sinputFC.in", workDir);

   /* Opens the file with FC run informations. */
   FCinfo = CHECKfopen (inputFile, "r");

   /* Reads and assigns global variables. */
   CHECKfscanf (fscanf (FCinfo, "%s", sysLabel), inputFile);
   printf ("    System label:\t\t\t%s\n", sysLabel);

   CHECKfscanf (fscanf (FCinfo, "%d", &nAtoms), inputFile);
   printf ("    Number of atoms:\t\t\t%d\n", nAtoms);

   /* Different species from the system. */
   CHECKfscanf (fscanf (FCinfo, "%d", &nSpecies), inputFile);
   printf ("    Different system species:\t\t%d\n", nSpecies);
   species = CHECKmalloc (nSpecies * sizeof (element));
   name = CHECKmalloc (nSpecies * sizeof (*name));
   for (i = 0; i < nSpecies; i++) {
      CHECKfscanf (fscanf (FCinfo, "%d", &species[i].id), inputFile);
      CHECKfscanf (fscanf (FCinfo, "%d", &species[i].atom.Z), inputFile);
      name[i] = CHECKmalloc (10 * sizeof (char));
      CHECKfscanf (fscanf (FCinfo, "%s", name[i]), inputFile);
      printf ("\t\t\t\t\t %d %d %s\n", species[i].id,
	      species[i].atom.Z, name[i]);
   }

   CHECKfscanf (fscanf (FCinfo, "%d", &nspin), inputFile);
   printf ("    Spin polarization (1 or 2):\t\t%d\n", nspin);

   CHECKfscanf (fscanf (FCinfo, "%d", &FCfirst), inputFile);
   printf ("    First dynamic atom:\t\t\t%d\n", FCfirst);

   CHECKfscanf (fscanf (FCinfo, "%d", &FClast), inputFile);
   printf ("    Last dynamic atom:\t\t\t%d\n", FClast);

   nDyn = FClast - FCfirst + 1;
   printf ("    Number of dynamic atoms:\t\t%d\n", nDyn);

   CHECKfscanf (fscanf (FCinfo, "%lf", &FCdispl), inputFile);
   CHECKfscanf (fscanf (FCinfo, "%s", FCdisplUnit), inputFile);
   if (strcmp (FCdisplUnit, "bohr") == 0) {
      FCdispl = bohr2ang * FCdispl;
   }
   else if (strcmp (FCdisplUnit, "ang") != 0) {
      fprintf (stderr,
	       "\n ERROR: unrecognized displacement unit \"%s\" at FC\n",
	       FCdisplUnit);
      fprintf (stderr,
	       "       fdf input file! It must be \"Ang\" or \"Bohr\"!");
      fprintf (stderr, "\n\n");
      exit (EXIT_FAILURE);
   }
   printf ("    Atoms displacement:\t\t\t%.5f (ang)\n", FCdispl);

   /* Species of the dynamic atoms. */
   printf ("    Dynamic atoms species:\n");
   dynAtoms = CHECKmalloc (nDyn * sizeof (element));
   for (i = 0; i < nDyn; i++)
      CHECKfscanf (fscanf (FCinfo, "%d", &dynAtoms[i].id), inputFile);
   for (i = 0; i < nDyn; i++)
      for (j = 0; j < nSpecies; j++) {
	 /* Searches the correspondent atom specie. */
	 if (dynAtoms[i].id == species[j].id) {
	    dynAtoms[i].id = FCfirst + i;
	    dynAtoms[i].atom.Z = species[j].atom.Z;
	    dynAtoms[i].atom.A = periodicTable[dynAtoms[i].atom.Z].A;
	    printf ("\t\t\t\t\t%d %d %.2f %s\n", dynAtoms[i].id,
		    dynAtoms[i].atom.Z, dynAtoms[i].atom.A, name[j]);
	    break ;
	 }
      }
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Closes the file 'inputFC.in'. */
   CHECKfclose (fclose (FCinfo), inputFile);

   /* Frees memory. */
   free (inputFile);
   free (species);
   for (i = 0; i < nSpecies; i++)
      free (name[i]);
   free (name);
   
} /* assignGlobalVar */


/* ********************************************************* */
/* Reads at '.ef' file the Fermi energy from the undisplaced */
/* system and the Fermi energy obtained after each           */
/* displacement.                                             */
static void readFermiEnergy ()
{
   register int i, len;
   char *efFile;
   FILE *EF;

   /* Sets the '.ef' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   efFile = CHECKmalloc ((len + 4) * sizeof (char));
   efFile[0] = '\0';
   sprintf (efFile, "%s%s.ef", FCdir, sysLabel);

   /* Opens the '.ef' file. */
   printf ("\n    reading \"%s\" file... ", efFile);
   EF = CHECKfopen (efFile, "r");

   /* Reads the Fermi energies. */
   ef = UTILdoubleVector (6 * nDyn + 1);
   for (i = 0; i < 6 * nDyn + 1; i++)
      CHECKfscanf (fscanf (EF, "%*d %lf", &ef[i]), efFile);

   /* Closes the file. */
   CHECKfclose (fclose (EF), efFile);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (efFile);

} /* readFermiEnergy */


/* ********************************************************* */
/* Reads the first orbital index of each atom in the unit    */
/* cell from the '.orb' file.                                */
static void readOrbitalIndex ()
{
   register int i, len;
   char *orbFile;
   FILE *ORB;

   /* Sets the '.orb' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   orbFile = CHECKmalloc ((len + 5) * sizeof (char));
   orbFile[0] = '\0';
   sprintf (orbFile, "%s%s.orb", FCdir, sysLabel);

   /* Opens the '.orb' file. */
   printf ("\n    reading \"%s\" file... ", orbFile);
   ORB = CHECKfopen (orbFile, "r");

   /* Reads the first orbital index of each atom. */
   orbIdx = CHECKmalloc ((nAtoms + 1) * sizeof (int));
   for (i = 0; i < nAtoms + 1; i++)
      CHECKfscanf (fscanf (ORB, "%d", &orbIdx[i]), orbFile);

   /* Closes the file. */
   CHECKfclose (fclose (ORB), orbFile);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Number of basis orbitals from unit cell. */
   no_u = orbIdx[nAtoms];
   printf ("\n    Basis orbitals per unit cell:\t\t%d\n", no_u);

   /* First orbital of the first dynamic atom. */
   printf ("    First orbital from dynamic atom:\t\t%d\n",
	   orbIdx[FCfirst-1] + 1);

   /* Dynamic atoms orbitals quantity. */
   printf ("    Number of dynamic atoms basis orbitals:\t%d\n",
	   orbIdx[FClast] - orbIdx[FCfirst-1]);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (orbFile);

} /* readOrbitalIndex */


/* ********************************************************* */
/* Runs a script that collects required informations from FC */
/* input 'fdf' file and writes them at 'inputFC.in' file.    */
/* Calls 'assignGlobalVar' function to assign global         */
/* variables.                                                */
void PHONreadFCfdf (char *exec, char *FCpath, char *FCinput,
		    int calcType, char *FCsplit, int *nDynTot,
		    int *nDynOrb, int *spinPol)
{
   register int i, len;
   char calc[2];
   char *scriptCall;

   /* For calling the script, one should use a string     */
   /* like this: "[script name with path] [FC directory]" */
   len = strlen(FCpath) + strlen(FCinput) + strlen(FCsplit);
   scriptCall = CHECKmalloc ((len + 19) * sizeof (char));
   scriptCall[0] = '\0';
   calc[0] = '\0';
   sprintf(&calc[0], "%d", calcType);
   strcpy (scriptCall, "buildInput.sh ");
   strcat (scriptCall, FCpath);
   strcat (scriptCall, " ");
   strcat (scriptCall, FCinput);
   strcat (scriptCall, " ");
   strcat (scriptCall, calc);
   strcat (scriptCall, " ");
   strcat (scriptCall, FCsplit);

   /* Assigns the work directory global variable. */
   len = strlen (exec);
   workDir = CHECKmalloc ((len - 9) * sizeof (char));
   for (i = 0; i < len - 10; i++)
      workDir[i] = exec[i];
   workDir[i] = '\0';

   /* Assigns the force constants directory global variable. */
   len = strlen (FCpath);
   if (FCpath[len - 1] == '/') {
      FCdir = CHECKmalloc ((len + 1) * sizeof (char));
      strcpy (FCdir, FCpath);
   }
   else {
      FCdir = CHECKmalloc ((len + 2) * sizeof (char));
      strcpy (FCdir, FCpath);
      FCdir[len] = '/';
      FCdir[len + 1] = '\0';
   }

   /* Runs the script that collects required informations from */
   /* the input 'fdf' file and puts at 'inputFC.in' file.      */
   i = system (scriptCall);
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   if (i != 0) {
      fprintf (stderr, " ERROR: problem when trying to run");
      fprintf (stderr, " 'buildInput.sh' script!\n\n");
      exit (EXIT_FAILURE);
   }

   /* Reads the file 'inputFC.in' and assigns global variables. */
   assignGlobalVar ();

   if (calcType == 1) { /* 'full' calculation */

      /* Reads the Fermi energies. */
      printf ("\n Gets Fermi energies from (un)displaced systems:\n");
      setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
      readFermiEnergy ();

      /* Gets the first orbital index of each atom. */
      printf ("\n Gets basis orbitals dimensions info:\n");
      setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
      readOrbitalIndex ();

      /* Returns the total number of orbitals for 'e-ph' coupling matrix. */
      *nDynOrb = orbIdx[FClast] - orbIdx[FCfirst-1];

   }

   /* Returns the dimension of 'FC' matrix. */
   *nDynTot = 3 * nDyn;
   *spinPol = nspin;

   /* Removes 'inputFC.in' file. */
   scriptCall[0] = '\0';
   strcpy (scriptCall, "rm ");
   strcat (scriptCall, exec);
   scriptCall[strlen (scriptCall)-10] = '\0'; /* removes 'exec' */
   strcat (scriptCall, "inputFC.in");
   i = system (scriptCall);

   /* Frees memory. */
   free (scriptCall);

} /* PHONreadFCfdf */


/* ********************************************************* */
/* Removes egg-box effect by imposing force conservation.    */
static void rmEggBox (double *fullFCM)
{
   register int i, j, k, dyn;
   double sum;

   for (j = 0; j < 3 * nDyn; j++) {
      dyn = FCfirst - 1 + j / 3; /* dynamic atom */
      for (k = 0; k < 3; k++) { /* coordinates (x,y,z) */
	 fullFCM[idx(dyn*3+k,j,3*nAtoms)] = 0.0;
	 sum = 0.0;
	 for (i = 0; i < nAtoms; i++) /* sum over all atoms coord 'k' */
	    sum = sum + fullFCM[idx(i*3+k,j,3*nAtoms)];
	 fullFCM[idx(dyn*3+k,j,3*nAtoms)] = - 1.0 * sum;
      }
   }

} /* rmEggBox */


/* ********************************************************* */
/* Symmetrizes and mass-scales the 'FC' matrix.              */
static void symmetrizesAndMassScale (double *FC)
{
   register int i, j;
   double aux;

   for (j = 0; j < 3 * nDyn; j++)
      for (i = j; i < 3 * nDyn; i++) {
	 /* Computes the mean. */
	 aux = (FC[idx(i,j,3*nDyn)] + FC[idx(j,i,3*nDyn)]) / 2.0;

	 /* Mass-scales. */
	 aux = aux / sqrt (dynAtoms[i/3].atom.A * dynAtoms[j/3].atom.A);

	 /* Simetrizes. */
	 FC[idx(i,j,3*nDyn)] = FC[idx(j,i,3*nDyn)] = aux;
      }

} /* symmetrizesAndMassScale */


/* ********************************************************* */
/* Reads the SIESTA force constants matrix and computes      */
/* phonon modes and frequencies with finite differences.     */
void PHONfreq (double *EigVec, double *EigVal)
{
   register int i, j, len;
   double cst;
   double *fullFCneg, *fullFCpos;
   char check[25];
   char *FCMfile;
   FILE *FCM;

   /* Sets the SIESTA FC matrix file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   FCMfile = CHECKmalloc ((len + 4) * sizeof (char));
   sprintf (FCMfile, "%s%s.FC", FCdir, sysLabel);

   /* Opens the SIESTA FC matrix file. */
   printf ("\n Reading %s file... ", FCMfile);
   FCM = CHECKfopen (FCMfile, "r");

   /* Checks if the file starts correctly. */
   if ((fgets(check,24,FCM)==NULL) ||
       (strncmp(check,"Force constants matrix",22)!=0)) {
      fprintf (stderr,
	       " ERROR: the file %s is not written correctly!\n\n",
	       FCMfile);
      exit (EXIT_FAILURE);
   }

   /* FC matrix with all atoms                             */
   /* [negative (minus) and positive (plus) displacement]. */
   fullFCneg = CHECKmalloc (3 * nAtoms * 3 * nDyn * sizeof (double));
   fullFCpos = CHECKmalloc (3 * nAtoms * 3 * nDyn * sizeof (double));

   /* Reads the SIESTA FC matrix file       */
   /* (obs.: matrix in column-major order). */
   for (j = 0; j < 3 * nDyn; j++) {
      /* Negative displacement. */
      for (i = 0; i < 3 * nAtoms; i++)
	 CHECKfscanf (fscanf (FCM, "%lE",
			      &fullFCneg[idx(i,j,3*nAtoms)]), FCMfile);

      /* Positive displacement. */
      for (i = 0; i < 3 * nAtoms; i++)
	 CHECKfscanf (fscanf (FCM, "%lE",
			      &fullFCpos[idx(i,j,3*nAtoms)]), FCMfile);
   }

   /* Closes the SIESTA FC matrix file. */
   CHECKfclose (fclose (FCM), FCMfile);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Removes egg-box effect. */
   printf ("\n Removing egg-box effect... ");
   rmEggBox (fullFCneg);
   rmEggBox (fullFCpos);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* The SIESTA FC matrix file contains the values of force    */
   /* by displacement [eV/Ang^2]. In order to compute the       */
   /* Hessian (force constants matrix) with finite differences, */
   /* one has only to subtract the value corresponding to the   */
   /* negative displacement by the value corresponding to the   */
   /* positive displacement and divide by 2.                    */
   printf ("\n Reducing to dynamic atoms dimension");
   printf (" and computing finite differences... ");
   len = (FCfirst - 1) * 3;
   for (j = 0; j < 3 * nDyn; j++)
      for (i = len; i < FClast * 3; i++)
	 EigVec[idx(i-len,j,3*nDyn)] =
	    (fullFCneg[idx(i,j,3*nAtoms)] +
	     fullFCpos[idx(i,j,3*nAtoms)]) / 2.0;
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Symmetrizes and mass-scales the matrix. */
   printf ("\n \"Symmetrizing\" and mass-scaling the FC matrix... ");
   symmetrizesAndMassScale (EigVec);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Prints on screen the reduced,          */
   /* mass-scaled and symmetrized FC matrix. */
   printf ("\n Reduced force constants matrix:\n\n");
   for (i = 0; i < 3 * nDyn; i++) {
      for (j = 0; j < 3 * nDyn; j++)
	 printf (" % .5e  ", EigVec[idx(i,j,3*nDyn)]);
      printf ("\n");
   }
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Computes all eigenvalues and eigenvectors. */
   CHECKdsyevd (3 * nDyn, EigVec, EigVal);

   /* Prints on screen the phonon energies. */
   printf ("\n Phonon energies (eV):\n\n");
   cst = hbar * sqrt (1.0e20 * eV2joule / amu2kg);
   for (i = 0; i < 3 * nDyn; i++) {
      if (EigVal[i] < 0.0) {
	 printf ("  %d % .5e\n", i + 1, - cst * sqrt(-EigVal[i]));
	 EigVal[i] = 0.0;
      } else {
	 EigVal[i] = cst * sqrt(EigVal[i]);
	 printf ("  %d % .5e\n", i + 1, EigVal[i]);
      }
   }
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Prints on screen the phonon modes. */
   printf ("\n Normalized phonon modes (Ang*amu^0.5):\n\n");
   for (i = 0; i < 3 * nDyn; i++)
      printf (" %7d       ", i+1);
   printf ("\n");
   for (i = 0; i < 3 * nDyn; i++) {
      for (j = 0; j < 3 * nDyn; j++)
	 printf (" % .5e  ", EigVec[idx(i,j,3*nDyn)]);
      printf ("\n");
   }
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (fullFCneg);
   free (fullFCpos);
   free (FCMfile);

} /* PHONfreq */


/* ********************************************************* */
/* Reads the SIESTA 'xyz' file and write 'xyz's files for    */
/* each computed phonon mode.                                */
void PHONjmolVib (double *EigVec)
{
   register int i, j, k, len;
   int aux;
   xyzcoord *coord;
   char *XYZfile, *JMOLfile;
   FILE *XYZ, *JMOL;

   /* Sets the SIESTA 'xyz' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   XYZfile = CHECKmalloc ((len + 5) * sizeof (char));
   sprintf (XYZfile, "%s%s.xyz", FCdir, sysLabel);

   /* Opens the SIESTA 'xyz' file. */
   printf ("\n Reading %s file... ", XYZfile);
   XYZ = CHECKfopen (XYZfile, "r");

   /* Checks if the file starts correctly. */
   CHECKfscanf (fscanf (XYZ, "%d", &aux), XYZfile);
   if (aux != nAtoms) {
      fprintf (stderr,
	       " ERROR: the file %s is not written correctly!\n\n",
	       XYZfile);
      exit (EXIT_FAILURE);
   }

   /* Allocs structure for atomic species and coordinates. */
   coord = CHECKmalloc (nAtoms * sizeof (xyzcoord));

   /* Reads the SIESTA 'xyz' file. */
   for (i = 0; i < nAtoms; i++)
      CHECKfscanf (fscanf (XYZ, "%s %le %le %le",
			   coord[i].name, &coord[i].x,
			   &coord[i].y, &coord[i].z), XYZfile);

   /* Closes the SIESTA 'xyz' file. */
   CHECKfclose (fclose (XYZ), XYZfile);
   printf ("ok!\n\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Sets the path for the output Jmol files. */
   len = strlen (workDir);
   len += strlen (sysLabel);
   JMOLfile = CHECKmalloc ((len + 12) * sizeof (char));

   for (i = 0, j = 1; i < 3*nDyn; i++) {
         
	 /* Opens the JMOL 'xyz' output file. */
	 sprintf (JMOLfile, "%s%sJMOL%d.xyz", workDir, sysLabel, j);
	 printf (" Writing %s file... ", JMOLfile);
	 JMOL = CHECKfopen (JMOLfile, "w");

	 fprintf (JMOL, "   %d\n\n", nAtoms);

	 /* Writes only the coordinates. */
	 for (k = 0; k < FCfirst - 1; k++)
	    fprintf (JMOL, "%s\t% e\t% e\t% e\n", coord[k].name,
		     coord[k].x, coord[k].y, coord[k].z);

	 /* Writes the coordinates and the normalized phonon mode. */
	 for (k = FCfirst - 1; k < FClast; k++)
	    fprintf (JMOL, "%s\t% e\t% e\t% e\t% e\t% e\t% e\n",
		     coord[k].name, coord[k].x, coord[k].y, coord[k].z,
		     EigVec[idx(3*(k-FCfirst+1),i,3*nDyn)],
		     EigVec[idx(3*(k-FCfirst+1)+1,i,3*nDyn)],
		     EigVec[idx(3*(k-FCfirst+1)+2,i,3*nDyn)]);

	 /* Writes only the coordinates. */
	 for (k = FClast; k < nAtoms; k++)
	    fprintf (JMOL, "%s\t% e\t% e\t% e\n", coord[k].name,
		     coord[k].x, coord[k].y, coord[k].z);

	 /* Closes the JMOL 'xyz' output file. */
	 CHECKfclose (fclose (JMOL), JMOLfile);
	 printf ("ok!\n");
	 setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

	 j++;

   } /* for (i = 0, j = 1; ...  */

   /* Frees memory. */
   free (coord);
   free (XYZfile);
   free (JMOLfile);

} /* PHONjmolVib */


/* ********************************************************* */
/* Reads the Hamiltonian and the overlap matrices from       */
/* '.gHS' file.                                              */
static void readHSfile (char *HSfile, double *H, double *S, int efIdx)
{
   register int i, j, k, s, foo;
   int CKno_u, CKnspin, maxnhtot;
   double *Hsparse, *Ssparse;
   int *numh, *listh;
   FILE *gHS;

   /* Opens the '.gHS' binary file. */
   printf ("    reading \"%s\" file... ", HSfile);
   gHS = CHECKfopen (HSfile, "rb");

   /* Reads matrices dimensions and checks file consistency. */
   fseek (gHS, 0L, SEEK_SET);
   CHECKfread (fread (&CKno_u, sizeof (int), 1, gHS), 1, HSfile);
   CHECKfread (fread (&CKnspin, sizeof (int), 1, gHS), 1, HSfile);
   CHECKfread (fread (&maxnhtot, sizeof (int), 1, gHS), 1, HSfile);
   if ((CKno_u != no_u) || (CKnspin != nspin)) {
      fprintf (stderr,
	       " ERROR: the file %s is not written correctly!\n\n",
	       HSfile);
      exit (EXIT_FAILURE);
   }

   /* Allocates memory. */
   Hsparse = UTILdoubleVector (nspin * maxnhtot);
   Ssparse = UTILdoubleVector (maxnhtot);
   numh = UTILintVector (no_u);
   listh = UTILintVector (maxnhtot);

   /* Reads the number of nonzero elements of each row of H. */
   CHECKfread (fread (numh, no_u * sizeof (int), 1, gHS), 1, HSfile);

   /* Reads the nonzero Hamiltonian-matrix */
   /* element column indexes for each row. */
   for (i = 0, k = 0; i < no_u; i++) {
      CHECKfread (fread (&(listh[k]),
			 numh[i]*sizeof(int), 1, gHS), 1, HSfile);
      k = k + numh[i];
   }

   /* Reads the Hamiltonian matrix in sparse form. */
   for (j = 0, k = 0; j < nspin; j++)
      for (i = 0; i < no_u; i++) {
   	 CHECKfread (fread (&(Hsparse[k]),
			    numh[i]*sizeof(double), 1, gHS), 1, HSfile);
   	 k = k + numh[i];
      }

   /* Reads the overlap matrix in sparse form. */
   for (i = 0, k = 0; i < no_u; i++) {
      CHECKfread (fread (&(Ssparse[k]),
			 numh[i]*sizeof (double), 1, gHS), 1, HSfile);
      k = k + numh[i];
   }

   /* Assigns nonzero elements from 'Hsparse' and 'Ssparse'.    */
   /* The index 'k' runs over 'Hsparse' and 'Ssparse' matrices. */
   for (i = 0, k = 0; i < no_u; i++)
      for (j = 0; j < numh[i]; j++, k++) {
   	 foo = (listh[k] - 1) % no_u; /* column index */
   	 H[idx3d(i,foo,0,no_u,no_u)] += rydberg2eV * Hsparse[k];
   	 S[idx(i,foo,no_u)] += Ssparse[k];
      }
   for (s = 1; s < nspin; s++) /* 'H' is bigger if 'nspin' > 1 */
      for (i = 0, k = 0; i < no_u; i++)
   	 for (j = 0; j < numh[i]; j++, k++) {
  	    foo = (listh[k] - 1) % no_u; /* column index */
   	    H[idx3d(i,foo,s,no_u,no_u)] += rydberg2eV * Hsparse[k];
   	 }

   /* "Shifts" the Fermi energy to 0. */
   for (s = 0; s < nspin; s++)
      for (i = 0; i < no_u; i++)
   	 for (j = 0; j < no_u; j++)
   	    H[idx3d(i,j,s,no_u,no_u)] -= ef[efIdx] * S[idx(i,j,no_u)];

   /* Closes file. */
   CHECKfclose (fclose (gHS), HSfile);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (numh);
   free (listh);
   free (Hsparse);
   free (Ssparse);

} /* readHSfile */


/* ********************************************************* */
/* Computes Hamiltonian derivative matrix by reading the     */
/* Hamiltonian and overlap matrices from '.gHs' files for    */
/* the displaced system.                                     */
static void deltaH (double *dH, double *S0)
{
   register int i, j, k, s, len;
   double *Hm, *Hp, *S;
   char *Hfile;

   /* Allocates memory. */
   Hm = UTILdoubleVector (3 * nDyn * nspin * no_u * no_u);
   Hp = UTILdoubleVector (3 * nDyn * nspin * no_u * no_u);
   S = CHECKmalloc (no_u * no_u * sizeof (double));

   /* Sets the '.gHS' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   Hfile = CHECKmalloc ((len + 9) * sizeof (char));

   for (k = 0; k < 3 * nDyn; k++) {
      /* 'H(-Q)' */
      Hfile[0] = '\0';
      sprintf (Hfile, "%s%s_%.3d.gHS", FCdir, sysLabel, 2 * k + 1);
      UTILresetDoubleVector (no_u * no_u, S);
      readHSfile (Hfile, &Hm[idx3d(0,0,k*nspin,no_u,no_u)], S, 2*k+1);

      /* 'H(Q)' */
      Hfile[0] = '\0';
      sprintf (Hfile, "%s%s_%.3d.gHS", FCdir, sysLabel, 2 * k + 2);
      UTILresetDoubleVector (no_u * no_u, S);
      readHSfile (Hfile, &Hp[idx3d(0,0,k*nspin,no_u,no_u)], S, 2*k+2);

      /* 'dH = {H(Q) - (ef(Q)-ef0)*S0 - [H(-Q)-(ef(-Q)-ef0)*S0]} / 2Q' */
      /* or, simplifying: 'dH = {H(Q)-H(-Q)-[ef(Q)-ef(-Q)]*S0} / 2Q'   */
      for (s = 0; s < nspin; s++)
      	 for (i = 0; i < no_u; i++)
      	    for (j = 0; j < no_u; j++)
      	       dH[idx3d(i,j,k*nspin+s,no_u,no_u)] = 
		  (Hp[idx3d(i,j,k*nspin+s,no_u,no_u)]
		   - Hm[idx3d(i,j,k*nspin+s,no_u,no_u)]
		   - (ef[2*k+2] - ef[2*k+1]) * S0[idx(i,j,no_u)])
		  / (2.0 * FCdispl);

   }

   /* Frees memory. */
   free (Hfile);
   free (Hp);
   free (Hm);
   free (S);

} /* deltaH */


/* ********************************************************* */
/* Reads the overlap matrix from '.onlyS' binary file.       */
static void readOnlyS (char *Sfile, double *S)
{
   register int i, j, k, foo;
   int CKno_u, maxnhtot;
   double *Ssparse, *Sbig;
   int *numh, *listh;
   FILE *OnlyS;

   /* Opens the '.onlyS' binary file. */
   printf ("    reading \"%s\" file... ", Sfile);
   OnlyS = CHECKfopen (Sfile, "rb");

   /* Reads matrices dimensions and checks file consistency. */
   fseek (OnlyS, 0L, SEEK_SET);
   CHECKfread (fread (&CKno_u, sizeof (int), 1, OnlyS), 1, Sfile);
   CHECKfread (fread (&maxnhtot, sizeof (int), 1, OnlyS), 1, Sfile);
   if (CKno_u != 2 * no_u) {
      fprintf (stderr,
	       " ERROR: the file %s is not written correctly!\n\n",
	       Sfile);
      exit (EXIT_FAILURE);
   }

   /* Allocates memory. */
   Ssparse = UTILdoubleVector (maxnhtot);
   Sbig = UTILdoubleVector (2 * no_u * 2 * no_u);
   numh = UTILintVector (2 * no_u);
   listh = UTILintVector (maxnhtot);

   /* Reads the number of nonzero elements of each row of S. */
   CHECKfread (fread (numh, 2*no_u*sizeof(int), 1, OnlyS), 1, Sfile);

   /* Reads the nonzero overlap-matrix     */
   /* element column indexes for each row. */
   for (i = 0, k = 0; i < 2 * no_u; i++) {
      CHECKfread (fread (&(listh[k]), numh[i]*sizeof(int),
			 1, OnlyS), 1, Sfile);
      k = k + numh[i];
   }

   /* Reads the overlap matrix in sparse form. */
   for (i = 0, k = 0; i < 2 * no_u; i++) {
      CHECKfread (fread (&(Ssparse[k]), numh[i]*sizeof(double),
			 1, OnlyS), 1, Sfile);
      k = k + numh[i];
   }

   /* Assigns nonzero elements from 'Ssparse'.  */
   /* The index 'k' runs over 'Ssparse' matrix. */
   for (i = 0, k = 0; i < 2 * no_u; i++)
      for (j = 0; j < numh[i]; j++, k++) {
   	 foo = (listh[k] - 1) % (2 * no_u); /* column index */
   	 Sbig[idx(i,foo,2*no_u)] += Ssparse[k];
      }

   /* Picks only the terms <i|j'> from dynamic atoms. */
   for (j = 0; j < no_u; j++)
      for (i = 0; i < no_u; i++)
   	 S[idx(i,j,no_u)] = (Sbig[idx(i,j+no_u,2*no_u)]
			     + Sbig[idx(j+no_u,i,2*no_u)]) / 2.0;

   /* Closes file. */
   CHECKfclose (fclose (OnlyS), Sfile);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (numh);
   free (listh);
   free (Ssparse);
   free (Sbig);

} /* readOnlyS */


/* ********************************************************* */
/* Computes 'dS = [ <i|j(Q)> - <i|j(-Q)> ] / 2Q' by reading  */
/* the overlap matrix from '.onlyS' files for the displaced  */
/* system.                                                   */
static void deltaS (double *dS)
{
   register int i, j, k, len;
   double *Sm, *Sp;
   char *Sfile;

   /* Allocates memory and initializes with '0'. */
   Sm = UTILdoubleVector (3 * no_u * no_u);
   Sp = UTILdoubleVector (3 * no_u * no_u);

   /* Sets the '.onlyS' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   Sfile = CHECKmalloc ((len + 9) * sizeof (char));

   for (k = 0; k < 3; k++) {
      /* '<i|j(-Q)>' */
      Sfile[0] = '\0';
      sprintf (Sfile, "%s%s_%d.onlyS", FCdir, sysLabel, 2 * k + 1);
      readOnlyS (Sfile, &Sm[idx3d(0,0,k,no_u,no_u)]);

      /* '<i|j(Q)>' */
      Sfile[0] = '\0';
      sprintf (Sfile, "%s%s_%d.onlyS", FCdir, sysLabel, 2 * k + 2);
      readOnlyS (Sfile, &Sp[idx3d(0,0,k,no_u,no_u)]);

      /* 'dS = [ <i|j(Q)> - <i|j(-Q)> ] / 2Q' */
      for (i = 0; i < no_u; i++)
	 for (j = 0; j < no_u; j++)
	    dS[idx3d(i,j,k,no_u,no_u)] =
	       (Sp[idx3d(i,j,k,no_u,no_u)] - Sm[idx3d(i,j,k,no_u,no_u)])
	       / (2.0 * FCdispl);

   }

   /* Frees memory. */
   free (Sfile);
   free (Sp);
   free (Sm);

} /* deltaS */


/* ********************************************************* */
/* Applies a correction due to the change in basis orbitals  */
/* with displacement:                                        */
/*         'dH = dH - dS*S0^-1*H0 - H0*S0^-1*(dS)^T'         */
static void dHCorrection (double *dH, double *H0, double *S0)
{
   register int i, j, k, s, coord, h;
   int *ipiv;
   double alpha, beta;
   double *dStot, *dS, *invS0, *Aux;

   /* Computes 'dS = [ <i|j(Q)> - <i|j(-Q)> ] / 2Q'. */
   dStot = UTILdoubleVector (3 * no_u * no_u);
   deltaS (dStot);

   /* Initializes 'invS0' with 'S0'. */
   invS0 = UTILdoubleVector (no_u * no_u);
   UTILcopyVector (invS0, S0, no_u * no_u);

   /* Computes 'invS0 = S0^-1'. */
   ipiv = CHECKmalloc (no_u * sizeof (int));
   CHECKdgetrf (no_u, invS0, ipiv); /* triangular matrix factorization */
   CHECKdgetri (no_u, invS0, ipiv); /* matrix inversion */

   /* Computes 'dH = dH - dS*S0^-1*H0 - H0*S0^-1*(dS)^T' */
   /* for each displacement direction.                   */
   printf ("\n    correcting Hamiltonian derivatives elements... ");
   dS = CHECKmalloc (no_u * no_u * sizeof (double));
   Aux = CHECKmalloc (no_u * no_u * sizeof (double));
   for (k = 0; k < nDyn; k++) {
      UTILresetDoubleVector (no_u * no_u, dS);
      for (coord = 0; coord < 3; coord++) { /* xyz */

	 /* The only non-zero elements from 'dS' matrix are those    */
	 /* corresponding to the orbitals from the dynamic atom 'k'. */
	 for (j = orbIdx[FCfirst+k-1]; j < orbIdx[FCfirst+k]; j++)
	    for (i = 0; i < no_u; i++)
	       dS[idx(i,j,no_u)] = dStot[idx3d(i,j,coord,no_u,no_u)];

	 for (s = 0; s < nspin; s++) {

	    /* 'dH' index. */
	    h = (k * 3 + coord) * nspin+s;

	    /* 'Aux = - dS * S0^-1' */
	    alpha = - 1.0;
	    beta = 0.0; /* as 'beta = 0' -> no need to reset 'Aux' */
	    dgemm ("T","N", &no_u, &no_u, &no_u, &alpha, dS,
	    	   &no_u, invS0, &no_u, &beta, Aux, &no_u);

	    /* 'dH = dH + Aux * H0' */
	    alpha = beta = 1.0;
	    dgemm ("N","N", &no_u, &no_u, &no_u, &alpha, Aux, &no_u,
		   &H0[idx3d(0,0,s,no_u,no_u)], &no_u, &beta,
		   &dH[idx3d(0,0,h,no_u,no_u)], &no_u);

	    /* 'Aux = H0 * S0^-1' */
	    beta = 0.0;
	    dgemm ("N","N", &no_u, &no_u, &no_u, &alpha,
		   &H0[idx3d(0,0,s,no_u,no_u)],  &no_u,
		   invS0, &no_u, &beta, Aux, &no_u);

	    /* 'dH = dH - Aux * (dS)^T' */
	    alpha = - 1.0;
	    beta = 1.0;
	    dgemm ("N","N", &no_u, &no_u, &no_u, &alpha, Aux, &no_u, dS,
	    	   &no_u, &beta, &dH[idx3d(0,0,h,no_u,no_u)], &no_u);
	 }
      }
   }
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (ipiv);
   free (dStot);
   free (dS);
   free (invS0);
   free (Aux);

} /* dHCorrection */


/* ********************************************************* */
/* Having calculated the phonon energies ('EigVal') and      */
/* modes ('EigVec') and the Hamiltonian derivatives ('dH'),  */
/* it computes the elements of the electron-phonon coupling  */
/* matrix.                                                   */
static void eph (double *EigVec, double *EigVal,
		 double *dH, double *Meph)
{
   register int i, j, k, s, l, firstOrb, nOrb;
   double cst;

   /* Computes each element of 'Meph'. */
   printf ("    computing the electron-phonon coupling elements... ");
   firstOrb = orbIdx[FCfirst - 1]; /* first orb of first dyn atom */
   nOrb = orbIdx[FClast] - orbIdx[FCfirst - 1]; /* number of dyn orbs */
   cst = hbar * sqrt (1.0e20 * eV2joule / amu2kg);
   for (l = 0; l < 3 * nDyn; l++) /* modes */
      for (k = 0; k < 3 * nDyn; k++) /* coordinates */
	 for (s = 0; s < nspin; s++) /* spin */
	    for (j = 0; j < nOrb; j++)
	       for (i = 0; i < nOrb; i++)
		  Meph[idx3d(i,j,l*nspin+s,nOrb,nOrb)] +=
		     dH[idx3d(firstOrb+i,firstOrb+j,k*nspin+s,no_u,no_u)]
		     * EigVec[idx(k,l,3*nDyn)] * cst
		     / sqrt (2 * dynAtoms[k/3].atom.A * EigVal[l]);
   printf ("ok!\n\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

} /* eph */


/* ********************************************************* */
/* For each phonon energy ('EigVal'), ouputs the             */
/* corresponding electron-phonon coupling matrix.            */
static void ephOut (double *EigVal, double *Meph)
{
   register int i, j, l, s, len;
   int nOrb, foo;
   char *ephFile, *ephFileB;
   FILE *EPH, *EPHb;

   /* Sets the '.Meph' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len += strlen (sysLabel);
   ephFile = CHECKmalloc ((len + 6) * sizeof (char));
   ephFile[0] = '\0';
   sprintf (ephFile, "%s%s.Meph", FCdir, sysLabel);
   ephFileB = CHECKmalloc ((len + 7) * sizeof (char));
   ephFileB[0] = '\0';
   sprintf (ephFileB, "%s%s.bMeph", FCdir, sysLabel);

   /* Opens the '.Meph' file. */
   printf (" Writing electron-phonon coupling matrix at \"%s\" file... ",
	   ephFile);
   EPH = CHECKfopen (ephFile, "w");
   EPHb = CHECKfopen (ephFileB, "wb");

   /* Writes the dimensions. */
   nOrb = orbIdx[FClast] - orbIdx[FCfirst - 1]; /* number of dyn orbs */
   fprintf (EPH, "%d  %d  %d  %d  %d\n\n", nspin, nDyn,
	    nOrb, orbIdx[FCfirst-1]+1, orbIdx[FCfirst-1]+nOrb);
   fwrite (&nspin, sizeof(int), 1, EPHb);
   fwrite (&nDyn, sizeof(int), 1, EPHb);
   fwrite (&nOrb, sizeof(int), 1, EPHb);
   foo = orbIdx[FCfirst-1]+1;
   fwrite (&foo, sizeof(int), 1, EPHb);
   foo = orbIdx[FCfirst-1]+nOrb;
   fwrite (&foo, sizeof(int), 1, EPHb);

   /* Prints the phonon frequencies. */
   for (l = 3 * nDyn - 1; l >= 0; l--) {
      fprintf (EPH, "%.10e  ", EigVal[l]);
      fwrite (&(EigVal[l]), sizeof(double), 1, EPHb);
   }
   fprintf (EPH, "\n\n");

   /* Prints the electron-phonon coupling matrix. */
   for (l = 3 * nDyn - 1; l >= 0; l--)
      if (EigVal[l] > 0.0)
	 for (s = 0; s < nspin; s++) {
	    for (i = 0; i < nOrb; i++) {
	       for (j = 0; j < nOrb; j ++)
		  fprintf (EPH, " % .15e",
			   Meph[idx3d(i,j,l*nspin+s,nOrb,nOrb)]);
	       fprintf (EPH, "\n");
	    }
	    fprintf (EPH, "\n");
	    fwrite (&(Meph[idx3d(0,0,l*nspin+s,nOrb,nOrb)]),
		    sizeof(double), nOrb*nOrb, EPHb);
	 }

   /* /\* for (l = 0; l < 3 * nDyn; l++) *\/ */
   /* for (l = 3 * nDyn - 1; l >= 0; l--) */
   /*    for (s = 0; s < nspin; s++) { */
   /* 	 /\* Writes spin and phonon frequency. *\/ */
   /* 	 fprintf (EPH, "%d %.10E\n\n", s + 1, EigVal[l]); */
   /* 	 for (i = 0; i < nOrb; i++) { */
   /* 	    for (j = 0; j < nOrb; j ++) */
   /* 	       fprintf (EPH, " (% .15e,0.000000000000000e00)", */
   /*                   Meph[idx3d(i,j,l*nspin+s,nOrb,nOrb)]); */
   /* 	    fprintf (EPH, "\n"); */
   /* 	 } */
   /* 	 fprintf (EPH, "\n"); */
   /*    } */

   /* Closes the file. */
   CHECKfclose (fclose (EPH), ephFile);
   CHECKfclose (fclose (EPHb), ephFileB);
   printf ("ok!\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */

   /* Frees memory. */
   free (ephFile);
   free (ephFileB);

} /* ephOut */


/* ********************************************************* */
/* Computes the electron-phonon coupling matrices.           */
void PHONephCoupling (double *EigVec, double *EigVal, double *Meph)
{
   register int len;
   double *H0, *S0, *dH;
   char *HSfile;

   /* Sets the '.only(X)' file name with 'FCdir' path. */
   len = strlen (FCdir);
   len = len + strlen (sysLabel);
   HSfile = CHECKmalloc ((len + 9) * sizeof (char));

   /* Reads 'H0' and 'S0' matrices (non-displaced system). */
   printf ("\n 'H0' and 'S0' matrices (non-displaced system):\n\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   H0 = UTILdoubleVector (nspin * no_u * no_u);
   S0 = UTILdoubleVector (no_u * no_u);
   HSfile[0] = '\0';
   sprintf (HSfile, "%s%s_%.3d.gHS", FCdir, sysLabel, 0);
   readHSfile (HSfile, H0, S0, 0);

   /* Computes 'dH={H(Q)-(ef(Q)-ef0)*S0-[H(-Q)-(ef(-Q)-ef0)*S0]}/2Q'. */
   printf ("\n 'H' matrix derivative:\n\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   dH = UTILdoubleVector (3 * nDyn * nspin * no_u * no_u);
   deltaH (dH, S0);

   /* Applies a correction due to the change in basis orbitals with */
   /* displacements: 'dH = dH - dS * S^-1 * H0 - H0 * S^-1 * dS'.   */
   printf ("\n 'dH' correction due to the changes in basis orbitals:");
   printf ("\n\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   dHCorrection (dH, H0, S0);

   /* Computes the electron-phonon coupling matrices. */
   printf ("\n Computes electron-phonon coupling matrix.\n\n");
   setvbuf (stdout, NULL, _IONBF, 0); /* print now! */
   eph (EigVec, EigVal, dH, Meph);

   /* Outputs the electron-phonon coupling matrices. */
   ephOut (EigVal, Meph);

   /* Frees memory. */
   free (H0);
   free (S0);
   free (dH);
   free (HSfile);

} /* PHONephCoupling */


/* ************************ Drafts ************************* */

