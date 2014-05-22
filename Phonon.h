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
/**   File: Phonon.h                                        **/
/**                                                         **/
/**   Versions: 1 - 10/10/2012                              **/
/**             2 - 08/01/2013                              **/
/**                                                         **/
/**  *****************************************************  **/
/**  Interface of phonon analysis to:                       **/
/**   - read the output data from a siesta force            **/
/**   constants (FC) calculation runned with the flag       **/
/**   'PB.FCwriteHS .true.' at the fdf input file;          **/
/**   - build the force constants matrix and computes the   **/
/**   phonon frequencies and modes by finite differences;   **/
/**   - compute the hamiltonian and overlap matrices        **/
/**   derivatives with respect to the dynamic atoms         **/
/**   displacements and apply suitable corrections and      **/
/**   calculate the electron-phonon coupling matrix.        **/
/**  *****************************************************  **/


/* For calling from fortran prograns. */
#ifdef FORTRAN
#define PHONheader phonheader_
#define PHONreadFCfdf phonreadfcfdf_
#define PHONfreq phonfreq_
#define PHONephCoupling phonephcoupling_
#endif

/* Collects required informations from FC input 'fdf' file. */
void PHONreadFCfdf (char *exec, char *FCpath, char *FCinput,
		    int calcType, char *FCsplit, int *nDynTot,
		    int *nDynOrb, int *spinPol);

/* Computes phonon frequencies and modes. */
void PHONfreq (double *EigVec, double *EigVal);

/* Writes a 'xyz' file for each computed phonon mode. */
void PHONjmolVib (double *EigVec);

/* Computes electron-phonon coupling matrices. */
void PHONephCoupling (double *EigVec, double *EigVal, double *Meph);
