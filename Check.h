/**  *****************************************************  **/
/**             ** Phonon Vibration Analysis **             **/
/**                                                         **/
/**                    **  Version 2  **                    **/
/**                                                         **/
/**   By: Pedro Brandimarte (brandimarte@gmail.com) and     **/
/**       Alexandre Reily Rocha (reilya@ift.unesp.br)       **/
/**                                                         **/
/**  *****************************************************  **/
/**  Interface for alternative calls of well known          **/
/**  functions from 'stdio.h' and 'stdlib.h', and also      **/
/**  functions from 'LAPACK' package, to avoid code         **/
/**  repetition when checking errors.                       **/
/**  *****************************************************  **/


/* Calls Lapack rotine 'dsyevd' for computing all eigenvalues and     */
/* eigenvectors of a real symmetric matrix and checks if it succeeds. */
void CHECKdsyevd (int n, double *M, double *eigval);

/* Calls Lapack rotine 'dgetrf' for triangular     */
/* matrix factorization and checks if it succeeds. */
void CHECKdgetrf (int n, double *M, int *ipiv);

/* Calls Lapack rotine 'dgetri' for computing the */
/* inverse matrix and checks if it succeeds.      */
void CHECKdgetri (int n, double *M, int *ipiv);

/* Allocates a block of bytes if there are     */
/* enough memory, otherwise exits the program. */
void *CHECKmalloc (unsigned int nbytes);

/* Change the size of a block of bytes if there */
/* are enough memory or exits the program.      */
void *CHECKrealloc (void *ptr1, unsigned int nbytes);

/* Opens the file named 'filename' in order to    */
/* execute a 'mode' operation and verifies error. */
FILE *CHECKfopen (const char *filename, const char *mode);

/* Verifies the returned value of a call to the    */
/* 'fclose' function to close the file 'filename'. */
void CHECKfclose (int info, const char *filename);

/* Verifies the returned value of a call to the 'fscanf' */
/* function to read data from the file named 'filename'. */
void CHECKfscanf (int info, const char *filename);

/* Verifies the returned value of a call to the 'sscanf'   */
/* function to read 'nvar' variables from an string 'str'. */
void CHECKsscanf (int info, int nvar, char *str);

/* Verifies the returned value of a call to the 'fread' function  */
/* to read 'count' blocks of data from the file named 'filename'. */
void CHECKfread (unsigned int info, unsigned int count, const char *filename);


/* ************************ Drafts ************************* */

