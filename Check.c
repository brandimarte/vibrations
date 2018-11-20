/**  *****************************************************  **/
/**             ** Phonon Vibration Analysis **             **/
/**                                                         **/
/**                    **  Version 2  **                    **/
/**                                                         **/
/**   By: Pedro Brandimarte (brandimarte@gmail.com) and     **/
/**       Alexandre Reily Rocha (reilya@ift.unesp.br)       **/
/**                                                         **/
/**  *****************************************************  **/
/**  Implementation of alternative calls of well known      **/
/**  functions from 'stdio.h' and 'stdlib.h', and also      **/
/**  functions from 'LAPACK' package, to avoid code         **/
/**  repetition when checking errors.                       **/
/**  *****************************************************  **/

#include <stdio.h>
#include <stdlib.h>
#include "Extern.h"
#include "Check.h"


/* ********************************************************* */
/* Lapack rotine: computes all eigenvalues and eigenvectors  */
/* of a real symmetric matrix 'M' using "divide and conquer" */
/* algorithm. Notice that in this program 'N = M = LDA = n'. */
void CHECKdsyevd (int n, double *M, double *eigval)
{
   int lwork, liwork, li, info = 0;
   double l;
   int *iwork = NULL;
   double *work = NULL;

   /* Workspace query: calculates the optimal sizes of 'work' and 'iwork'. */
   lwork = liwork = -1;
   dsyevd ("V", "U", &n, M, &n, eigval, &l, &lwork, &li, &liwork, &info);
   lwork = (int) l;
   liwork = li;
   iwork = CHECKmalloc (liwork * sizeof (int));
   work = CHECKmalloc (lwork * sizeof (double));

   /* Computes eigenvalues and eigenvectors. */
   dsyevd ("V", "U", &n, M, &n, eigval, work, &lwork, iwork, &liwork, &info);
   if (info < 0) {
      fprintf (stderr, "\n In lapack dsyevd: \n");
      fprintf (stderr, " The %d-th argument had an illegal value\n", -info);
      exit (EXIT_FAILURE);
   }
   if (info > 0) {
      fprintf (stderr, "\n In lapack dsyevd: \n");
      fprintf (stderr, " The algorithm has failed to compute an eigenvalue\n");
      fprintf (stderr, " while working on the submatrix lying in rows and \n");
      fprintf (stderr, " columns %d through %d!\n\n", info/(n + 1), info%(n + 1));
      exit (EXIT_FAILURE);
   }

   /* Frees memory. */
   free (work);
   free (iwork);

} /* CHECKdsyevd */


/* ********************************************************** */
/* Lapack rotine: forms a triangular matrix factorization     */
/* (trf) from a general matrix (ge) of double precision real  */
/* (d). Notice that in this program 'N = M = LDA = n'.        */
void CHECKdgetrf (int n, double *M, int *ipiv)
{
   int info = 0;

   dgetrf (&n, &n, M, &n, ipiv, &info);
   if (info < 0) {
      fprintf (stderr, "\n In lapack dgetrf: \n");
      fprintf (stderr, "\n The %d-th argument had an illegal value.\n\n", info);
      exit (EXIT_FAILURE);
   }
   else if (info > 0) {
      fprintf (stderr, "\n In lapack dgetrf: \n");
      fprintf (stderr, " U(%d,%d) is exactly zero. The factorization has\n", info, info);
      fprintf (stderr, " been completed, but the factor U is exactly\n");
      fprintf (stderr, " singular, and division by zero will occur if\n");
      fprintf (stderr, " it is used to solve a system of equations.\n\n");
   }

} /* CHECKdgetrf */


/* ********************************************************** */
/* Lapack rotine: computes the inverse matrix using the       */
/* factorization (tri) from a general matrix (ge) of double   */
/* precision real (d). Notice that in this program            */
/* 'N = LDA = lwork = n'.                                     */
void CHECKdgetri (int n, double *M, int *ipiv)
{
   int info = 0;
   double *work = CHECKmalloc (n * sizeof (double));

   dgetri (&n, M, &n, ipiv, work, &n, &info);
   if (info < 0) {
      fprintf (stderr, "\n In lapack dgetri: \n");
      fprintf (stderr, " The %d-th argument had an illegal value\n", info);
      exit (EXIT_FAILURE);
   }
   else if (info > 0) {
      fprintf (stderr, "\n In lapack dgetri: \n");
      fprintf (stderr, " U(%d,%d) is exactly zero. The matrix is singular\n", info, info);
      fprintf (stderr, " and its inverse could not be computed.\n\n");
      exit (EXIT_FAILURE);
   }

   /* Frees memory. */
   free (work);

} /* CHECKdgetri */


/* ********************************************************* */
/* Allocates a block of bytes if there are enough memory.    */
/* Otherwise returns an error message and exits the program. */
void *CHECKmalloc (unsigned int nbytes)
{
   void *ptr;
   ptr = malloc (nbytes);

   if (ptr == NULL) {    
      fprintf (stderr, "\n\n Insufficient memory.\n\n");
      exit(EXIT_FAILURE);
   }

   return ptr;

} /* CHECKmalloc */


/* ********************************************************* */
/* Reallocates a block of bytes pointed by 'ptr' and change  */
/* its size to 'nbytes' if there are enough memory.          */
/* Otherwise returns an error message and exits the program. */
void *CHECKrealloc (void *ptr1, unsigned int nbytes)
{
   void *ptr2;
   ptr2 = realloc (ptr1, nbytes);

   if (ptr2 == NULL) {    
      fprintf (stderr, "\n\n Insufficient memory.\n\n");
      exit (EXIT_FAILURE);
   }

   return ptr2;

} /* CHECKrealloc */


/* ********************************************************* */
/* Opens the file named 'filename' in order to execute an    */
/* operation specified by 'mode' (operations of the 'fopen'  */
/* function from 'stdio.h') and verifies error.              */
FILE *CHECKfopen (const char *filename, const char *mode)
{
   FILE *pfile;

   pfile = fopen (filename, mode);
   if (pfile == NULL) {
      fprintf (stderr, "\n\n Error: Unable to open the file '%s'!\n\n", filename);
      exit (EXIT_FAILURE);
   }

   return pfile;

} /* CHECKfopen */


/* ********************************************************* */
/* Receives the integer 'info' returned from a call to the   */
/* 'fclose' function to close a File stream from the file    */
/* named 'filename' and verifies error.                      */
void CHECKfclose (int info, const char *filename)
{

   if (info != 0) {
      fprintf (stderr, "\n\n Error: Unable to close the file '%s'!\n\n", filename);
      exit (EXIT_FAILURE);
   }

} /* CHECKfclose */


/* ********************************************************* */
/* Receives the integer 'info' returned from a call to the   */
/* 'fscanf' function to read data from the file named        */
/* 'filename' and checks if 'info = EOF', which indicates    */
/* that an input failure happened before any data could be   */
/* successfully read.                                        */
void CHECKfscanf (int info, const char *filename)
{
   if (info == EOF) {
      fprintf (stderr, "\n\n Error: Unable to read the file '%s'!\n\n", filename);
      exit (EXIT_FAILURE);
   }
   
} /* CHECKfscanf */


/* ********************************************************* */
/* Receives the integer 'info' returned from a call to the   */
/* 'sscanf' function to read data from an string 'str' and   */
/* checks if 'info' differs from the number of variables     */
/* 'nvar' that should be read.                               */
void CHECKsscanf (int info, int nvar, char *str)
{
   if (info != nvar) {
      fprintf (stderr, "\n\n Error: Problem reading the string '%s'!\n\n", str);
      exit (EXIT_FAILURE);
   }
   
} /* CHECKsscanf */


/* ********************************************************* */
/* Receives the integer 'info' returned from a call to the   */
/* 'fread' function to read 'count' blocks of data from the  */
/* file named 'filename' and checks if 'info' differs from   */
/* the 'count' parameter, which indicates either an error    */
/* ocurred or the 'End Of File' was reached.                 */
void CHECKfread (unsigned int info, unsigned int count, const char *filename)
{
   if (info != count) {
      fprintf (stderr, "\n\n Error: Read something strange at file '%s'!\n\n", filename);
      /* exit (EXIT_FAILURE); */
   }
   
} /* CHECKfread */


/* ************************ Drafts ************************* */

