/**  *****************************************************  **/
/**             ** Phonon Vibration Analysis **             **/
/**                                                         **/
/**                    **  Version 1  **                    **/
/**                                                         **/
/**                         IF/USP                          **/
/**                                                         **/
/**   Advisor: Prof. Dr. Alexandre Reily Rocha              **/
/**                                                         **/
/**   Author: Pedro Brandimarte Mendonca                    **/
/**                                                         **/
/**   File: Utils.h                                         **/
/**                                                         **/
/**  *****************************************************  **/
/**  Interface for some math-physics definitions and        **/
/**  matrix (vector) algebra utilities.                     **/
/**  *****************************************************  **/


/**  ******************** Definitions ********************  **/

/* Constants (from CODATA - 2012). */
#define hbar 6.58211928e-16 /* eV*s */
#define amu2kg 1.660538921e-27 /* kg */
#define eV2joule 1.602176565e-19 /* eV -> joule */
#define bohr2ang 0.52917721092 /* angstrom */
#define rydberg2eV 13.60569253 /* Ry*hc -> eV */

/* C matrix indexation (row-major order) */
/* #define idx(i, j, ncol) ((i) * (ncol) + j) */

/* Fortran matrix indexation (column-major order) */
#define idx(i, j, nrow) (i + (j) * (nrow))
#define idx3d(i, j, k, nrow, ncol) (i + (j) * (nrow) + (k) * (nrow) * (ncol))


/**  *********************** Types ***********************  **/

/* Complex number. */
typedef struct {
   double re; /* real part */
   double im; /* imaginary part */
} cdouble;


/**  *********** Matrix and Vectors Utilities ************  **/

/* Allocates and initialize a vector 'V[n]' of integers. */
void *UTILintVector (int n);

/* Allocates and initializes a vector 'V[n]' of doubles. */
void *UTILdoubleVector (int n);

/* Resets a vector 'V[n]' of doubles with 0.0's. */
void UTILresetDoubleVector (int n, double *V);

/* Copies the elements form a vector 'Orig[n]' to 'Dest[n]'. */
void UTILcopyVector (double *Dest, double *Orig, int n);

/* Checks if the matrix 'M[n][n]' is symetric. */
void UTILcheckSym (char *name, double *M, int n, double lim);

/* Computes the Euclidean norm of a double precision vector 'V[n]'. */
double UTILnorm (int n, double *V);

/* ************************ Drafts ************************* */

