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
/**   File: Extern.h                                        **/
/**                                                         **/
/**   Versions: 1 - 10/10/2012                              **/
/**             2 - 08/01/2013                              **/
/**                                                         **/
/**  *****************************************************  **/
/**  Interface for some external functions from 'LAPACK'    **/
/**  and 'BLAS' packages.                                   **/
/**  *****************************************************  **/


/* For compiling with old version of Intel MKL. */
#ifdef OLD
#define dsyevd dsyevd_
#define dgetrf dgetrf_
#define dgetri dgetri_
#define dgemm dgemm_
#endif


/* Lapack rotine: computes all eigenvalues and       */
/* (optionally) all eigenvectors of a real symmetric */
/* matrix using "divide and conquer" algorithm.      */
void dsyevd (char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
	     double *work, int *lwork, int *iwork, int *liwork, int *info);

/* Lapack rotine: forms a triangular matrix factorization (trf) */
/* from a general matrix (ge) of double precision real (d).     */
void dgetrf (int *m, int *n, double *a, int *lda, int *ipiv, int *info);

/* Lapack rotine: computes the inverse matrix using the factorization */
/* (tri) from a general matrix (ge) of double precision real (d).     */
void dgetri (int *n, double *a, int *lda, int *ipiv,
	    double *work, int *lwork, int *info);

/* Blas rotine: computes a scalar-matrix-matrix product */
/* and adds the result to a scalar-matrix product.      */
void dgemm (char *transa, char *transb, int *m, int *n, int *k,
	    double *alpha, double *a, int *lda, double *b,
	    int *ldb, double *beta, double *c, int *ldc);


/* ************************ Drafts ************************* */

