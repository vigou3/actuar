/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute the exponential of a matrix.
 *
 *  AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Christophe
 *  Dutang
 */

#include <R_ext/BLAS.h>

/* Matrix exponential exp(x), where x is an (n x n) matrix. Result z
 * is an (n x n) matrix. */
static void expm(double *x, int n, double *z)
{
    z = x;			/* for testing purpose ;-) */
}

/* Product x * exp(M) * y, where x is an (1 x n) vector, M is an (n x
 * n) matrix and y is an (n x 1) vector. Result z is a scalar. */
static void expmprod(double *x, double *M, double *y, int n, double *z)
{
    char *transa = "N";
    double *tmp, *eM, one = 1.0, zero = 0.0;

    /* Compute exp(M) */
    expm(M, n, eM);

    /* Product      tmp   := x     * exp(M)
     * (Dimensions: 1 x n    1 x n   n x n)
     *
     * !!! DELETE WHEN DONE !!!
     * DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC);
     * void F77_NAME(dgemm)(const char *transa, const char *transb, const int *m,
     * 		const int *n, const int *k, const double *alpha,
     * 		const double *a, const int *lda,
     * 		const double *b, const int *ldb,
     * 		const double *beta, double *c, const int *ldc);
     */
     F77_CALL(dgemm)(transa, transa, &one, &n, &n, &one,
		     x, &one, eM, &n, &zero, tmp, &one);

    /* Product      z     := tmp   * y
     * (Dimensions: 1 x 1    1 x n   n x 1)
     *
     * !!! DELETE WHEN DONE !!!
     * DDOT(N,DX,INCX,DY,INCY)
     * double F77_NAME(ddot)(const int *n, const double *dx, const int *incx,
     *			const double *dy, const int *incy);
     */
     z = F77_CALL(ddot)(&n, tmp, &one, y, &one);
}
