/*  actuar: Actuarial Functions and Heavy Tailed Distributions
 *
 *  Function to compute the recursive part of the Panjer formula
 *  to approximate the aggregate claim amount distribution of
 *  a portfolio over a period.
 *
 *  AUTHORS: Tommy Ouellet, Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "locale.h"

#define CAD5R(e)  CAR(CDR(CDR(CDR(CDR(CDR(e))))))
#define CAD6R(e)  CAR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))
#define CAD7R(e)  CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))
#define CAD8R(e)  CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))))
#define CAD9R(e)  CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e))))))))))
#define CAD10R(e) CAR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(CDR(e)))))))))))

#define INITSIZE 100		/* default size for prob. vector */

SEXP actuar_do_panjer(SEXP args)
{
    SEXP p0, p1, fs0, sfx, a, b, conv, tol, maxit, echo, sfs;
    double *fs, *fx, cumul;
    int upper, m, k, n, x = 1;
    double norm;                /* normalizing constant */
    double term;                /* constant in the (a, b, 1) case */

    /*  The length of vector fs is not known in advance. We opt for a
     *  simple scheme: allocate memory for a vector of size 'size',
     *  double the size when the vector is full. */
    int size = INITSIZE;
    fs = (double *) S_alloc(size, sizeof(double));

    /*  All values received from R are then protected. */
    PROTECT(p0 = coerceVector(CADR(args), REALSXP));
    PROTECT(p1 = coerceVector(CADDR(args), REALSXP));
    PROTECT(fs0 = coerceVector(CADDDR(args), REALSXP));
    PROTECT(sfx = coerceVector(CAD4R(args), REALSXP));
    PROTECT(a = coerceVector(CAD5R(args), REALSXP));
    PROTECT(b = coerceVector(CAD6R(args), REALSXP));
    PROTECT(conv = coerceVector(CAD7R(args), INTSXP));
    PROTECT(tol = coerceVector(CAD8R(args), REALSXP));
    PROTECT(maxit = coerceVector(CAD9R(args), INTSXP));
    PROTECT(echo = coerceVector(CAD10R(args), LGLSXP));

    /* Initialization of some variables */
    fx = REAL(sfx);             /* severity distribution */
    upper = length(sfx) - 1;    /* severity distribution support upper bound */
    fs[0] = REAL(fs0)[0];       /* value of Pr[S = 0] (computed in R) */
    cumul = REAL(fs0)[0];       /* cumulative probability computed */
    norm = 1 - REAL(a)[0] * fx[0]; /* normalizing constant */
    n = INTEGER(conv)[0];	   /* number of convolutions to do */

    /* If printing of recursions was asked for, start by printing a
     * header and the probability at 0. */
    if (LOGICAL(echo)[0])
        Rprintf("x\tPr[S = x]\tCumulative probability\n%d\t%.8g\t%.8g\n",
                0, fs[0], fs[0]);

    /* (a, b, 0) case (if p0 is NULL) */
    if (isNull(CADR(args)))
        do
        {
            /* Stop after 'maxit' recursions and issue warning. */
            if (x > INTEGER(maxit)[0])
            {
                warning(_("maximum number of recursions reached before the probability distribution was complete"));
                break;
            }

            /* If fs is too small, double its size */
            if (x >= size)
            {
                fs = (double *) S_realloc((char *) fs, size << 1, size, sizeof(double));
                size = size << 1;
            }

	    m = x;
	    if (x > upper) m = upper; /* upper bound of the sum */

            /* Compute probability up to the scaling constant */
            for (k = 1; k <= m; k++)
                fs[x] += (REAL(a)[0] + REAL(b)[0] * k / x) * fx[k] * fs[x - k];
            fs[x] = fs[x]/norm;   /* normalization */
            cumul += fs[x];       /* cumulative sum */

            if (LOGICAL(echo)[0])
                Rprintf("%d\t%.8g\t%.8g\n", x, fs[x], cumul);

            x++;
        } while (cumul < REAL(tol)[0]);
    /* (a, b, 1) case (if p0 is non-NULL) */
    else
    {
        /* In the (a, b, 1) case, the recursion formula has an
         * additional term involving f_X(x). The mathematical notation
         * assumes that f_X(x) = 0 for x > m (the maximal value of the
         * distribution). We need to treat this specifically in
         * programming, though. */
	double fxm;

        /* Constant term in the (a, b, 1) case. */
        term = (REAL(p1)[0] - (REAL(a)[0] + REAL(b)[0]) * REAL(p0)[0]);

        do
        {
            /* Stop after 'maxit' recursions and issue warning. */
            if (x > INTEGER(maxit)[0])
            {
                warning(_("maximum number of recursions reached before the probability distribution was complete"));
                break;
            }

            if (x >= size)
            {
                fs = (double *) S_realloc((char *) fs, size << 1, size, sizeof(double));
                size = size << 1;
            }

	    m = x;
	    if (x > upper)
	    {
		m = upper;	/* upper bound of the sum */
		fxm = 0.0;	/* i.e. no additional term */
	    }
	    else
		fxm = fx[m];	/* i.e. additional term */

            for (k = 1; k <= m; k++)
                fs[x] += (REAL(a)[0] + REAL(b)[0] * k / x) * fx[k] * fs[x - k];
            fs[x] = (fs[x] + fxm * term) / norm;
            cumul += fs[x];

            if (LOGICAL(echo)[0])
                Rprintf("%d\t%.8g\t%.8g\n", x, fs[x], cumul);

            x++;
        } while (cumul < REAL(tol)[0]);
    }

    /* If needed, convolve the distribution obtained above with itself
     * using a very simple direct technique. Since we want to
     * continue storing the distribution in array 'fs', we need to
     * copy the vector in an auxiliary array at each convolution. */
    if (n)
    {
	int i, j, ox;
	double *ofs;		/* auxiliary array */

	/* Resize 'fs' to its final size after 'n' convolutions. Each
	 * convolution increases the length from 'x' to '2 * x - 1'. */
	fs = (double *) S_realloc((char *) fs, (1 << n) * (x - 1) + 1, size, sizeof(double));

	/* Allocate enough memory in the auxiliary array for the 'n'
	 * convolutions. This is just slightly over half the final
	 * size of 'fs'. */
	ofs = (double *) S_alloc((1 << (n - 1)) * (x - 1) + 1, sizeof(double));

	for (k = 0; k < n; k++)
	{
	    memcpy(ofs, fs, x * sizeof(double)); /* keep previous array */
	    ox = x;		/* previous array length */
	    x = (x << 1) - 1;	/* new array length */
	    for(i = 0; i < x; i++)
		fs[i] = 0.0;
	    for(i = 0; i < ox; i++)
		for(j = 0; j < ox; j++)
		    fs[i + j] += ofs[i] * ofs[j];
	}
    }

    /*  Copy the values of fs to a SEXP which will be returned to R. */
    PROTECT(sfs = allocVector(REALSXP, x));
    memcpy(REAL(sfs), fs, x * sizeof(double));

    UNPROTECT(11);
    return(sfs);
}
