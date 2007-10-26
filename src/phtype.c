/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute density, cumulative distribution and moment
 *  generating functions, raw moments and to simulate random variates
 *  for Phase-type distributions. See ../R/PhaseType.R for details.
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 */

#include <R.h>
#include <Rmath.h>
#include "actuar.h"
#include "locale.h"
#include "dpq.h"

double dphtype(double x, double *pi, double *T, int m, int give_log)
{
    /*  Density function is
     *
     *	pi      * exp(x * T) * t
     *  (1 x m)   (m x m)      (m x 1)
     *
     *  with t = -T * e and e is a 1-vector.
     */

    if (!R_FINITE(x) || x < 0.0)
	return R_D__0;

    int i, j, jm;
    double *t, *tmp;

    /* Build vector t (equal to minus the row sums of matrix T) and
     * matrix tmp = x * T. */
    t = (double *) S_alloc(m, sizeof(double)); /* initialized to 0 */
    tmp = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m; i++)
	for (j = 0; j < m; j++)
	{
	    jm = j * m;
	    t[i] -= T[i + jm];
	    tmp[i + jm] = x * T[i + jm];
	}

    return R_D_val(expmprod(pi, tmp, t, m));
}

double pphtype(double q, double *pi, double *T, int m, int lower_tail,
	       int log_p)
{
    /*  Cumulative distribution function is
     *
     *	1 - pi      * exp(q * T) * e
     *      (1 x m)   (m x m)      (m x 1)
     *
     *  where e is a 1-vector.
     */

    if (q <= 0)
	return R_DT_0;

    int i;
    double *e, *tmp;

    /* Create the 1-vector and multiply each element of T by q. */
    e = (double *) R_alloc(m, sizeof(double));
    for (i = 0; i < m; i++)
	e[i] = 1;
    tmp = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m * m; i++)
	tmp[i] = q * T[i];

    return R_DT_Cval(expmprod(pi, tmp, e, m));
}

double rphtype(double *pi, double *T, int m)
{
    /* constants */
	int i, j; 
	int currentState, state;
	int is2ptdistr = 0; /* is two point distribution */
	int indexmassvalue = 0; /* index such that probstate[indexmassvalue]=1 */
	int nbjump;
	double currentrowsum = 0.0;
	double result ;
	int mplus1 = m+1;
	
	/* arrays */
	double *t = (double *) S_alloc(m, sizeof(double)); /* initialized to 0 */
	int *nbvisit = (int *) R_alloc(mplus1, sizeof(int)); /*vector of visit number of the Markov chain*/
	double *probstate = (double *) R_alloc(m, sizeof(double)); /*associated probability vector of state 'state'*/
	double *valstate = (double *) R_alloc(m, sizeof(double)); /*associated value vector of state 'state'*/
	/* transition matrix Qhat (a (m+1)x(m+1) matrix) and its associated temporary variables */
	double *Qhat = (double *) R_alloc( mplus1 * mplus1, sizeof(double)); /*transition matrix Q^*/
	double *exitrate = (double *) R_alloc( mplus1, sizeof(double));; /* exit rates */
		
	
	/* Build vector t (equal to minus the row sums of matrix T) */
	for (i = 0; i < m; i++)
		for (j = 0; j < m; j++)
			t[i] -= T[i + j * m];
	
	
	/* special case where T is 1x1 matrix */
	if(m == 1)
	{
		double rate = -T[0];
		GetRNGstate();
		result = exp_rand() / rate;
		PutRNGstate();
		return result;
	}
	
	/* initialize Qhat :
	 * i row index, j column index */
	for(i = 0; i < mplus1; i++)
	{ 
		exitrate[i] = 0.0; //reset
		
		/* compute the ith row of intensity matrix Q */
		for(j = 0; j < mplus1; j++) 
		{
			if(i == 0)
				Qhat[j + i * mplus1] = 0.0;
			if(i > 0 && j == 0)
				Qhat[j + i * mplus1] = t[i-1];	
			if(i > 0 && j > 0)	
				Qhat[j + i * mplus1] = T[i-1 + (j-1) * m]; /* warning: T is stored column by column */
			if(i != j)
				exitrate[i] += Qhat[j + i * mplus1];						
		}
		
		if(exitrate[i] + Qhat[i + i * mplus1] != 0)
			error(_("T must be a sub intensity matrix"));
		
		/* compute the ith row of transition matrix Qhat */
		if(i > 0)
		{					
			if(exitrate[i] == 0)
				error(_("T must be a sub intensity matrix"));						
			
			currentrowsum = 0;
			for(j = 0; j < mplus1; j++) 
			{
				/* divide the non diagonal term by the exitrate of 
				   state i and set to 0 the diagonal term */
				if(j < mplus1-1 && i < mplus1-1)
				{
					if(i != j)
						Qhat[j + i * mplus1] /= exitrate[i];
					else
						Qhat[j + i * mplus1] = 0.0;
				}	
				/* to avoide numerical errors, set the last term of each
				   column such that the rowsum is equal to 1 */	
				if(j == mplus1-1 && i < mplus1-1)
				{				
					Qhat[j + i * mplus1] = 1 - currentrowsum;
				}
				/* special case of the last line */
				if(i == mplus1-1)
				{
					if(j < mplus1-2)
						Qhat[j + i * mplus1] /= exitrate[i];
					if(j == mplus1-2)
						Qhat[j + i * mplus1] = 1 - currentrowsum;
					if(j == mplus1-1)
						Qhat[j + i * mplus1] = 0.0;			
				}
				
				currentrowsum += Qhat[j + i * mplus1];
			}
		}
	}
	
	nbjump=0;
	
	/* phase type simulation algorithm simulates the underlying 
	 * markov chain on space {0,1,...,m}, with 0 the absorbing 
	 * state and Qhat the transition matrix */
	
	result = 0.0;
		
	/* init the vector nbvisit couting the number of visits of 
	 * the markov chain. nbvisit[j] is defined as the visit 
	 * number in state 'j' */
	for(j = 0; j < mplus1; j++)
		nbvisit[j] = 0; 
		
	/* choose the initial state in {1,...,m} according to vector Pi */
	/* init the probstate vector */
	for(j = 0; j < m; j++)
		probstate[j] = pi[j];		
	/* cumsum vector probstate */
	for(j = 1; j < m; j++)
		probstate[j] += probstate[j-1];
	/* set the last mass prob value to 1 */	
	if(probstate[m-1] != 1)
		probstate[m-1] = 1;			
	/* set valstate vector to {1,...,m} */
	for(j = 0; j < m; j++)
		valstate[j] = j+1;
			
	/* check if probstate is a two point distribution and set
	 * indexmassvalue such that probstate[indexmassvalue] = 1. 
	 * typically, two point distribution appears with Erlang 
	 * distribution */
	is2ptdistr = (probstate[0] == 1) || (probstate[0] == 0);
	indexmassvalue = 0;
	for(j = 0; j < m; j++)
	{
		if(probstate[j] == 1 )
		{
			indexmassvalue = j;			
			break;
		}	
		else			
			is2ptdistr *= (probstate[j] == 0);	
	}	
	/* simulate the initial state */
	if(is2ptdistr)
		currentState = (int) genDiscretVariable(1, probstate+indexmassvalue, valstate+indexmassvalue);
	else
		currentState = (int) genDiscretVariable(m, probstate, valstate );
	
	/* simulate the Markov chain */		
	while(currentState != 0)
	{
		/* increment the visit number */
		nbvisit[currentState]++;
			
		/* extract the probability vector associated with 'currentState'.
		 * currentState is the row index, j the column index of transition
		 * matrix Qhat */
		for(j = 0; j < mplus1; j++)
		{
			if(currentState > j)
			{
				valstate[j] = j;
				probstate[j] = Qhat[j + currentState * mplus1]; 
				if(j > 0)			
					probstate[j] += probstate[j-1];
			}
			if(currentState < j)
			{
				valstate[j-1] = j;
				probstate[j-1] = Qhat[j + currentState * mplus1];
				if(j > 1)			
					probstate[j-1] += probstate[j-2];
			}
		}
			
			
		/* check if probstate is a two point distribution */
		is2ptdistr = (probstate[0] == 1) || (probstate[0] == 0);
		indexmassvalue = 1;
		for(j = 0; j < m; j++)
		{
			if(probstate[j] == 1 )
			{
				indexmassvalue = j;			
				break;
			}	
			else			
				is2ptdistr *= (probstate[j] == 0);	
		}
			
		/* change state */
		if(is2ptdistr)
			currentState = (int) genDiscretVariable(1, probstate+indexmassvalue, valstate+indexmassvalue );
		else
			currentState = (int) genDiscretVariable(m, probstate, valstate );
		
		nbjump++;
	}
		
		
	/* simulate the visit duration of the markov process in the 
	 * visited states of space {1,...,m} (0 is excluded since it 
	 * is the absorbing state). */
	for(state = 1; state < mplus1; state++)
	{
		if(nbvisit[state] > 0)
			result += genErlangVariable(nbvisit[state], exitrate[state]);
	}	
	
	return result;
}


double mphtype(double order, double *pi, double *T, int m, int give_log)
{
    /*  Raw moment is
     *
     *	order!  * pi      * (-T)^(-order) * e
     *  (1 x 1)   (1 x m)   (m x m)         (m x 1)
     *
     * where e is a 1-vector. Below, the moment is computed as
     * (-1)^order * order! * sum(pi * T^(-order))
     */

    if (order < 0.0 || (int) order != order)
	return R_NaN;

    int i, j;
    double tmp = 0.0, *Tpow;

    /* Compute the power of T */
    Tpow = (double *) R_alloc(m * m, sizeof(double));
    matpow(T, m, (int) -order, Tpow);

    /* Compute vector tmp = sum(pi * Tpow) */
    for (i = 0; i < m; i++)
	for (j = 0; j < m; j++)
	    tmp += pi[j] * Tpow[i * m + j];

    /* Multiply by -1 if order is odd */
    return R_D_val((int) order % 2 ?
		   -gammafn(order + 1.0) * tmp :
		   gammafn(order + 1.0) * tmp);
}

double mgfphtype(double x, double *pi, double *T, int m, int give_log)
{
    /*  Moment generating function is
     *
     *	pi      * (-x * I - T)^(-1) * t
     *  (1 x m)   (m x m)             (m x 1)
     *
     *  with t = -T * e, e a 1-vector and I the identity matrix.
     */

    if (!R_FINITE(x) || x < 0.0)
	return R_D__0;

    int i, j, jm;
    double z = 0.0, *t, *tmp1, *tmp2;

    /* Build vector t (equal to minux the row sums of matrix T) and
     * matrix tmp1 = -x * I - T. */
    t = (double *) S_alloc(m, sizeof(double)); /* initialized to 0 */
    tmp1 = (double *) R_alloc(m * m, sizeof(double));
    for (i = 0; i < m; i++)
	for (j = 0; j < m; j++)
	{
	    jm = j * m;
	    t[i] -= T[i + jm];
	    tmp1[i + jm] = (i == j) ? -x - T[i + jm] : -T[i + jm];
	}

    /* Compute tmp2 = tmp1^(-1) * t */
    tmp2 = (double *) R_alloc(m, sizeof(double));
    solve(tmp1, t, m, 1, tmp2);

    /* Compute z = pi * tmp2*/
    for (i = 0; i < m; i++)
	z += pi[i] * tmp2[i];

    return R_D_val(z);
}
