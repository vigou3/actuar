/*  ===== actuar: an R package for Actuarial Science =====
 *
 *	Function to generate n pseudo random variable of a phase
 *	type distribution PH(pi,T,m)
 *
 *	Costs of this algorithm
 *	- time : O(n.m.mean)
 *	- memory : O(m^2)
 *	where 'mean' the expected number of jumps of the underlying Markov chain.
 *	'mean'	= m if PH(pi, T, m) is equivalent to an generalized Erlang distribution
 *			= 1 if PH(pi, T, m) is equivalent to a mixture of exponential distributions
 *
 *	AUTHOR: Christophe Dutang and
 *			Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          
 */
 

 // TODO
 // jumper number, mean
 // http://www.ingentaconnect.com/content/routledg/sact/2003/00002003/00000004/art00002?crawler=true
 
#include "randphasetype.h"
#include "actuar.h"

SEXP randphasetype( SEXP n, SEXP pi, SEXP T, SEXP tzero)
{
	/* extract arguments */
	PROTECT(n = coerceVector(n, INTSXP));
	int nbrand = *INTEGER(n);
	
	PROTECT(pi = coerceVector( duplicate(pi), REALSXP));
	R_len_t lenPi = GET_LENGTH(pi);
	double *valPi = REAL( pi );
	
	PROTECT(tzero = coerceVector(tzero, REALSXP));
	R_len_t lenTzero = GET_LENGTH(tzero);
	double *valTzero = REAL( tzero );
	
	SEXP dimT;
	PROTECT(dimT = allocVector(INTSXP, 2));
	PROTECT(dimT = GET_DIM(T));
	
	SEXP duplicateT = PROTECT( duplicate(T) );
	double *valT = REAL( duplicateT );
		
	/* row number and column number of T */
	R_len_t nrowT = INTEGER(dimT)[0];
	R_len_t ncolT = INTEGER(dimT)[1];
	int ncolplus1T = ncolT+1;
	int nrowplus1T = nrowT+1;
	
	/* check arguments */
	if(nrowT != ncolT || ncolT < 0 || nrowT < 0)
		error(_("T must be a non null 'nxn' matrix"));
	if(lenPi != ncolT || lenTzero != lenPi)
		error(_("vectors Pi and tzero must be compatible with size of matrix T"));		
	
	/* variable declaration */
	int i, j, isim; //loop index
	int currentState, state;
	//vector of visit number of the Markov chain
	int *nbvisit = (int *) R_alloc(ncolplus1T, sizeof(int));
	//associated probability vector of state 'state'
	double *probstate = (double *) R_alloc(ncolT, sizeof(double)); 
	//associated value vector of state 'state'
	double *valstate = (double *) R_alloc(ncolT, sizeof(double));
	int is2ptdistr = 0; //is two point distribution
	int indexmassvalue = 0; //index such that probstate[indexmassvalue]=1
	
	int nbjump;
	
	//transition matrix Qhat and its associated temporary variables
	double *Qhat = (double *) R_alloc( ncolplus1T * ncolplus1T, sizeof(double)); //transition matrix Q^
	double *exitrate = (double *) R_alloc( ncolplus1T, sizeof(double));; // exit rates
	double currentrowsum = 0.0;
		
		
	
	//result
	double *result = (double *) R_alloc(nbrand, sizeof(double)); // vector of results
	SEXP resultinR;	// result in R type
	PROTECT(resultinR = allocVector(REALSXP, nbrand)); // allocate a real vector of length x
	result = REAL( resultinR ); // plug the C pointer 'result' on the R type 'resultinR'
	
	
	/* special case where T is 1x1 matrix */
	if(ncolT == 1)
	{
		double rate = -valT[0];
		
		GetRNGstate();
		for(isim = 0; isim < nbrand; isim++)
		{
			result[isim] = exp_rand() / rate;
		}
		PutRNGstate();
		
		UNPROTECT(7);
	
		return resultinR;
	}
	
	/* initialize Qhat */
	//i row index, j column index
	for(i = 0; i < nrowplus1T; i++)
	{ 
		exitrate[i] = 0.0; //reset
		
		//for a line i, compute the ith row of intensity matrix Q
		for(j = 0; j < ncolplus1T; j++) 
		{
			if(i == 0)
				Qhat[j + i * ncolplus1T] = 0.0;
			if(i > 0 && j == 0)
				Qhat[j + i * ncolplus1T] = valTzero[i-1];	
			if(i > 0 && j > 0)	
				Qhat[j + i * ncolplus1T] = valT[i-1 + (j-1) * ncolT]; //warning: T is stored column by column
			if(i != j)
			{
				exitrate[i] += Qhat[j + i * ncolplus1T];			
				//Rprintf("%d %d %f %f\n",i,j, Qhatii, Qhat[j + i * ncolplus1T]);
			}
		}
		
		//Rprintf("%f - %f | %d \n", Qhat[i + i* ncolplus1T] ,Qhatii ,i);
		
		if(exitrate[i] + Qhat[i + i * ncolplus1T] != 0)
			error(_("T must be a sub intensity matrix"));
		
		//for a line i, compute the ith row of transition matrix Qhat
		if(i > 0)
		{					
			if(exitrate[i] == 0)
				error(_("T must be a sub intensity matrix"));						
			
			currentrowsum = 0;
			for(j = 0; j < ncolplus1T; j++) 
			{
				//divide the non diagonal term by the exitrate of state i
				//and set to 0 the diagonal term
				if(j < ncolplus1T-1 && i < ncolplus1T-1)
				{
					if(i != j)
						Qhat[j + i * ncolplus1T] /= exitrate[i];
					else
						Qhat[j + i * ncolplus1T] = 0.0;
				}	
				//to avoide numerical errors, set the last term of columns such
				//that the rowsum is equal to 1	
				if(j == ncolplus1T-1 && i < ncolplus1T-1)
				{				
					Qhat[j + i * ncolplus1T] = 1 - currentrowsum;
				}
				//special case of the last line
				if(i == ncolplus1T-1)
				{
					if(j < ncolplus1T-2)
						Qhat[j + i * ncolplus1T] /= exitrate[i];
					if(j == ncolplus1T-2)
						Qhat[j + i * ncolplus1T] = 1 - currentrowsum;
					if(j == ncolplus1T-1)
						Qhat[j + i * ncolplus1T] = 0.0;			
				}
				
				currentrowsum += Qhat[j + i * ncolplus1T];
				//Rprintf("%d %d %.30g \n",i,j, currentrowsum);		
			}
			//Rprintf("%.20g\n", currentrowsum);
		}
	}
	
	/*
	Rprintf("matrix Q^\n");
	for(i = 0; i < nrowplus1T; i++) 
	{
		for(j = 0; j < ncolplus1T; j++) 
			Rprintf("%.20g\t", Qhat[j + i* ncolplus1T] );
		
		Rprintf("\n");		
	}
	Rprintf("\n");
	Rprintf("Q^ stockage reel\n");
	for(i = 0; i < nrowplus1T*nrowplus1T; i++) 
	{
		Rprintf("%f,\t", Qhat[i] );
	}
	Rprintf("\n");*/

	
	nbjump=0;//jump number
	
	//phase type simulation algorithm by simulating the underlying markov chain
	//on {0,1,...,m}, with 0 the absorbing state and Qhat the transition matrix
	for(isim = 0; isim < nbrand; isim++)
	{
		//init the isim th variable
		result[isim] = 0.0;
		
		//init the vector nbvisit couting the number of visits of the markov chain
		for(j = 0; j < ncolplus1T; j++)
			nbvisit[j] = 0; //by definition nbvisit[j] is the visit number in state 'i'
		
		//choose the initial state in {1,...,m} according to vector Pi
		for(j = 0; j < lenPi; j++)
			probstate[j] = valPi[j];		
		//cumsum vector probstate
		for(j = 1; j < ncolT; j++)
			probstate[j] += probstate[j-1];
		//set the last mass prob value to 1	
		if(probstate[ncolT-1] != 1)
			probstate[ncolT-1] = 1;			
		//set valstate vector to {1,...,m}
		for(j = 0; j < ncolT; j++)
			valstate[j] = j+1;
			
		//check probstate defines a two point distribution
		is2ptdistr = (probstate[0] == 1) || (probstate[0] == 0);
		indexmassvalue = 0;
		for(j = 0; j < ncolT; j++)
		{
			if(probstate[j] == 1 )
			{
				indexmassvalue = j;			
				break;
			}	
			else			
				is2ptdistr *= (probstate[j] == 0);	
		}	
		/*for(i = 0; i < ncolT; i++)		
			{
				Rprintf("%f ",probstate[i] );
			}*/	
			//Rprintf("%d is2point dim %d /// %f %f\n",is2ptdistr,indexmassvalue,probstate[indexmassvalue], valstate[indexmassvalue]);
			
			
		if(is2ptdistr)
			currentState = (int) genDiscretVariable(1, probstate+indexmassvalue, valstate+indexmassvalue);
		else
			currentState = (int) genDiscretVariable(ncolT, probstate, valstate );
		
		//Rprintf("etat init %d \n",currentState);
		
		while(currentState != 0)
		{
			//increment the visit number
			nbvisit[currentState]++;
			
			//extract the probability vector associated with 'currentState'.
			//currentState is the row index, j the column index of matrix Qhat
			for(j = 0; j < ncolplus1T; j++)
			{
				if(currentState > j)
				{
					valstate[j] = j;
					probstate[j] = Qhat[j + currentState * ncolplus1T]; 
					if(j > 0)			
						probstate[j] += probstate[j-1];
				}
				if(currentState < j)
				{
					valstate[j-1] = j;
					probstate[j-1] = Qhat[j + currentState * ncolplus1T];
					if(j > 1)			
						probstate[j-1] += probstate[j-2];
				}
			}
			
			
		/*
			for(j = 0; j < ncolplus1T; j++)
			{
				Rprintf("%d ", nbvisit[j]);			
			}
			Rprintf("\n");*/
			
			//check probstate defines a two point distribution such that
			//probstate[indexmassvalue]
			is2ptdistr = (probstate[0] == 1) || (probstate[0] == 0);
			indexmassvalue = 1;
			for(j = 0; j < ncolT; j++)
			{
				if(probstate[j] == 1 )
				{
					indexmassvalue = j;			
					break;
				}	
				else			
					is2ptdistr *= (probstate[j] == 0);	
			}
			/*
			for(i = 0; i < ncolT; i++)		
			{
				Rprintf("%f ",probstate[i] );
			}	*/
			//Rprintf("%d is2point dim %d /// %f %f\n",is2ptdistr,indexmassvalue,probstate[indexmassvalue], valstate[indexmassvalue]);
			
			//change state
			if(is2ptdistr)
				currentState = (int) genDiscretVariable(1, probstate+indexmassvalue, valstate+indexmassvalue );
			else
				currentState = (int) genDiscretVariable(ncolT, probstate, valstate );
			
			nbjump++;
			//Rprintf("%d ",nbjump);
			//Rprintf("%d \n",currentState);
		}
		
		
		
		
		//simulate the visit duration of the markov process in the visited states
		//of {1,...,m} (0 is excluded since it is the absorbing state).
		for(state = 1; state < ncolplus1T; state++)
		{
			if(nbvisit[state] > 0)
				result[isim] += genErlangVariable(nbvisit[state], exitrate[state]);
		}	
		
	}
	
	
	
	Rprintf("nb moyen sauts : %f",(double) nbjump / nbrand);
	
	
	/* //test le generateur de v.a. discrete
	int nb = 0;
	int nbZERO,nbUN,nbDEUX,nbTROIS = 0;
	int rand;
	
	while(nb < 10001)
	{
		rand = (int) genDiscretVariable(ncolT, probstate, valstate );
		switch(rand)
		{	
			case 0: nbZERO++;break;
			case 1: nbUN++;break;
			case 2: nbDEUX++;break;
			case 3: nbTROIS++;break;
			default: printf("zog\n");break;
		}
		nb++;
	}	
	Rprintf("0 : %d, 1 : %d, 2 : %d, 3 : %d, nb %d\n",nbZERO,nbUN,nbDEUX,nbTROIS,nb);
	*/
	
	/* //test le generateur d'erlang
	int nb=1;
	double sum=0.0;
	double sumsq=0.0;
	double temp;
	while(nb < 100001)
	{
		temp = genErlangVariable( 3, 2.0);
		sum += temp;
		sumsq += temp*temp;
		nb++;
	}
	Rprintf("\nmoyenne %f theo %f empiriq || variance %f %f\n",3/2.0,sum/nb,3/4.0,sumsq/nb-sum/nb*sum/nb);
	*/		
			
			
		
	UNPROTECT(7);
	
	return resultinR;
}

/*	generate a random variable U of a discrete distribution in a finite state E, where 'dim' is the
*	dimension of space E, 'prob' is the cumulative mass probability vector and 'value' its 
*	corresponding value vector. We use the inverse method
*/
double genDiscretVariable(int dim, double *prob, double *value )
{
	int i; //loop index
	int argOK = 1; //pseudo boolean
	
	/*
	for(i = 0; i < dim; i++)		
	{
		Rprintf("%f ",prob[i] );
	}	
	Rprintf("\n");*/
	
	//check arguments not needed
	for(i = 0; i < dim; i++)		
	{
		//Rprintf("%d- %d\n",i,prob[i] >= 0.0);
		argOK *= (prob[i] >= 0);
		if(i > 0)
		{
			argOK *= (prob[i] > prob[i-1]);
		//	Rprintf("%d- %d\n",i,prob[i] > prob[i-1]);		
		}
		if(i > 0)
		{
			argOK *= (value[i] > value[i-1]);		
		//	Rprintf("%d- %d\n",i,value[i] > value[i-1]);	
			}
			
		if(i == dim-1)
		{
			argOK *= (prob[i] == 1.0);	
		//	Rprintf("%d- %d - %.30g\n",i,prob[i] == 1.0,prob[i] );	
			}
	}
	
	if(!argOK)
		error(_("wrong argument for genDiscretVariable function"));
	
	GetRNGstate();
	double u = unif_rand();
	PutRNGstate();
	//Rprintf("u : %f \n",u);

	//Rprintf("%d / %f / %f\t",i,prob[0],value[0]);
	
	if(u <= prob[0])
		return value[0];
	
	for(i = 1; i < dim; i++)		
	{
	//	Rprintf("%d / %f / %f \t",i,prob[i],value[i]);
		if(prob[i-1] < u && u <= prob[i])
			return value[i];
	}
	//Rprintf("\n");
}

/*	generate a random variable U of an Erlang distribution E(n,rate)
*/
double genErlangVariable(int n, double rate)
{	
	double res = 0.0;
	int i;
	GetRNGstate();
	for(i = 0; i < n; i++)
	{
		res += exp_rand() / rate;
	//	Rprintf("%d - ",i);
	}
	//Rprintf("\n");
	PutRNGstate();
	return res;
}