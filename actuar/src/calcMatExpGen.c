/*  ===== actuar: an R package for Actuarial Science =====
 *
 *	Function to compute the matrix exponential 
 *
 *		u %*% e^( x[k] * T ) %*% v, for k in 1:length(x)	(1)
 *
 *	where u,v are vectors of R^m, T a non null m x m real matrix
 *	x a vector of R^n.
 *
 *	We use an algorithm base on Octace's expm optimized for the 
 *	equation (1).
 *
 *	Costs of this algorithm
 *	- time : O(n.m^3)
 *	- memory : O(m^2)
 *
 *	AUTHOR: Christophe Dutang and
 *			Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          
 */
 
//see R_ext/BLAS.h for fortran functions
 
 
#include "calcMatExpGen.h"
#include "actuar.h"

// compute u %*% e^( x[k] * T ) %*% v, for k in 1:length(x)
SEXP calcMatExp2 (SEXP x, SEXP u, SEXP T, SEXP v)
{
	R_len_t lenX, lenV, lenU; //length of x and v
	SEXP dimT; //size of square matrix T
	
	//extract argument sizes
	PROTECT(x = coerceVector(x, REALSXP));
	PROTECT(u = coerceVector(u, REALSXP));
	PROTECT(v = coerceVector(v, REALSXP));
	lenX = length(x);
	lenV = length(v);
	lenU = length(u);
	PROTECT(dimT = allocVector(INTSXP, 2));
	PROTECT(dimT = getAttrib(T, R_DimSymbol));
	
	/* row number and column number of T */
	R_len_t nrowT = INTEGER(dimT)[0];
	R_len_t ncolT = INTEGER(dimT)[1];
	int ncolplus1T = ncolT+1;
	int ncolTsqr = ncolT * ncolT;
	
	
	/* check arguments */
	if(nrowT != ncolT || ncolT < 0 || nrowT < 0)
		error(_("T must be a non null 'nxn' matrix"));
	if(lenV != lenU || lenV != ncolT)
		error(_("vectors U and V must be compatible with size of matrix T"));
				
		
	/* temporary working variables */
	SEXP duplicateT ; //initialized later 
	double *valT ; //initialized later
	double *valX = REAL( x );
	double *valU = REAL( u );
	double *valV = REAL( v );	
		
	/* declaration */	
	int i, j, k; //loop index
	int traceshiftT = 0;
	int iloperm, ihiperm; //arguments for dgebal permutation
	int iloscal, ihiscal; //arguments for dgebal scaling	
	double * perm = (double *) R_alloc(ncolT, sizeof(double)); //permutation array
	double * scale = (double *) R_alloc(ncolT, sizeof(double)); //scale array
	int exitcode; //exitcode used for different fortran routine
	double * work = (double *) R_alloc(ncolTsqr, sizeof(double)); //workspace array
	double infnormT; //infinite norm of matrix T
	int sqrpowscal ; //the square power used for scaling T
	double * npp = (double *) R_alloc(ncolTsqr, sizeof(double)); //numerator power Pade'
	double * dpp = (double *) R_alloc(ncolTsqr, sizeof(double)); //denominator power Pade'
	double minus1powj = -1; // (-1)^j 
	double one = 1.0; //useful for fortran routine dgemm
	double zero = 0.0;
	int * pivot = (int *) R_alloc(ncolT, sizeof(int)); // pivot vector
	int * invperm = (int *) R_alloc(ncolT, sizeof(int)); // inverse permutation vector
	
	double * result = (double *) R_alloc(lenX, sizeof(double)); // vector of results
	SEXP resultinR;	// result in R type
	PROTECT(resultinR = allocVector(REALSXP, lenX)); // allocate a real vector of length x
	result = REAL( resultinR ); // plug the C pointer 'result' on the R type 'resultinR'
	
	//plot
	/*
	Rprintf("matrix T\n");
	for(i = 0; i < nrowT; i++) 
	{
		for(j = 0; j < ncolT; j++) 
			Rprintf("%f\t", valT[i + j* ncolT] );
		
		Rprintf("\n");		
	}
	Rprintf("\n");
	*/
		
	
	/* init result */
	for(k = 0; k < lenX; k++) 
		result[k] = 0.0;
	
	/* special case where T is 1x1 matrix */
	if(ncolT == 1)
	{
		// set T to its initial value
		duplicateT = PROTECT( duplicate(T) );
		valT = REAL( duplicateT );
		
		//compute u * e^( x[k] * T ) * v
		for(k = 0; k < lenX; k++)
			result[k] = valU[0] * exp(valX[k] * valT[0]) * valV[0];
			
		// clean up			
		UNPROTECT(7);			
		
		return resultinR;
	}
	
	/* compute u %*% e^( x[k] * T ) %*% v, for k in 1:length(x) */
	for(k = 0; k < lenX; k++)
	{
		/* set T to its initial value */
		duplicateT = PROTECT( duplicate(T) );
		valT = REAL( duplicateT );
		
		/* multiply matrix T by x[k] */
		for(i = 0; i < nrowT; i++) 
			for(j = 0; j < ncolT; j++) 
				valT[i + j* ncolT]  *= valX[k];
		
	
		/* step 1 of preconditioning : trace normalisation */
		// compute the average trace of matrix T 
		
		for(i = 0; i < ncolT; i++) 
			traceshiftT += valT[i * ncolplus1T];
			
		traceshiftT /= ncolT;
		
		if(traceshiftT > 0)
		{
			for(i = 0; i < ncolT; i++) 
				valT[i * ncolplus1T] -= traceshiftT;
		}
		
	
	
		/* step 2 of preconditioning : balancing */
		/* call fortran routine dgebal with arguments
			JOB     (input) CHARACTER*1
			Specifies the operations to be performed on A:
			=  'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0 for
			i = 1,...,N; = 'P':  permute only;
			= 'S':  scale only;
			= 'B':  both permute and scale.			
			N       (input) INTEGER
			The order of the matrix A.  N >= 0.
			A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
			ILO     (output) INTEGER
			IHI     (output)  INTEGER ILO and IHI are set to integers such
			that on exit A(i,j) = 0 if i > j and j =  1,...,ILO-1  or  I  =
			IHI+1,...,N.  
			SCALE   (output) DOUBLE PRECISION array, dimension (N)
			Details  of  the permutations and scaling factors applied to A.
			If P(j) is the index of the row and  column  interchanged  with
			row  and column j and D(j) is the scaling factor applied to row
			and column j, then SCALE(j) = P(j)    for  j  =  1,...,ILO-1  =
			D(j)    for j = ILO,...,IHI = P(j)    for j = IHI+1,...,N.  The
			order in which the interchanges are made is N to IHI+1, then  1
			to ILO-1.
			INFO    (output) INTEGER, =0:  successful exit; < 0:  if INFO = -i, 
			the i-th argument had an illegal value.
			*/
			
		F77_CALL(dgebal) ("P", &ncolT, valT, &ncolT, &iloperm, &ihiperm, perm, &exitcode);
		if(exitcode) 
			error(_("LAPACK routine dgebal returned info code %d"), exitcode);
		
		F77_CALL(dgebal) ("S", &ncolT, valT, &ncolT, &iloscal, &ihiscal, scale, &exitcode);
		if(exitcode) 
			error(_("LAPACK routine dgebal returned info code %d"), exitcode);			
		
		/* step 3 of preconditioning : scaling */			
		/* call fortran routine dlange to compute the infinite norm with arguments 
			NORM    (input) CHARACTER*1
			Specifies  the  value  to  be  returned  in DLANGE as described
			above.
			M       (input) INTEGER
			The number of rows of the matrix A.  M  >=  0.   When  M  =  0,
			DLANGE is set to zero.
			N       (input) INTEGER
			The  number  of  columns of the matrix A.  N >= 0.  When N = 0,
			DLANGE is set to zero.
			A       (input) DOUBLE PRECISION array, dimension (LDA,N)
			The m by n matrix A.
			LDA     (input) INTEGER
			The leading dimension of the array A.  LDA >= max(M,1).
			WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
			where LWORK >= M when NORM = 'I'; otherwise, WORK is not refer-
			enced.
			*/		
				
		infnormT = F77_CALL(dlange) ("I", &nrowT, &ncolT, valT, &ncolT, work);
		if(infnormT > 0) 
			sqrpowscal = imax2( (int) 1 + log(infnormT)/log(2.0) , 0 );
		else
			sqrpowscal = 0;				
		
		if(sqrpowscal > 0)
		{
			double scalefactor = R_pow_di(2 , sqrpowscal);
			for(i = 0; i < ncolTsqr; i++)
				valT[i] /= scalefactor;
		}	
		
		/* Pade' approximation (p=q=8) : T^8, T^7, T^6, ..., T^1*/			
		// init npp and dpp
		
		for(i = 0; i < ncolTsqr; i++)
		{
			npp[i] = 0.0;
			dpp[i] = 0.0;
		}
		for(j = 7; j>=0; j--)
		{
			/* call fortran routine dgemm with arguments */
			/*
			 TRANSA -  CHARACTER*1.  On entry, TRANSA specifies the form of op( A )
			 to be used in the matrix multiplication as follows:
			 if = 'N' or 'n',  op( A ) = A. if = 'T' or 't',  op( A ) = A'.
			 if = 'C' or 'c',  op( A ) = A'.
			 TRANSB  -   CHARACTER*1.  On entry, TRANSB specifies the form of op( B )
			 to be used in the matrix multiplication as TRANSA
			 M      -  INTEGER. On entry, M specifies  the number  of rows  of the  
			 matrix op( A ) and of the matrix C. M >=0.
			 N      -  INTEGER. On entry, N specifies the number  of columns of the 
			 matrix op( B ) and the number of columns of the matrix C. N >=0.
			 K      -  INTEGER. On entry, K specifies  the number of columns of the 
			 matrix op( A ) and the number of rows of the matrix op( B ). K >=0.
			 ALPHA  -  DOUBLE PRECISION.
			 On entry, ALPHA specifies the scalar alpha.  
			 A      -  DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
			 k when TRANSA = 'N' or 'n', and is m otherwise. Before entry with TRANSA 
			 = 'N' or 'n', the leading m by k part of the array A must contain the matrix
			 A, otherwise the leading k by m  part of the array A must contain the matrix A.
			 LDA    -  INTEGER. On  entry, LDA specifies the first dimension of A as 
			 declared in the calling program. When TRANSA = 'N' or 'n' then  
			 LDA must be at least max( 1, m ), otherwise  LDA must be at least max( 1, k ).  
			 B      -  DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
			 n  when  TRANSB = 'N' or 'n', and  is k otherwise. Before entry  with TRANSB  
			 = 'N' or 'n', the leading k by n part of the array B must contain the matrix  
			 B, otherwise the leading n by k part of the array B must contain the matrix B.
			 Unchanged on exit.
			 LDB    -  INTEGER. On entry, LDB specifies the first dimension of B as declared  
			 in the calling program. When  TRANSB = 'N' or 'n' then LDB must be at least  
			 max( 1, k ), otherwise LDB must be at least max( 1, n ).  
			 BETA   -  DOUBLE PRECISION. On  entry, BETA specifies the scalar beta. When BETA  
			 is supplied as zero then C need not be set on input.
			 C      -  DOUBLE PRECISION array of DIMENSION ( LDC, n ). Before  entry, the 
			 leading m by n  part of the array  C must contain the matrix  C, except when 
			 beta is zero, in which case C need not be set on entry. On exit, the array C is
			 overwritten by the m by n  matrix ( alpha*op( A  )*op(  B  )  + beta*C ).
			 LDC    -  INTEGER. On  entry, LDC specifies the first dimension of C as declared in
			 the  calling program. LDC must be at least max( 1,m ). 
			 */
			 
			// npp = m * npp + padec88[j] * m 
			F77_CALL(dgemm) ("N", "N", &ncolT, &ncolT, &ncolT, &one, valT, &ncolT, npp, &ncolT, &zero, work, &ncolT);
			
			for(i = 0; i < ncolTsqr; i++)
				npp[i] = work[i] + padec88[j] * valT[i];
			
			// dpp = m * dpp + (-1)^j * padec88[j] * m 
			F77_CALL(dgemm) ("N", "N", &ncolT, &ncolT, &ncolT, &one, valT, &ncolT, dpp, &ncolT, &zero, work, &ncolT);
			
			for(i = 0; i < ncolTsqr; i++)
				dpp[i] = work[i] + minus1powj * padec88[j] * valT[i];
			
			minus1powj *= -1;	
		}
		
		// T^0 
		for(i = 0; i < ncolTsqr; i++) 
			dpp[i] *= -1.0;
			
		for(j = 0; j < ncolT; j++)
		{
			npp[j * ncolplus1T] += 1.0;
			dpp[j * ncolplus1T] += 1.0;		
		}									
		
		/* Compute pade approximation = (dpp)^-1 * npp. */																		
		/* call fortran routine dgetrf  with arguments
			M       (input) INTEGER, The number of rows of the matrix A.  M >= 0.
			N       (input) INTEGER, The number of columns of the matrix A.  N >= 0.
			A       (input/output) DOUBLE PRECISION array, dimension (LDA,N),  On entry, 
			the M-by-N matrix to be factored.  On exit, the factors L and U from 
			the factorization A = P*L*U; the unit diagonal elements of L are not stored.
			LDA     (input) INTEGER, The leading dimension of the array A.  LDA >= max(1,M).
			IPIV    (output) INTEGER array, dimension (min(M,N)), The  pivot indices; 
		for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).
			INFO    (output) INTEGER, = 0:  successful exit ; < 0:  if INFO = -i, the i-th 
			argument had an illegal value ; > 0:  if INFO = i, U(i,i) is exactly  zero.  
			The factorization has  been  completed, but the factor U is exactly singular, 
			and division by zero will occur if it is used to solve a system of equations.
			*/	
			
		F77_CALL(dgetrf) (&ncolT, &ncolT, dpp, &ncolT, pivot, &exitcode);
		if(exitcode) 
			error(_("LAPACK routine dgetrf returned info code %d"), exitcode);
		
		/* call fortran routine dgetrf  with arguments
			TRANS   (input) CHARACTER*1, Specifies the form of the system of equations:
			= 'N':  A * X = B  (No transpose) ; = 'T':  A'* X = B  (Transpose) ;
		= 'C':  A'* X = B  (Conjugate transpose = Transpose)
			N       (input) INTEGER, The order of the matrix A.  N >= 0.
			NRHS    (input) INTEGER, The  number of right hand sides, i.e., the number of columns of
			the matrix B.  NRHS >= 0.
			A       (input) DOUBLE PRECISION array, dimension (LDA,N), The factors L and U 
			from the factorization A =  P*L*U  as  computed by DGETRF.
			LDA     (input) INTEGER, The leading dimension of the array A.  LDA >= max(1,N).
			IPIV    (input) INTEGER array, dimension (N), The pivot indices from DGETRF; for 1<=i<=N, 
			row i of the matrix was interchanged with row IPIV(i). 
			B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS), On entry, 
			the right hand side matrix B.  On exit, the  solution matrix X.
			LDB     (input) INTEGER, The leading dimension of the array B.  LDB >= max(1,N).
			INFO    (output) INTEGER, = 0:  successful exit; < 0: if INFO = -i, 
			the i-th argument had an illegal value
			*/	
			
		F77_CALL(dgetrs) ("N", &ncolT, &ncolT, dpp, &ncolT, pivot, npp, &ncolT, &exitcode);
		if(exitcode) 
			error(_("LAPACK routine dgetrs returned info code %d"), exitcode);
		
		// copy ncolTsqr character from npp to valT, i.e. valt <- npp 
		Memcpy(valT, npp, ncolTsqr);		
		
		/* undo preconditioning step 3: scaling */
		
		while(sqrpowscal --)
		{
			F77_CALL(dgemm) ("N", "N", &ncolT, &ncolT, &ncolT, &one, valT, &ncolT, valT, &ncolT, &zero, work, &ncolT);
			Memcpy(valT, work, ncolTsqr); 		// valT <- work
		}
		
		/* undo preconditioning step 2: inverse scaling */
		
		for(j = 0; j < ncolT; j++)
			for(i = 0; i < ncolT; i++)
				valT[i + j * ncolT] *= scale[i]/scale[j];
		
		// construct balancing permutation vector 
		for(i = 0; i < ncolT; i++)
			invperm[i] = i; //identity permutation
		
		// leading permutations applied in forward order 
		for(i = 0; i < (iloperm - 1); i++)
		{
			int permutedindex = (int) (perm[i]) - 1;
			int tmp = invperm[i];
			invperm[i] = invperm[permutedindex];
			invperm[permutedindex] = tmp;	
		}
		
		// trailing permutations applied in reverse order 
		for(i = ncolT - 1; i >= ihiperm; i--)
		{
			int permutedindex = (int) (perm[i]) - 1;
			int tmp = invperm[i];
			invperm[i] = invperm[permutedindex];
			invperm[permutedindex] = tmp;	
		}
		
		// construct inverse balancing permutation vector 
		Memcpy(pivot, invperm, ncolT); //pivot <- invperm
		for(i = 0; i < ncolT; i++)
			invperm[pivot[i]] = i;
		
		// apply inverse permutation 
		Memcpy(work, valT, ncolTsqr); //work <- valT
		for(j = 0; j < ncolT; j++)
			for(i = 0; i < ncolT; i++)
				valT[i + j * ncolT] = work[invperm[i] + invperm[j] * ncolT];
		
		/* undo preconditioning step 1 : trace denormalization */
		
		if(traceshiftT > 0)
		{
			double exptraceshift = exp(traceshiftT);
			for(i = 0; i < ncolTsqr; i++) 
				valT[i] *= exptraceshift;
		}
		
		/* compute u %*% e^( x[k] * T ) %*% v, for k in 1:lenX */
		//Rprintf("%d-%f\t",k,result[k] );
		for(i = 0; i < nrowT; i++) 
		{
			double temp = 0;
			for(j = 0; j < ncolT; j++) 
			{
				temp += valT[i + j* ncolT] * valV[j];
			}
			result[k] += temp * valU[i];
		//	Rprintf("%f\t",result[k] );	
		}
		//Rprintf("\t%f\n",result[k]);
		
		//plot
		/*
		Rprintf("matrix exp(x[k]*T) %d-\n",k);																																												
		for(i = 0; i < nrowT; i++) 
		{
			for(j = 0; j < ncolT; j++) 
				Rprintf("%f\t", valT[i + j* ncolT] );
			
			Rprintf("\n");		
		}
		Rprintf("\n");		
		*/	
	}	
				
	// clean up			
	UNPROTECT(6+k);			
		
	return resultinR;
}