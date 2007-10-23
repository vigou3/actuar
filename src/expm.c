/*  ===== actuar: an R package for Actuarial Science =====
 *
 *  Functions to compute the exponential of a matrix.
 *
 *  AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Christophe
 *  Dutang
 */

#include "expm.h"

/* Matrix exponential exp(x), where x is an (n x n) matrix. Result z
 * is an (n x n) matrix. */
static void expm(double *x, int n, double *z)
{
    /* declaration */	
	int nsqr = n*n;
	int nplus1 = n+1;	
	int i, j, k; //loop index
	int isUpperTriangular =0; //pseudo boolean: is x upper triangular ?
	int traceshiftx = 0; 
	int iloperm, ihiperm; //arguments for dgebal permutation
	int iloscal, ihiscal; //arguments for dgebal scaling	
	double *perm = (double *) R_alloc(n, sizeof(double)); //permutation array
	double *scale = (double *) R_alloc(n, sizeof(double)); //scale array
	int exitcode; //exitcode used for different fortran routine
	double *work = (double *) R_alloc(nsqr, sizeof(double)); //workspace array
	double infnorm; //infinite norm of matrix x
	int sqrpowscal ; //the square power used for scaling x
	double *npp = (double *) R_alloc(nsqr, sizeof(double)); //numerator power Pade'
	double *dpp = (double *) R_alloc(nsqr, sizeof(double)); //denominator power Pade'
	double minus1powj = -1; // (-1)^j 
	double one = 1.0; //useful for fortran routine dgemm
	double zero = 0.0;
	int *pivot = (int *) R_alloc(n, sizeof(int)); // pivot vector
	int *invperm = (int *) R_alloc(n, sizeof(int)); // inverse permutation vector
	
	/* special case where x is 1x1 matrix */
	if(n == 1)
	{
		z = exp(x); 	
		return;
	}
	
	
	
	//start
	z = x;
	
	//plot
	
	Rprintf("matrix z\n");
	for(i = 0; i < n; i++) 
	{
		for(j = 0; j < n; j++) 
			Rprintf("%f\t", z[i + j* n] );
		
		Rprintf("\n");		
	}
	Rprintf("\n");
	
	
	/* check if matrix x is upper triangular */
	// time cost O(n^2)
	isUpperTriangular = 1;
	//i row index and j column index
	for(j = 0; j < n-1; j++)
		for(i = j+1; i < n; i++) 
			isUpperTriangular  *= x[i + j* n] == 0; //x is stored column by column
	
			
	Rprintf("%d\n",isUpperTriangular);	
	
	/* step 1 of preconditioning : trace normalisation */
	// compute the average trace of matrix x
		
	for(i = 0; i < n; i++) 
		traceshiftT += x[i * nplus1];
			
	traceshiftT /= n;
		
	if(traceshiftT > 0)
	{
		for(i = 0; i < n; i++) 
			z[i * nplus1] -= traceshiftT;
	}
		
	
	/* step 2 of preconditioning : balancing with dgebal*/
	/* call fortran routine dgebal with arguments
		JOB			(input) CHARACTER*1
		Specifies the operations to be performed on A:
		=  'N':  none:  simply set ILO = 1, IHI = N, SCALE(I) = 1.0 for
		i = 1,...,N; = 'P':  permute only; = 'S':  scale only; 
		= 'B':  both permute and scale.			
		N			(input) INTEGER
		The order of the matrix A.  N >= 0.
		A			(input/output) DOUBLE PRECISION array, dimension (LDA,N)
		ILO			(output) INTEGER
		IHI			(output)  INTEGER ILO and IHI are set to integers such
		that on exit A such that A(i,j) = 0 if i > j and j <=ILO-1 or 
		j>=IHI+1,...,N. 
		SCALE		(output) DOUBLE PRECISION array, dimension (N)
		Details  of  the permutations and scaling factors applied to A.
		If P(j) is the index of the row and  column  interchanged  with
		row  and column j and D(j) is the scaling factor applied to row
		and column j, then SCALE(j) = P(j)    for  j  =  1,...,ILO-1  =
		D(j)		for j = ILO,...,IHI = P(j)    for j = IHI+1,...,N.  The
		order in which the interchanges are made is N to IHI+1, then  1
		to ILO-1.
		INFO		(output) INTEGER, =0:  successful exit; < 0:  if INFO = -i, 
		the i-th argument had an illegal value.
	*/
		
	//no need to permute T if T is upper triangular
	if(!isUpperTriangular)
	{			
		printf("blii\n");
		F77_CALL(dgebal) ("P", &n, z, &n, &iloperm, &ihiperm, perm, &exitcode);
		if(exitcode) 
			error(_("LAPACK routine dgebal returned info code %d when permuting T"), exitcode);
	}
	else
	{	//identity permutation
		iloperm = 1; 
		ihiperm = n;
	}	
	
	//scaling T in order to make the 1-norms of each row of T and its corresponding column nearly equal
	F77_CALL(dgebal) ("S", &n, z, &n, &iloscal, &ihiscal, scale, &exitcode);
	if(exitcode) 
		error(_("LAPACK routine dgebal returned info code %d when scaling T"), exitcode);			
		
				
	/* step 3 of preconditioning : scaling (a priori always needed) */			
	/* call fortran routine dlange to compute the infinite norm with arguments 
		NORM		(input) CHARACTER*1
		Specifies  the  value  to  be  returned  in DLANGE as described above.
		M			(input) INTEGER
		The number of rows of the matrix A.  M  >=  0.   When  M  =  0,
		DLANGE is set to zero.
		N			(input) INTEGER
		The  number  of  columns of the matrix A.  N >= 0.  When N = 0,
		DLANGE is set to zero.
		A			(input) DOUBLE PRECISION array, dimension (LDA,N)
		The m by n matrix A.
		LDA			(input) INTEGER
		The leading dimension of the array A.  LDA >= max(M,1).
		WORK		(workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
		where LWORK >= M when NORM = 'I'; otherwise, WORK is not referenced.
	*/		
				
	infnorm = F77_CALL(dlange) ("I", &n, &n, z, &n, work);
	if(infnorm > 0) 
		sqrpowscal = imax2( (int) 1 + log(infnorm)/log(2.0) , 0 );
	else
		sqrpowscal = 0;				
		
	if(sqrpowscal > 0) //scaling
	{
		double scalefactor = R_pow_di(2 , sqrpowscal);
		for(i = 0; i < nsqr; i++)
			z[i] /= scalefactor;
	}
		
			
		
	/* Pade' approximation (p=q=8) : compute T^8, T^7, T^6, ..., T^1*/			
	// init npp and dpp
	
	for(i = 0; i < nsqr; i++)
	{
		npp[i] = 0.0;
		dpp[i] = 0.0;
	}
	for(j = 7; j>=0; j--)
	{
	/* call fortran routine dgemm with arguments to compute: alpha*op(A)*op(B) + beta*C
		 TRANSA			CHARACTER*1.  On entry, TRANSA specifies the form of op( A )
		 to be used in the matrix multiplication as follows:
		 if = 'N' or 'n',  op( A ) = A. if = 'T' or 't',  op( A ) = A'.
		 if = 'C' or 'c',  op( A ) = A'.
		 TRANSB			CHARACTER*1.  On entry, TRANSB specifies the form of op( B )
		 to be used in the matrix multiplication as TRANSA
		 M				INTEGER. On entry, M specifies  the number  of rows  of the  
		 matrix op( A ) and of the matrix C. M >=0.
		 N				INTEGER. On entry, N specifies the number  of columns of the 
		 matrix op( B ) and the number of columns of the matrix C. N >=0.
		 K				INTEGER. On entry, K specifies  the number of columns of the 
		 matrix op( A ) and the number of rows of the matrix op( B ). K >=0.
		 ALPHA			DOUBLE PRECISION.
		 On entry, ALPHA specifies the scalar alpha.  
		 A				DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
		 k when TRANSA = 'N' or 'n', and is m otherwise. Before entry with TRANSA 
		 = 'N' or 'n', the leading m by k part of the array A must contain the matrix
		 A, otherwise the leading k by m  part of the array A must contain the matrix A.
		 LDA			INTEGER. On  entry, LDA specifies the first dimension of A as 
		 declared in the calling program. When TRANSA = 'N' or 'n' then  
		 LDA must be at least max( 1, m ), otherwise  LDA must be at least max( 1, k ).  
		 B				DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
		 n  when  TRANSB = 'N' or 'n', and  is k otherwise. Before entry  with TRANSB  
		 = 'N' or 'n', the leading k by n part of the array B must contain the matrix  
		 B, otherwise the leading n by k part of the array B must contain the matrix B.
		 Unchanged on exit.
		 LDB			INTEGER. On entry, LDB specifies the first dimension of B as declared  
		 in the calling program. When  TRANSB = 'N' or 'n' then LDB must be at least  
		 max( 1, k ), otherwise LDB must be at least max( 1, n ).  
		 BETA			DOUBLE PRECISION. On  entry, BETA specifies the scalar beta. When BETA  
		 is supplied as zero then C need not be set on input.
		 C				DOUBLE PRECISION array of DIMENSION ( LDC, n ). Before  entry, the 
		 leading m by n  part of the array  C must contain the matrix  C, except when 
		 beta is zero, in which case C need not be set on entry. On exit, the array C is
		 overwritten by the m by n  matrix ( alpha*op( A  )*op(  B  )  + beta*C ).
		 LDC			INTEGER. On  entry, LDC specifies the first dimension of C as declared in
		 the  calling program. LDC must be at least max( 1,m ). 
	 */
			 
		// compute npp = z * npp + padec88[j] * z 
			// first matrix operation: work <- z %*% npp
		F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, npp, &n, &zero, work, &n);
				
			// second : npp <- work + padec88[j] * z
		for(i = 0; i < nsqr; i++)
			npp[i] = work[i] + padec88[j] * z[i];
			
		// compute dpp = z * dpp + (-1)^j * padec88[j] * z 
			// first matrix operation: work <- T %*% dpp
		F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, dpp, &n, &zero, work, &n);
			
			// second : dpp <- work + (-1)^j * padec88[j] * z
		for(i = 0; i < nsqr; i++)
			dpp[i] = work[i] + minus1powj * padec88[j] * valT[i];
			
		// (-1)^j
		minus1powj *= -1;	
	}
		
	// compute T^0 
	for(i = 0; i < nsqr; i++) 
		dpp[i] *= -1.0;
			
	for(j = 0; j < n; j++)
	{
		npp[j * nplus1] += 1.0;
		dpp[j * nplus1] += 1.0;		
	}									
		
	/* Compute pade approximation = (dpp)^-1 * npp. */
			
	// first compute LU factorisation  A = P * L * U, permutation P is written in IPIV argument																	
	/* call fortran routine dgetrf  with arguments
		M			(input) INTEGER, The number of rows of the matrix A.  M >= 0.
		N			(input) INTEGER, The number of columns of the matrix A.  N >= 0.
		A			(input/output) DOUBLE PRECISION array, dimension (LDA,N),  On entry, 
		the M-by-N matrix to be factored.  On exit, the factors L and U from 
		the factorization A = P*L*U; the unit diagonal elements of L are not stored.
		LDA			(input) INTEGER, The leading dimension of the array A.  LDA >= max(1,M).
		IPIV		(output) INTEGER array, dimension (min(M,N)), The  pivot indices; 
		for 1 <= i <= min(M,N), row i of the matrix was interchanged with row IPIV(i).
		INFO		(output) INTEGER, = 0:  successful exit ; < 0:  if INFO = -i, the i-th 
		argument had an illegal value ; > 0:  if INFO = i, U(i,i) is exactly  zero.  
		The factorization has  been  completed, but the factor U is exactly singular, 
		and division by zero will occur if it is used to solve a system of equations.
	*/	
			
	F77_CALL(dgetrf) (&n, &n, dpp, &n, pivot, &exitcode);
	if(exitcode) 
		error(_("LAPACK routine dgetrf returned info code %d"), exitcode);
		
	// second compute A * X = B using the LU factorization	
	/* call fortran routine dgetrs  with arguments
		TRANS		(input) CHARACTER*1, Specifies the form of the system of equations:
		= 'N':  A * X = B  (No transpose) ; = 'T':  A'* X = B  (Transpose) ;
		= 'C':  A'* X = B  (Conjugate transpose = Transpose)
		N			(input) INTEGER, The order of the matrix A.  N >= 0.
		NRHS		(input) INTEGER, The  number of right hand sides, i.e., the number of columns of
		the matrix B.  NRHS >= 0.
		A			(input) DOUBLE PRECISION array, dimension (LDA,N), The factors L and U 
		from the factorization A =  P*L*U  as  computed by DGETRF.
		LDA			(input) INTEGER, The leading dimension of the array A.  LDA >= max(1,N).
		IPIV		(input) INTEGER array, dimension (N), The pivot indices from DGETRF; for 1<=i<=N, 
		row i of the matrix was interchanged with row IPIV(i). 
		B			(input/output) DOUBLE PRECISION array, dimension (LDB,NRHS), On entry, 
		the right hand side matrix B.  On exit, the  solution matrix X.
		LDB			(input) INTEGER, The leading dimension of the array B.  LDB >= max(1,N).
		INFO		(output) INTEGER, = 0:  successful exit; < 0: if INFO = -i, 
		the i-th argument had an illegal value
	*/	
			
	F77_CALL(dgetrs) ("N", &n, &n, dpp, &n, pivot, npp, &n, &exitcode);
	if(exitcode) 
		error(_("LAPACK routine dgetrs returned info code %d"), exitcode);
		
	// copy nsqr character from npp to z, i.e. z <- npp 
	Memcpy(z, npp, nsqr);		
	
	
	/* undo preconditioning step 3: scaling */		
	while(sqrpowscal --)
	{
		// compute matrix operation: work <- z %*% z
		F77_CALL(dgemm) ("N", "N", &n, &n, &n, &one, z, &n, z, &n, &zero, work, &n);
		Memcpy(z, work, nsqr); 		// z <- work
	}
		
	/* undo preconditioning step 2: inverse balancing */
	//debalancing
	for(j = 0; j < n; j++)
		for(i = 0; i < n; i++)
			valT[i + j * n] *= scale[i]/scale[j];
		
	//inverse permuation if T is not upper triangular and 'perm' is not the identity permutation
	if( (iloperm != 1 || ihiperm != n) && !isUpperTriangular)
	{
		// construct balancing permutation vector 
		for(i = 0; i < n; i++)
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
		for(i = n - 1; i >= ihiperm; i--)
		{
			int permutedindex = (int) (perm[i]) - 1;
			int tmp = invperm[i];
			invperm[i] = invperm[permutedindex];
			invperm[permutedindex] = tmp;	
		}
			
		// construct inverse balancing permutation vector 
		Memcpy(pivot, invperm, n); //pivot <- invperm
		for(i = 0; i < n; i++)
			invperm[pivot[i]] = i;
			
		// apply inverse permutation 
		Memcpy(work, z, nsqr); //work <- z
		for(j = 0; j < n; j++)
			for(i = 0; i < n; i++)
				z[i + j * n] = work[invperm[i] + invperm[j] * n];
	}	
		
	/* undo preconditioning step 1 : trace denormalization */
	if(traceshiftT > 0)
	{
		double exptraceshift = exp(traceshiftT);
		for(i = 0; i < nsqr; i++) 
			z[i] *= exptraceshift;
	}		
	
		//plot
		
		Rprintf("matrix exp(x) -\n");																																												
		for(i = 0; i < n; i++) 
		{
			for(j = 0; j < n; j++) 
				Rprintf("%f\t", z[i + j* ncolT] );
			
			Rprintf("\n");		
		}
		Rprintf("\n");		
													
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
