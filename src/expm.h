#include <R.h>
#include <Rinternals.h>
#include "locale.h"
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <Rmath.h>


// manual of lapack routine 
// http://www.mathkeisan.com/UsersGuide/man/

/* for matrix exponential calculation. */
/* Pade' approximation of e^A is defined by
		Rpq(A) = Dpq(A))^-1 * Npq(A) ,
   where
		Npq(A) = sum_{i=0}^p n_{pqj} A^j and Dpq(A) = sum_{i=0}^p d_{pqj} (-A)^j
   with n_{pqj} = [ (p+q-j)!p! ]/[ (p+q)!j!(p-j)! ] and 
		d_{pqj} = [ (p+q-j)!q! ]/[ (p+q)!j!(q-j)! ]
   The following constants are these n_{88j} and d_{88j}, i.e.
		padec[j-1] = n_{88j} = (-1)^j * d_{88j}, 1 <= j <= p.
 */
const static double padec88 [] = 
{
  5.0000000000000000e-1,
  1.1666666666666667e-1,
  1.6666666666666667e-2,
  1.6025641025641026e-3,
  1.0683760683760684e-4,
  4.8562548562548563e-6,
  1.3875013875013875e-7,
  1.9270852604185938e-9,
};

