/* Utilities for `dpq' handling (density/probability/quantile) */

/* give_log in "d";  log_p in "p" & "q" : */
#define give_log log_p

#define R_D__0	(log_p ? R_NegInf : 0.)
#define R_D__1	(log_p ? 0. : 1.)
#define R_D_d0	(give_log ? R_NegInf : 0.)
#define R_D_d1	(give_log ? 0. : 1.)
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)
#define R_D_qIv(p)	(log_p	? exp(p) : (p))
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))

/*Boundaries*/
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    error(_("NAN"));				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == R_NegInf)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    error(_("NAN"));				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }


#ifndef M_PI
#define M_PI		3.141592653589793238462643383280	/* pi */
#endif

#ifndef M_2PI
#define M_2PI		6.283185307179586476925286766559	/* 2*pi */
#endif

#define R_NaN  (0.0/0.0)
