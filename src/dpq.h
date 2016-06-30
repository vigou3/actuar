/*  ===== actuar: An R Package for Actuarial Science =====
 *
 * Utilities for `dpq' handling (density/probability/quantile)
 *
 * These (except the last one) are copied from nmath/dpq.h in the R
 * sources with the names changed from "R_" to "ACT_".
 *
 *  AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>
 *          with much indirect help from the R Core Team
 */


/* give_log in "d" & "mgf";  log_p in "p" & "q" : */
#define give_log log_p

#define ACT_D__0  (log_p ? R_NegInf : 0.)
#define ACT_D__1  (log_p ? 0. : 1.)
#define ACT_DT_0  (lower_tail ? ACT_D__0 : ACT_D__1)
#define ACT_DT_1  (lower_tail ? ACT_D__1 : ACT_D__0)

/* Use 0.5 - p + 0.5 to perhaps gain 1 bit of accuracy */
#define ACT_D_Lval(p)     (lower_tail ? (p) : (0.5 - (p) + 0.5))  /*  p  */
#define ACT_D_Cval(p)     (lower_tail ? (0.5 - (p) + 0.5) : (p))  /*  1 - p */

#define ACT_D_val(x)      (log_p  ? log(x) : (x))         /*  x  in pF(x,..) */
#define ACT_D_qIv(p)      (log_p  ? exp(p) : (p))         /*  p  in qF(p,..) */
#define ACT_D_exp(x)      (log_p  ?  (x)   : exp(x))      /* exp(x) */
#define ACT_D_Clog(p)     (log_p  ? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */

#define ACT_DT_val(x)     (lower_tail ? ACT_D_val(x)  : ACT_D_Clog(x))
#define ACT_DT_Cval(x)    (lower_tail ? ACT_D_Clog(x) : ACT_D_val(x))

// log(1 - exp(x))  in more stable form than log1p(- R_D_qIv(x)) :
#define ACT_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))

/*Boundaries*/
#define ACT_Q_P01_boundaries(p, _LEFT_, _RIGHT_)          \
    if (log_p) {                                        \
        if(p > 0)                                       \
            return R_NaN;				\
        if(p == 0) /* upper bound*/                     \
            return lower_tail ? _RIGHT_ : _LEFT_;       \
        if(p == R_NegInf)                               \
            return lower_tail ? _LEFT_ : _RIGHT_;       \
    }                                                   \
    else { /* !log_p */                                 \
        if(p < 0 || p > 1)                              \
            return R_NaN;				\
        if(p == 0)                                      \
            return lower_tail ? _LEFT_ : _RIGHT_;       \
        if(p == 1)                                      \
            return lower_tail ? _RIGHT_ : _LEFT_;       \
    }

/* Infinite limit in "lev" */
#define ACT_DLIM__0(x, y)   (R_FINITE(x) ? R_pow(x, y) : 0.)
