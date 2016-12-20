#include <Rinternals.h>

/*Error messages */
#define R_MSG_NA        _("NaNs produced")

/* Interfaces to routines from package expint */
double(*actuar_gamma_inc)(double,double);

/* Functions accessed from .External() */
SEXP actuar_do_dpq(SEXP args);
SEXP actuar_do_dpq0(int code, SEXP args);
SEXP actuar_do_dpq1(int code, SEXP args);
SEXP actuar_do_dpq2(int code, SEXP args);
SEXP actuar_do_dpq3(int code, SEXP args);
SEXP actuar_do_dpq4(int code, SEXP args);
SEXP actuar_do_dpq5(int code, SEXP args);

SEXP actuar_do_random(SEXP args);
SEXP actuar_do_random1(int code, SEXP args, SEXPTYPE type);
SEXP actuar_do_random2(int code, SEXP args, SEXPTYPE type);
SEXP actuar_do_random3(int code, SEXP args, SEXPTYPE type);
SEXP actuar_do_random4(int code, SEXP args, SEXPTYPE type);

SEXP actuar_do_dpqphtype(SEXP args);
SEXP actuar_do_dpqphtype2(int code, SEXP args);

SEXP actuar_do_randomphtype(SEXP args);
SEXP actuar_do_randomphtype2(int code, SEXP args, SEXPTYPE type);

SEXP actuar_do_hierarc(SEXP args);
SEXP actuar_do_panjer(SEXP args);

/* Utility functions */
/*   Matrix algebra */
void actuar_expm(double *x, int n, double *z);
double actuar_expmprod(double *x, double *M, double *y, int n);
void actuar_matpow(double *x, int n, int k, double *z);
void actuar_solve(double *A, double *B, int n, int p, double *z);

/*   Special integrals */
double betaint(double x, double a, double b, int foo);
double betaint_raw(double x, double a, double b);

/*   Sampling */
int SampleSingleValue(int n, double *p);

/*   One parameter distributions */
double mexp(double order, double scale, int give_log);
double levexp(double limit, double scale, double order, int give_log);
double mgfexp(double t, double scale, int give_log);

double dinvexp(double x, double scale, int give_log);
double pinvexp(double q, double scale, int lower_tail, int log_p);
double qinvexp(double p, double scale, int lower_tail, int log_p);
double rinvexp(double scale);
double minvexp(double order, double scale, int give_log);
double levinvexp(double limit, double scale, double order, int give_log);

double dlogarithmic(double x, double p, int give_log);
double plogarithmic(double x, double p, int lower_tail, int log_p);
double qlogarithmic(double x, double p, int lower_tail, int log_p);
double rlogarithmic(double p);

double dztpois(double x, double lambda, int give_log);
double pztpois(double q, double lambda, int lower_tail, int log_p);
double qztpois(double p, double lambda, int lower_tail, int log_p);
double rztpois(double lambda);

double dztgeom(double x, double prob, int give_log);
double pztgeom(double q, double prob, int lower_tail, int log_p);
double qztgeom(double p, double prob, int lower_tail, int log_p);
double rztgeom(double prob);

/*   Two parameter distributions */
double munif(double order, double min, double max, int give_log);
double levunif(double limit, double min, double max, double order, int give_log);
double mgfunif(double t, double min, double max, int give_log);

double mnorm(double order, double mean, double sd, int give_log);
double mgfnorm(double t, double mean, double sd, int give_log);

double mbeta(double order, double shape1, double shape2, int give_log);
double levbeta(double limit, double shape1, double shape2, double order, int give_log);

double mgamma(double order, double shape, double scale, int give_log);
double levgamma(double limit, double shape, double scale, double order, int give_log);
double mgfgamma(double t, double shape, double scale, int give_log);

double mchisq(double order, double df, double ncp, int give_log);
double levchisq(double limit, double df, double ncp, double order, int give_log);
double mgfchisq(double t, double df, double ncp, int give_log);

double dinvgamma(double x, double scale, double shape, int give_log);
double pinvgamma(double q, double scale, double shape, int lower_tail, int log_p);
double qinvgamma(double p, double scale, double shape, int lower_tail, int log_p);
double rinvgamma(double scale, double shape);
double minvgamma(double order, double scale, double shape, int give_log);
double levinvgamma(double limit, double scale, double shape, double order, int give_log);
double mgfinvgamma(double t, double shape, double scale, int give_log);

double dinvparalogis(double x, double shape, double scale, int give_log);
double pinvparalogis(double q, double shape, double scale, int lower_tail, int log_p);
double qinvparalogis(double p, double shape, double scale, int lower_tail, int log_p);
double rinvparalogis(double shape, double scale);
double minvparalogis(double order, double shape, double scale, int give_log);
double levinvparalogis(double limit, double shape, double scale, double order, int give_log);

double dinvpareto(double x, double shape, double scale, int give_log);
double pinvpareto(double q, double shape, double scale, int lower_tail, int log_p);
double qinvpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rinvpareto(double shape, double scale);
double minvpareto(double order, double shape, double scale, int give_log);
double levinvpareto(double limit, double shape, double scale, double order, int log_p);

double dinvweibull(double x, double scale, double shape, int give_log);
double pinvweibull(double q, double scale, double shape, int lower_tail, int log_p);
double qinvweibull(double p, double scale, double shape, int lower_tail, int log_p);
double rinvweibull(double scale, double shape);
double minvweibull(double order, double scale, double shape, int give_log);
double levinvweibull(double limit, double scale, double shape, double order, int give_log);

double dlgamma(double x, double shapelog, double ratelog, int give_log);
double plgamma(double q, double shapelog, double ratelog, int lower_tail, int log_p);
double qlgamma(double p, double shapelog, double ratelog, int lower_tail, int log_p);
double rlgamma(double ratelog, double shapelog);
double mlgamma(double order, double shapelog, double ratelog, int give_log);
double levlgamma(double limit, double shapelog, double ratelog, double order, int give_log);

double dllogis(double x, double shape, double scale, int give_log);
double pllogis(double q, double shape, double scale, int lower_tail, int log_p);
double qllogis(double p, double shape, double scale, int lower_tail, int log_p);
double rllogis(double shape, double scale);
double mllogis(double order, double shape, double scale, int give_log);
double levllogis(double limit, double shape, double scale, double order, int give_log);

double mlnorm(double order, double logmean, double logsd, int give_log);
double levlnorm(double limit, double logmean, double logsd, double order, int give_log);

double dparalogis(double x, double shape, double scale, int give_log);
double pparalogis(double q, double shape, double scale, int lower_tail, int log_p);
double qparalogis(double p, double shape, double scale, int lower_tail, int log_p);
double rparalogis(double shape, double scale);
double mparalogis(double order, double shape, double scale, int give_log);
double levparalogis(double limit, double shape, double scale, double order, int give_log);

double dpareto(double x, double shape, double scale, int give_log);
double ppareto(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto(double shape, double scale);
double mpareto(double order, double shape, double scale, int give_log);
double levpareto(double limit, double shape, double scale, double order, int give_log);

double dpareto1(double x, double shape, double scale, int give_log);
double ppareto1(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto1(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto1(double shape, double scale);
double mpareto1(double order, double shape, double scale, int give_log);
double levpareto1(double limit, double shape, double scale, double order, int give_log);

double mweibull(double order, double scale, double shape, int give_log);
double levweibull(double limit, double scale, double shape, double order, int give_log);

double minvGauss(double order, double nu, double lambda, int give_log); /* deprecated v2.0-0 */
double levinvGauss(double limit, double nu, double lambda, double order, int give_log); /* idem */
double mgfinvGauss(double t, double nu, double lambda, int give_log); /* idem */

double dgumbel(double x, double alpha, double beta, int give_log);
double pgumbel(double q, double alpha, double beta, int lower_tail, int log_p);
double qgumbel(double p, double alpha, double beta, int lower_tail, int log_p);
double rgumbel(double alpha, double beta);
double mgumbel(double order, double alpha, double beta, int give_log);
double mgfgumbel(double t, double alpha, double beta, int give_log);

double dinvgauss(double x, double mu, double phi, int give_log);
double pinvgauss(double q, double mu, double phi, int lower_tail, int log_p);
double qinvgauss(double q, double mu, double phi, int lower_tail, int log_p,
		 double tol, int maxit, int echo);
double rinvgauss(double mu, double phi);
double minvgauss(double order, double mean, double phi, int give_log);
double levinvgauss(double limit, double mean, double phi, double order, int give_log);
double mgfinvgauss(double t, double mean, double phi, int give_log);

double dztnbinom(double x, double size, double prob, int give_log);
double pztnbinom(double q, double size, double prob, int lower_tail, int log_p);
double qztnbinom(double p, double size, double prob, int lower_tail, int log_p);
double rztnbinom(double size, double prob);

double dztbinom(double x, double size, double prob, int give_log);
double pztbinom(double q, double size, double prob, int lower_tail, int log_p);
double qztbinom(double p, double size, double prob, int lower_tail, int log_p);
double rztbinom(double size, double prob);

double dzmlogarithmic(double x, double p, double p0m, int give_log);
double pzmlogarithmic(double x, double p, double p0m, int lower_tail, int log_p);
double qzmlogarithmic(double x, double p, double p0m, int lower_tail, int log_p);
double rzmlogarithmic(double p, double p0m);

double dzmpois(double x, double lambda, double p0m, int give_log);
double pzmpois(double q, double lambda, double p0m, int lower_tail, int log_p);
double qzmpois(double p, double lambda, double p0m, int lower_tail, int log_p);
double rzmpois(double lambda, double p0m);

double dzmgeom(double x, double prob, double p0m, int give_log);
double pzmgeom(double q, double prob, double p0m, int lower_tail, int log_p);
double qzmgeom(double p, double prob, double p0m, int lower_tail, int log_p);
double rzmgeom(double prob, double p0m);
double rzmgeom2(double prob, double p0m);

double dpoisinvgauss(double x, double mu, double phi, int give_log);
double ppoisinvgauss(double q, double mu, double phi, int lower_tail, int log_p);
double qpoisinvgauss(double p, double mu, double phi, int lower_tail, int log_p);
double rpoisinvgauss(double mu, double phi);

/*   Three parameter distributions */
double dburr(double x, double shape1, double shape2, double scale, int give_log);
double pburr(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qburr(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rburr(double shape1, double shape2, double scale);
double mburr(double order, double shape1, double shape2, double scale, int give_log);
double levburr(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dgenpareto(double x, double shape1, double shape2, double scale, int give_log);
double pgenpareto(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qgenpareto(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rgenpareto(double shape1, double shape2, double scale);
double mgenpareto(double order, double shape1, double shape2, double scale, int give_log);
double levgenpareto(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dinvburr(double x, double shape1, double shape2, double scale, int give_log);
double pinvburr(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qinvburr(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rinvburr(double shape1, double shape2, double scale);
double minvburr(double order, double shape1, double shape2, double scale, int give_log);
double levinvburr(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dinvtrgamma(double x, double shape1, double shape2, double scale, int give_log);
double pinvtrgamma(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qinvtrgamma(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rinvtrgamma(double shape1, double shape2, double scale);
double minvtrgamma(double order, double shape1, double shape2, double scale, int give_log);
double levinvtrgamma(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dtrgamma(double x, double shape1, double shape2, double scale, int give_log);
double ptrgamma(double q, double shape1, double shape2, double scale, int lower_tail, int log_p);
double qtrgamma(double p, double shape1, double shape2, double scale, int lower_tail, int log_p);
double rtrgamma(double shape1, double shape2, double scale);
double mtrgamma(double order, double shape1, double shape2, double scale, int give_log);
double levtrgamma(double limit, double shape1, double shape2, double scale, double order, int give_log);

double dzmnbinom(double x, double size, double prob, double p0m, int give_log);
double pzmnbinom(double q, double size, double prob, double p0m, int lower_tail, int log_p);
double qzmnbinom(double p, double size, double prob, double p0m, int lower_tail, int log_p);
double rzmnbinom(double size, double prob, double p0m);

double dzmbinom(double x, double size, double prob, double p0m, int give_log);
double pzmbinom(double q, double size, double prob, double p0m, int lower_tail, int log_p);
double qzmbinom(double p, double size, double prob, double p0m, int lower_tail, int log_p);
double rzmbinom(double size, double prob, double p0m);

/*   Four parameter distributions */
double dgenbeta(double x, double shape1, double shape2, double shape3, double scale, int give_log);
double pgenbeta(double q, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double qgenbeta(double p, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double rgenbeta(double shape1, double shape2, double shape3, double scale);
double mgenbeta(double order, double shape1, double shape2, double shape3, double scale, int give_log);
double levgenbeta(double limit, double shape1, double shape2, double shape3, double scale, double order, int give_log);

double dtrbeta(double x, double shape1, double shape2, double shape3, double scale, int give_log);
double ptrbeta(double q, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double qtrbeta(double p, double shape1, double shape2, double shape3, double scale, int lower_tail, int log_p);
double rtrbeta(double shape1, double shape2, double shape3, double scale);
double mtrbeta(double order, double shape1, double shape2, double shape3, double scale, int give_log);
double levtrbeta(double limit, double shape1, double shape2, double shape3, double scale, double order, int give_log);

/*   Phase-type distributions */
double dphtype(double x, double *pi, double *T, int m, int give_log);
double pphtype(double x, double *pi, double *T, int m, int lower_tail, int log_p);
double rphtype(double *pi, double **Q, double *rates, int m);
double mphtype(double order, double *pi, double *T, int m, int give_log);
double mgfphtype(double x, double *pi, double *T, int m, int give_log);


/* Definitions for the tables linking the first group of functions to
 * the second one. Tables found in names.c. One table for
 * {d,p,q,m,lev} functions and one for the {r} functions since we
 * need one more argument: the type of the result. */
typedef struct {
    char *name;
    SEXP (*cfun)(int, SEXP);
    int code;
} DPQTAB;
extern DPQTAB dpq_tab[];

typedef struct {
    char *name;
    SEXP (*cfun)(int, SEXP, SEXPTYPE);
    int code;
    SEXPTYPE type;
} RANDOMTAB;
extern RANDOMTAB random_tab[];
