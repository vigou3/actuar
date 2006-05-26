/* Functions accessed from .External() */
SEXP do_dpq(SEXP args);
SEXP do_dpq1(int code, SEXP args);
SEXP do_dpq2(int code, SEXP args);
SEXP do_dpq3(int code, SEXP args);
SEXP do_dpq4(int code, SEXP args);

SEXP do_random(SEXP args);
SEXP do_random1(int code, SEXP args);
SEXP do_random2(int code, SEXP args);
SEXP do_random3(int code, SEXP args);
SEXP do_random4(int code, SEXP args);

/* Utility functions */
/*   One parameter distributions */
double dinvexp(double x, double scale, int give_log);
double rinvexp(double scale);

/*   Two parameter distributions */
double rllogis(double shape, double scale);
double rparalogis(double shape, double scale);

double dpareto(double x, double shape, double scale, int give_log);
double ppareto(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto(double shape, double scale);

double dinvpareto(double x, double shape, double scale, int give_log);
double rinvpareto(double shape, double scale);

double dpareto1(double x, double shape, double scale, int give_log);
double ppareto1(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto1(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto1(double shape, double scale);

double dinvparalogis(double x, double shape, double scale, int give_log);
double rinvparalogis(double shape, double scale);

double dinvweibull(double x, double scale, double shape, int give_log);
double rinvweibull(double scale, double shape);

/*   Three parameter distributions */
double dburr(double x, double shape1, double scale, double shape2, int give_log);
double pburr(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qburr(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rburr(double shape1, double scale, double shape2);

double dgenpareto(double x, double shape1, double scale, double shape2, int give_log);
double pgenpareto(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qgenpareto(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rgenpareto(double shape1, double scale, double shape2);

double dinvburr(double x, double shape1, double scale, double shape2, int give_log)
double rinvburr(double shape, double scale, double gamma);

/*   Four parameter distributions */
double dtrbeta(double x, double shape1, double scale, double shape2, double shape3, int give_log);
double ptrbeta(double q, double shape1, double scale, double shape2, double shape3, int lower_tail, int log_p);
double qtrbeta(double p, double shape1, double scale, double shape2, double shape3, int lower_tail, int log_p);
double rtrbeta(double shape1, double scale, double shape2, double shape3);


/* Definitions for the table linking the first group of functions to
 * the second one. Table found in names.c */
typedef struct {
    char *name; 
    SEXP (*cfun)(int, SEXP);
    int code;
} FUNTAB;
extern FUNTAB fun_tab[];
