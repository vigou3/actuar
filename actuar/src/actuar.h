/* Functions accessed from .External() */
SEXP do_dpq(SEXP args);
SEXP do_dpq1(int code, SEXP args);
SEXP do_dpq2(int code, SEXP args);
SEXP do_dpq3(int code, SEXP args);
SEXP do_dpq4(int code, SEXP args);
SEXP do_dpq5(int code, SEXP args);

SEXP do_random(SEXP args);
SEXP do_random1(int code, SEXP args);
SEXP do_random2(int code, SEXP args);
SEXP do_random3(int code, SEXP args);
SEXP do_random4(int code, SEXP args);

/* Utility functions */
/*   One parameter distributions */
double dinvexp(double x, double scale, int give_log);
double pinvexp(double q, double scale, int lower_tail, int log_p);
double qinvexp(double p, double scale, int lower_tail, int log_p);
double rinvexp(double scale);

/*   Two parameter distributions */
double dinvgamma(double x, double scale, double shape, int give_log);
double pinvgamma(double q, double scale, double shape, int lower_tail, int log_p);
double qinvgamma(double p, double scale, double shape, int lower_tail, int log_p);
double rinvgamma(double scale, double shape);

double dlgamma(double x, double shapelog, double ratelog, int give_log);
double plgamma(double q, double shapelog, double ratelog, int lower_tail, int log_p);
double qlgamma(double p, double shapelog, double ratelog, int lower_tail, int log_p);
double rlgamma(double ratelog, double shapelog);

double dinvparalogis(double x, double shape, double scale, int give_log);
double pinvparalogis(double q, double shape, double scale, int lower_tail, int log_p);
double qinvparalogis(double p, double shape, double scale, int lower_tail, int log_p);
double rinvparalogis(double shape, double scale);

double dinvpareto(double x, double shape, double scale, int give_log);
double pinvpareto(double q, double shape, double scale, int lower_tail, int log_p);
double qinvpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rinvpareto(double shape, double scale);

double dinvweibull(double x, double scale, double shape, int give_log);
double pinvweibull(double q, double scale, double shape, int lower_tail, int log_p);
double qinvweibull(double p, double scale, double shape, int lower_tail, int log_p);
double rinvweibull(double scale, double shape);

double dllogis(double x, double shape, double scale, int give_log);
double pllogis(double q, double shape, double scale, int lower_tail, int log_p);
double qllogis(double p, double shape, double scale, int lower_tail, int log_p);
double rllogis(double shape, double scale);

double dparalogis(double x, double shape, double scale, int give_log);
double pparalogis(double q, double shape, double scale, int lower_tail, int log_p);
double qparalogis(double p, double shape, double scale, int lower_tail, int log_p);
double rparalogis(double shape, double scale);

double dpareto(double x, double shape, double scale, int give_log);
double ppareto(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto(double shape, double scale);

double dpareto1(double x, double shape, double scale, int give_log);
double ppareto1(double q, double shape, double scale, int lower_tail, int log_p);
double qpareto1(double p, double shape, double scale, int lower_tail, int log_p);
double rpareto1(double shape, double scale);

/*   Three parameter distributions */
double dburr(double x, double shape1, double scale, double shape2, int give_log);
double pburr(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qburr(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rburr(double shape1, double scale, double shape2);
double mburr(double k, double shape1, double scale, double shape2, int give_log);
double levburr(double x, double shape1, double scale, double shape2, double order, int give_log);

double dgenpareto(double x, double shape1, double scale, double shape2, int give_log);
double pgenpareto(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qgenpareto(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rgenpareto(double shape1, double scale, double shape2);
double mgenpareto(double k, double shape1, double scale, double shape2, int give_log);
double levgenpareto(double x, double shape1, double scale, double shape2, double order, int give_log);

double dinvburr(double x, double shape1, double scale, double shape2, int give_log);
double pinvburr(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qinvburr(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rinvburr(double shape1, double scale, double shape2);
double minvburr(double k, double shape1, double scale, double shape2, int give_log);
double levinvburr(double x, double shape1, double scale, double shape2, double order, int give_log);

double dtrgamma(double x, double shape1, double scale, double shape2, int give_log);
double ptrgamma(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qtrgamma(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rtrgamma(double shape1, double scale, double shape2);
double mtrgamma(double k, double shape1, double scale, double shape2, int give_log);
double levtrgamma(double x, double shape1, double scale, double shape2, double order, int give_log);

double dinvtrgamma(double x, double shape1, double scale, double shape2, int give_log);
double pinvtrgamma(double q, double shape1, double scale, double shape2, int lower_tail, int log_p);
double qinvtrgamma(double p, double shape1, double scale, double shape2, int lower_tail, int log_p);
double rinvtrgamma(double shape1, double scale, double shape2);
double minvtrgamma(double k, double shape1, double scale, double shape2, int give_log);
double levinvtrgamma(double x, double shape1, double scale, double shape2, double order, int give_log);

/*   Four parameter distributions */
double dtrbeta(double x, double shape1, double scale, double shape2, double shape3, int give_log);
double ptrbeta(double q, double shape1, double scale, double shape2, double shape3, int lower_tail, int log_p);
double qtrbeta(double p, double shape1, double scale, double shape2, double shape3, int lower_tail, int log_p);
double rtrbeta(double shape1, double scale, double shape2, double shape3);
double mtrbeta(double k, double shape1, double scale, double shape2, double shape3, int give_log);
double levtrbeta(double x, double shape1, double scale, double shape2, double shape3, double order, int give_log);


/* Definitions for the table linking the first group of functions to
 * the second one. Table found in names.c */
typedef struct {
    char *name; 
    SEXP (*cfun)(int, SEXP);
    int code;
} FUNTAB;
extern FUNTAB fun_tab[];
