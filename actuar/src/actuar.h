/* Functions accessed from .External() */
SEXP do_random(SEXP args);
SEXP do_random1(int code, SEXP args);
SEXP do_random2(int code, SEXP args);
SEXP do_random3(int code, SEXP args);
SEXP do_random4(int code, SEXP args);

/* Utility functions */
double rbetatrans(double shape, double scale, double gamma, double tau);
double rparetogen(double shape, double scale, double tau);
double rburr(double shape, double scale, double gamma);
double rloglogistic(double shape, double scale);
double rparalogistic(double shape, double scale);
double rpareto(double shape, double scale);
double rspareto(double shape, double scale);

double rinverseburr(double tau, double scale, double gamma);
double rinverseparalogistic(double tau, double scale);
double rinversepareto(double tau, double scale);
double rinverseweibull(double scale, double tau);
double rinverseexp(double scale);

/* Definitions for the table linking the first group of functions to
 * the second one. Table found in names.c */
typedef struct {
    char *name; 
    SEXP (*cfun)(int, SEXP);
    int code;
} FUNTAB;
extern FUNTAB fun_tab[];

