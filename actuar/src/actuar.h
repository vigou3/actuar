/* Functions accessed from .External() */
SEXP do_random(SEXP args);
SEXP do_random1(int code, SEXP args);
SEXP do_random2(int code, SEXP args);
SEXP do_random3(int code, SEXP args);
SEXP do_random4(int code, SEXP args);

/* Utility functions */
double rpareto(double shape, double scale);

/* Definitions for the table linking the first group of functions to
 * the second one. Table found in names.c */
typedef struct {
    char *name; 
    SEXP (*cfun)(int, SEXP);
    int code;
} FUNTAB;
extern FUNTAB fun_tab[];




