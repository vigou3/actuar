/* Localization */
#include <R.h>
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("actuar", String)
#else
#define _(String) (String)
#endif
