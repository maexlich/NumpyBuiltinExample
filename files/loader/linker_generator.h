#include <stdbool.h>
#ifdef __linux
#define MAX_PATH 4096
#endif
int boinc_temporary_exit_wrapper(int, const char *, int);
double boinc_get_ram(double);


#define TRUE 1
#define FALSE 0