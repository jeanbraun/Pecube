/* dtime_compat.c — provides dtime_() for compilers lacking the GNU extension.
 * C symbol dtime_ matches exactly what flang generates for Fortran calls to dtime().
 * tarray[0] = user CPU time elapsed since last call (seconds)
 * tarray[1] = system time (not separately available via clock(), set to 0)
 */
#include <time.h>

static clock_t last_clock = 0;
static int initialized = 0;

float dtime_(float *tarray) {
    clock_t now = clock();
    float elapsed;
    if (!initialized) {
        last_clock = now;
        initialized = 1;
    }
    elapsed = (float)(now - last_clock) / CLOCKS_PER_SEC;
    last_clock = now;
    tarray[0] = elapsed;
    tarray[1] = 0.0f;
    return elapsed;
}