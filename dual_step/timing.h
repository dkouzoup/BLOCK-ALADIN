#ifndef TIMING_FUNCTIONS_H
#define TIMING_FUNCTIONS_H

/** Internally used floating point type */
typedef double real_t;

#if !(defined _DSPACE)
#if (defined _WIN32 || defined _WIN64) && !(defined __MINGW32__ || defined __MINGW64__)

/* Use Windows QueryPerformanceCounter for timing. */
#include <Windows.h>

/** A structure for keeping internal timer data. */
typedef struct timer
{
	LARGE_INTEGER tic;
	LARGE_INTEGER toc;
	LARGE_INTEGER freq;
} timer;


#elif (defined __APPLE__)

#include "unistd.h"
#include <mach/mach_time.h>

/** A structure for keeping internal timer data. */
typedef struct timer
{
	uint64_t tic;
	uint64_t toc;
	mach_timebase_info_data_t tinfo;
} timer;

#else

/* Use POSIX clock_gettime() for timing on non-Windows machines. */
#include <time.h>

#if __STDC_VERSION__ >= 199901L
/* C99 mode of operation. */

#include <sys/stat.h>
#include <sys/time.h>

typedef struct timer
{
	struct timeval tic;
	struct timeval toc;
} timer;

#else
/* ANSI C */

/** A structure for keeping internal timer data. */
typedef struct timer
{
	struct timespec tic;
	struct timespec toc;
} timer;

#endif /* __STDC_VERSION__ >= 199901L */

#endif /* (defined _WIN32 || defined _WIN64) */

/** A function for measurement of the current time. */
void tic(timer* t);

/** A function which returns the elapsed time. */
real_t toc(timer* t);

#endif

#endif /* TIMING_FUNCTIONS_H */
