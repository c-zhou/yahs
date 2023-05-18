/*  File: utils.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) R Durbin, 1996
 *-------------------------------------------------------------------
 * Description: core utility functions
 * Exported functions:
 * HISTORY:
 * Last edited: May 15 14:44 2023 (rd109)
 * * Feb 22 14:52 2019 (rd109): added fzopen()
 * Created: Thu Aug 15 18:32:26 1996 (rd)
 *-------------------------------------------------------------------
 */

#ifdef LINUX
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdarg.h>
#include <ctype.h>
#include "utils.h"

void die (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "FATAL ERROR: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;

  exit (-1) ;
}

void warn (char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "WARNING: ") ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
}

static char* commandLine = 0 ;

void storeCommandLine (int argc, char **argv)
{
  int i, totLen = 0 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  if (commandLine) free (commandLine) ;
  commandLine = new (totLen, char) ;
  strcpy (commandLine, argv[0]) ;
  for (i = 1 ; i < argc ; ++i) { strcat (commandLine, " ") ; strcat (commandLine, argv[i]) ; }
}

char *getCommandLine (void) { return commandLine ; }

long totalAllocated = 0 ;

void *myalloc (size_t size)
{
  void *p = (void*) malloc (size) ;
  if (!p) die ("myalloc failure requesting %d bytes", size) ;
  totalAllocated += size ;
  return p ;
}

void *mycalloc (size_t number, size_t size)
{
  void *p = (void*) calloc (number, size) ;
  if (!p) die ("mycalloc failure requesting %d objects of size %d", number, size) ;
  totalAllocated += size*number ;
  return p ;
}

char *fgetword (FILE *f)
{
  int n = 0 ;
  static char *buf = 0 ;
  int bufSize = 64 ;
  char *cp ;
  if (!buf) buf = (char*) myalloc (bufSize) ;
  cp = buf ;
  while (!feof (f) && (*cp = getc (f)))
    if (isgraph(*cp) && !isspace(*cp))
      { if (++n >= bufSize)
	  { bufSize *= 2 ;
	    if (!(buf = (char*) realloc (buf, bufSize)))
	      die ("fgetword realloc failure requesting %d bytes", bufSize) ;
	    cp = &buf[n] ;
	  }
	else
	  ++cp ;
      }
    else
      { while (*cp != '\n' && !feof(f) && (isspace(*cp) || !isgraph(*cp))) *cp = getc (f) ;
	ungetc (*cp, f) ;
	break ;
      }
  *cp = 0 ;
  return buf ;
}

#define WITH_ZLIB
#ifdef WITH_ZLIB
#include <zlib.h>
#endif

FILE *fzopen(const char *path, const char *mode)
{  /* very cool from https://stackoverflow.com/users/3306211/fernando-mut */
#if defined(WITH_ZLIB) && (defined(MACOS) || defined(LINUX))
  gzFile zfp = 0 ;			/* fernando said *zfp - makes me worry.... */

  if (strlen(path) > 3 && !strcmp(&path[strlen(path)-3], ".gz")) // only gzopen on .gz files
    zfp = gzopen(path,mode) ; // try gzopen

  if (!zfp) return fopen(path,mode);

  /* open file pointer */
#ifdef MACOS
  return funopen(zfp,
                 (int(*)(void*,char*,int))gzread,
                 (int(*)(void*,const char*,int))gzwrite,
                 (fpos_t(*)(void*,fpos_t,int))gzseek,
                 (int(*)(void*))gzclose) ;
#else
  { cookie_io_functions_t io_funcs ;
    io_funcs.read = gzread ;
    io_funcs.write = gzwrite ;
    io_funcs.seek = gzseek ;
    io_funcs.close = gzclose ;
    return fopencookie (zfp, mode, io_funcs) ;
  }
#endif

#else
  return fopen(path,mode);
#endif
}

char* fnameTag (char* root, char* tag)
{
  char *fileName = new (strlen (root) + strlen (tag) + 2, char) ;
  strcpy (fileName, root) ;
  strcat (fileName, ".") ;
  strcat (fileName, tag) ;
  return fileName ;
}

FILE *fopenTag (char* root, char* tag, char* mode)
{
  char *fileName = fnameTag (root, tag) ;
  FILE *f = fzopen (fileName, mode) ;
  free (fileName) ;
  return f ;
}

/***************** rusage for timing information ******************/

#include <sys/resource.h>
#ifndef RUSAGE_SELF     /* to prevent "RUSAGE_SELF redefined" gcc warning, fixme if this is more intricate */
#define RUSAGE_SELF 0
#endif

#ifdef RUSAGE_STRUCTURE_DEFINITIONS
struct rusage {
  struct timeval ru_utime; /* user time used */
  struct timeval ru_stime; /* system time used */
  long ru_maxrss;          /* integral max resident set size */
  long ru_ixrss;           /* integral shared text memory size */
  long ru_idrss;           /* integral unshared data size */
  long ru_isrss;           /* integral unshared stack size */
  long ru_minflt;          /* page reclaims */
  long ru_majflt;          /* page faults */
  long ru_nswap;           /* swaps */
  long ru_inblock;         /* block input operations */
  long ru_oublock;         /* block output operations */
  long ru_msgsnd;          /* messages sent */
  long ru_msgrcv;          /* messages received */
  long ru_nsignals;        /* signals received */
  long ru_nvcsw;           /* voluntary context switches */
  long ru_nivcsw;          /* involuntary context switches */
};

struct timeval {
  time_t       tv_sec;   /* seconds since Jan. 1, 1970 */
  suseconds_t  tv_usec;  /* and microseconds */
} ;
#endif /* RUSAGE STRUCTURE_DEFINITIONS */

static struct rusage rOld, rFirst ;

void timeUpdate (FILE *f)
{
  static bool isFirst = 1 ;
  struct rusage rNew ;
  int secs, usecs ;

  getrusage (RUSAGE_SELF, &rNew) ;
  if (!isFirst)
    { secs = rNew.ru_utime.tv_sec - rOld.ru_utime.tv_sec ;
      usecs =  rNew.ru_utime.tv_usec - rOld.ru_utime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "user\t%d.%06d", secs, usecs) ;
      secs = rNew.ru_stime.tv_sec - rOld.ru_stime.tv_sec ;
      usecs =  rNew.ru_stime.tv_usec - rOld.ru_stime.tv_usec ;
      if (usecs < 0) { usecs += 1000000 ; secs -= 1 ; }
      fprintf (f, "\tsystem\t%d.%06d", secs, usecs) ;
      fprintf (f, "\tmax_RSS\t%ld", rNew.ru_maxrss - rOld.ru_maxrss) ;
      fprintf (f, "\tmemory\t%li", totalAllocated) ;   
      fputc ('\n', f) ;
    }
  else
    { rFirst = rNew ;
      isFirst = false ;
    }

  rOld = rNew ;
}

void timeTotal (FILE *f) { rOld = rFirst ; timeUpdate (f) ; }

/********************* end of file ***********************/
