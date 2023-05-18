/*  File: ONEview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2019
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: May 18 10:39 2023 (rd109)
 * Created: Thu Feb 21 22:40:28 2019 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "ONElib.h"

#include <assert.h>

#include <string.h>		/* strcmp etc. */
#include <stdlib.h>		/* for exit() */

static char *commandLine (int argc, char **argv)
{
  int i, totLen = 0 ;
  for (i = 0 ; i < argc ; ++i) totLen += 1 + strlen(argv[i]) ;
  char *buf = new (totLen, char) ;
  strcpy (buf, argv[0]) ;
  for (i = 1 ; i < argc ; ++i) { strcat (buf, " ") ; strcat (buf, argv[i]) ; }
  return buf ;
}

typedef struct IndexListStruct {
  I64 i0, iN ;
  struct IndexListStruct *next ;
} IndexList ;

static IndexList *parseIndexList (char *s)
{
  IndexList *ol, *ol0 = ol = new0 (1, IndexList) ;
  while (*s)
    { while (*s >= '0' && *s <= '9') ol->i0 = ol->i0*10 + (*s++ - '0') ;
      if (*s == '-')
	{ ++s ; while (*s >= '0' && *s <= '9') ol->iN = ol->iN*10 + (*s++ - '0') ;
	  if (ol->iN <= ol->i0) die ("end index %" PRId64 " <= start index %" PRId64 "", ol->iN, ol->i0) ;
	}
      else
	ol->iN = ol->i0 + 1 ;
      if (*s == ',') { ol->next = new0 (1, IndexList) ; ol = ol->next ; ++s ; }
      else if (*s) die ("unrecognised character %c at %s in object list\n", *s, s) ;
    }
  return ol0 ; 
}

static void transferLine (OneFile *vfIn, OneFile *vfOut, size_t *fieldSize)
{ memcpy (vfOut->field, vfIn->field, fieldSize[(int)vfIn->lineType]) ;
  oneWriteLine (vfOut, vfIn->lineType, oneLen(vfIn), oneString(vfIn)) ;
  char *s = oneReadComment (vfIn) ; if (s) oneWriteComment (vfOut, "%s", s) ;
}

int main (int argc, char **argv)
{
  I64 i ;
  char *fileType = 0 ;
  char *outFileName = "-" ;
  char *schemaFileName = 0 ;
  bool isNoHeader = false, isHeaderOnly = false, isBinary = false, isVerbose = false ;
  IndexList *objList = 0, *groupList = 0 ;
  
  timeUpdate (0) ;

  char *command = commandLine (argc, argv) ;
  --argc ; ++argv ;		/* drop the program name */

  if (!argc)
    { fprintf (stderr, "ONEview [options] onefile\n") ;
      fprintf (stderr, "  -t --type <abc>           file type, e.g. seq, aln - required if no header\n") ;
      fprintf (stderr, "  -S --schema <schemafile>  schema file name\n") ;
      fprintf (stderr, "  -h --noHeader             skip the header in ascii output\n") ;
      fprintf (stderr, "  -H --headerOnly           only write the header (in ascii)\n") ;
      fprintf (stderr, "  -b --binary               write in binary (default is ascii)\n") ;
      fprintf (stderr, "  -o --output <filename>    output file name (default stdout)\n") ;
      fprintf (stderr, "  -i --index x[-y](,x[-y])* write specified objects\n") ;
      fprintf (stderr, "  -g --group x[-y](,x[-y])* write specified groups\n") ;
      fprintf (stderr, "  -v --verbose              write commentary including timing\n") ;
      fprintf (stderr, "index and group only work for binary files; '-i 0-10' outputs first 10 objects\n") ;
      exit (0) ;
    }
  
  while (argc && **argv == '-')
    if (!strcmp (*argv, "-t") || !strcmp (*argv, "--type"))
      { fileType = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-S") || !strcmp (*argv, "--schema"))
      { schemaFileName = argv[1] ;
	argc -= 2 ; argv += 2 ;
      }
    else if (!strcmp (*argv, "-h") || !strcmp (*argv, "--noHeader"))
      { isNoHeader = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-H") || !strcmp (*argv, "--headerOnly"))
      { isHeaderOnly = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-b") || !strcmp (*argv, "--binary"))
      { isBinary = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-v") || !strcmp (*argv, "--verbose"))
      { isVerbose = true ; --argc ; ++argv ; }
    else if (!strcmp (*argv, "-o") || !strcmp (*argv, "--output"))
      { outFileName = argv[1] ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-i") || !strcmp (*argv, "--index"))
      { objList = parseIndexList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-g") || !strcmp (*argv, "--group"))
      { groupList = parseIndexList (argv[1]) ; argc -= 2 ; argv += 2 ; }
    else die ("unknown option %s - run without arguments to see options", *argv) ;

  if (isHeaderOnly) isBinary = false ;
    
  if (argc != 1)
    die ("need a single data one-code file as argument") ;

  OneSchema *vs = 0 ;
  if (schemaFileName && !(vs = oneSchemaCreateFromFile (schemaFileName)))
      die ("failed to read schema file %s", schemaFileName) ;
  OneFile *vfIn = oneFileOpenRead (argv[0], vs, fileType, 1) ; /* reads the header */
  if (!vfIn) die ("failed to open one file %s", argv[0]) ;

  if ((objList || groupList) && !vfIn->isBinary)
    die ("%s is ascii - you can only access objects and groups by index in binary files", argv[0]) ;
  
  OneFile *vfOut = oneFileOpenWriteFrom (outFileName, vfIn, isBinary, 1) ;
  if (!vfOut) die ("failed to open output file %s", outFileName) ;

  if (isNoHeader) vfOut->isNoAsciiHeader = true ; // will have no effect if binary
  
  if (!isHeaderOnly)
    { oneAddProvenance (vfOut, "ONEview", "0.0", command) ;
      
      static size_t fieldSize[128] ;
      for (i = 0 ; i < 128 ; ++i)
	if (vfIn->info[i]) fieldSize[i] = vfIn->info[i]->nField*sizeof(OneField) ;
      
      if (objList)
	{ while (objList)
	    { if (!oneGotoObject (vfIn, objList->i0))
		die ("can't locate to object %" PRId64 "", objList->i0 ) ;
	      if (!oneReadLine (vfIn))
		die ("can't read object %" PRId64 "", objList->i0) ;
	      while (objList->i0 < objList->iN)
		{ transferLine (vfIn, vfOut, fieldSize) ;
		  if (!oneReadLine (vfIn)) break ;
		  if (vfIn->lineType == vfIn->objectType) ++objList->i0 ;
		}
	      objList = objList->next ;
	    }
	}
      else if (groupList)
	{ while (groupList)
	    { if (!oneGotoGroup (vfIn, groupList->i0))
		die ("can't locate to group %" PRId64 "", groupList->i0 ) ;
	      if (!oneReadLine (vfIn))
		die ("can't read group %" PRId64 "", groupList->i0) ;
	      while (groupList->i0 < groupList->iN)
		{ transferLine (vfIn, vfOut, fieldSize) ;
		  if (!oneReadLine (vfIn)) break ;
		  if (vfIn->lineType == vfIn->groupType) ++groupList->i0 ;
		}
	      groupList = groupList->next ;
	    }
	}
      else
	while (oneReadLine (vfIn))
	  transferLine (vfIn, vfOut, fieldSize) ;
    }
  
  oneFileClose (vfIn) ;
  oneFileClose (vfOut) ;
  oneSchemaDestroy (vs) ;
  
  free (command) ;
  if (isVerbose)
    timeTotal (stderr) ;

  exit (0) ;
}

/********************* end of file ***********************/
