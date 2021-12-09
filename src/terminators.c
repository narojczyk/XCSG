/*
 * terminators.c
 *
 * Declarations of terminator functions (ran when exit signal is received)
 */

#include <stdio.h>
#include <stdlib.h>
#include "data.h"

const char *fmt_free_allocated_memory_msg =
    " Exitting (%d), release allocated memory for %s\n";

/*
 * free_*(exit_status, pointer)
 *
 * Cleanup functions to releace allocated memory
 */
void free_dimers  (int stat, void *ptr)
{
#ifdef VERBOSE_MDOE
  fprintf(stdout, fmt_free_allocated_memory_msg, stat, "dimers");
#endif
  free((DIM3D *) ptr ) ;
  ptr = NULL;
}

void free_spheres (int stat, void *ptr)
{
#ifdef VERBOSE_MDOE
  fprintf(stdout, fmt_free_allocated_memory_msg, stat, "spheres");
#endif
  free((SPH *) ptr ) ;
  ptr = NULL;
}

void free_slits   (int stat, void *ptr)
{
#ifdef VERBOSE_MDOE
  fprintf(stdout, fmt_free_allocated_memory_msg, stat, "slits");
#endif
  free((SLI *) ptr ) ;
  ptr = NULL;
}

void free_channels(int stat, void *ptr)
{
#ifdef VERBOSE_MDOE
  fprintf(stdout, fmt_free_allocated_memory_msg, stat, "channels");
#endif
  free((CHA *) ptr ) ;
  ptr = NULL;
}

/* vim: set tw=80 ts=2 sw=2 et: */