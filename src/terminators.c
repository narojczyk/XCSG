/*
 * terminators.c
 *
 * Declarations of terminator functions (ran when exit signal is received)
 */

#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "config.h"


/*
 * releace_memory(exit_status, pointer)
 *
 * Cleanup function to releace allocated memory
 */
void releace_memory(int stat, void *ptr)
{
#ifdef VERBOSE_MDOE
  const char *fmt_free_allocated_memory_msg =
    " Exitting (%d), release allocated memory for %p\n";
#endif
  const char *fmt_free_allocated_memory_failed_msg =
    " Exitting (%d), cannot free memory for null-pointer\n";

  if(ptr != NULL){
  #ifdef VERBOSE_MDOE
    fprintf(stdout, fmt_free_allocated_memory_msg, stat, ptr);
  #endif
    free(ptr);
    ptr = NULL;
  }else{
    fprintf(stderr, fmt_free_allocated_memory_failed_msg, stat);
  }
}


/* vim: set tw=80 ts=2 sw=2 et: */
