/*
 * io.c
 *
 * File defines I/O functions used in the program
 */

#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
// #include <string.h>
// #include "config.h"
#include "data.h"
#include "io.h"


/*
 * export_spheres(f,sp,ns)
 * Export positions and other data for spheres
 */
int export_spheres(FILE *file, SPH *sp_tab, int ns)
{
  const char *exp_form = "%5d % .16le % .16le % .16le %.16le %d\n";
  int i;

  for(i=0; i<ns; i++){
    if(fprintf(file, exp_form, i, 
        sp_tab[i].r[0], sp_tab[i].r[1], sp_tab[i].r[2], 
        sp_tab[i].d, sp_tab[i].type ) == EOF){
      fprintf(stderr,"  [%s]: error: exporting positions failed\n", __func__);
      return EXIT_FAILURE;
    }
  }
  return 0;
}