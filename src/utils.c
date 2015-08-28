/*
 * utils.c
 *
 * File defines utility functions used in the setup and control of the program
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <string.h>
// #include "config.h"
#include "data.h"

extern const double zero;
extern const double one;
extern const double two;
extern const double pi;

/*
 * sph_set_fcc(sp_tab, ns, fcc_x)
 * Set fcc structure of ns spheres in a cubic system of fcc_x cells at the
 * edge. The edge length is assumed sqrt(2)
 */
int sph_set_fcc( SPH *sp_tab, int ns, int fcc_x)
{
  int i, x=0, y=0, z=0;
  double kr = sqrt(two);

  for(i=0; i<ns; i+=4){
    sp_tab[i  ].r[0] = x * kr;
    sp_tab[i  ].r[1] = y * kr;
    sp_tab[i  ].r[2] = z * kr;
    
    sp_tab[i+1].r[0] = x * kr +  kr/2;
    sp_tab[i+1].r[1] = y * kr;
    sp_tab[i+1].r[2] = z * kr +  kr/2;
    
    sp_tab[i+2].r[0] = x * kr;
    sp_tab[i+2].r[1] = y * kr +  kr/2;
    sp_tab[i+2].r[2] = z * kr +  kr/2;
    
    sp_tab[i+3].r[0] = x * kr +  kr/2;
    sp_tab[i+3].r[1] = y * kr +  kr/2;
    sp_tab[i+3].r[2] = z * kr;

    // Increment cell in x direction
    x++;
    // Increment cell in y direction
    if( x == fcc_x){
      x = 0;
      y++;
    };
    // increment cell in z direction
    if( y == fcc_x){
      y = 0;
      z++;
    };
  }
  // Security check
  if( (x+1)*(y+1)*(z+1)*4 > ns ){
    fprintf(stderr,"\n [%s] structure index out of range\n",__func__);
    fprintf(stderr," index %d > Ns\n",(x+1)*(y+1)*(z+1)*4);
    return EXIT_FAILURE;
  }
  
  return 0;
}

/* vim: set tw=80 ts=2 sw=2 et: */
