/*
 * utils.c
 *
 * File defines utility functions used in the setup and control of the program
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <string.h>
#include "config.h"

#ifdef USE_64BIT_MT19937
  #include "mt19937_64.h"
#else
  #include "mt19937.h"
#endif

#include "data.h"
#include "utils.h"

extern const double zero;
extern const double one;
extern const double two;
extern const double pi;

/*
 * sph_set_fcc(sph, ns, fcc_x)
 * Set fcc structure of ns spheres in a cubic system of fcc_x cells at the
 * edge. The edge length is assumed sqrt(2)
 */
int sph_set_fcc( SPH *sph, int ns, int fcc_x)
{
  int i, x=0, y=0, z=0;
  double box_edge = sqrt(two);
  double box_edge_half = box_edge/two;

  for(i=0; i<ns; i+=4){
    sph[i  ].r[0] = - fcc_x * box_edge_half + x * box_edge;
    sph[i  ].r[1] = - fcc_x * box_edge_half + y * box_edge;
    sph[i  ].r[2] = - fcc_x * box_edge_half + z * box_edge;
    
    sph[i+1].r[0] = - fcc_x * box_edge_half + x * box_edge + box_edge_half;
    sph[i+1].r[1] = - fcc_x * box_edge_half + y * box_edge;
    sph[i+1].r[2] = - fcc_x * box_edge_half + z * box_edge + box_edge_half;
    
    sph[i+2].r[0] = - fcc_x * box_edge_half + x * box_edge;
    sph[i+2].r[1] = - fcc_x * box_edge_half + y * box_edge + box_edge_half;
    sph[i+2].r[2] = - fcc_x * box_edge_half + z * box_edge + box_edge_half;
    
    sph[i+3].r[0] = - fcc_x * box_edge_half + x * box_edge + box_edge_half;
    sph[i+3].r[1] = - fcc_x * box_edge_half + y * box_edge + box_edge_half;
    sph[i+3].r[2] = - fcc_x * box_edge_half + z * box_edge;
    
    sph[i  ].d = one;
    sph[i+1].d = one;
    sph[i+2].d = one;
    sph[i+3].d = one;

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

/*
 * init_MT19937( s )
 *
 * This function initiates MT19937 rng (either 32 or 64bit, as selected
 * with includes) with a unsigned long int 's' seed.
 */
void init_MT19937(unsigned long int s)
{
  #ifdef USE_64BIT_MT19937
    init_genrand64( s );
  #else
    init_genrand( s );
  #endif
}

/* vim: set tw=80 ts=2 sw=2 et: */
