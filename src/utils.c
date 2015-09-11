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
#include "algebra.h"

extern const double zero;
extern const double one;
extern const double two;
extern const double pi;

/*
 * make_chanel()
 * 
 */
void make_chanel(DIM3D *dim, SPH *sph, int c[3], double box_x, int nd)
{
  int i;
  int is=-1;
  double chan[3] = {one*c[0], one*c[1], one*c[2]};
  double llc[3]={zero,zero,zero}, llc_r = -box_x/two + one;
  double p[3], pc[3], dist;
  
  /*llc[0] = llc_r+1e-15;
  llc[1] = llc_r+1e-15;
  llc[2] = llc_r+1e-15;*/
  
  printf("\n\t[%s]\n",__func__);
  printf("make channel in direction: %d %d %d\n",c[0], c[1], c[2]);
  
  for(i=1;i<2*nd;i++){
    if(sph[i].r[0]<llc_r && sph[i].r[1]<llc_r && sph[i].r[2]< llc_r){
      is=i; // obsolete
      llc[0] = sph[i].r[0];
      llc[1] = sph[i].r[1];
      llc[2] = sph[i].r[2];
    }
  }
  if(is != -1){
    printf("lowe left corner: (%d) %lf %lf %lf\n",is, 
           sph[is].r[0], sph[is].r[1], sph[is].r[2]);
    printf("shift by: %lf %lf %lf\n", llc[0], llc[1], llc[2]);
  }
  
  for(i=0;i<2*nd;i++){
    // Determine vector from point 'i' to lowe-left-corner atom of the cube
    p[0] = llc[0] - sph[i].r[0];
    p[1] = llc[1] - sph[i].r[1];
    p[2] = llc[2] - sph[i].r[2];
    
    if(sph[i].r[2]<-3.8){
      printf("\n%3d % .6lf % .6lf\t(r)\n",i,sph[i].r[0],sph[i].r[1]);
      printf("%3d % .6lf % .6lf\t(p)\n",i,p[0],p[1]);
    }
    
    // Boundary conditions
    
    p[0] = p[0] - box_x * round( p[0]/box_x );
    p[1] = p[1] - box_x * round( p[1]/box_x );
    p[2] = p[2] - box_x * round( p[2]/box_x );

    if(sph[i].r[2]<-3.8){
      printf("%3d % .6lf % .6lf\t(p-pbc)\n",i,p[0],p[1]);
    }
    
    pc[0]=zero; pc[1]=zero; pc[2]=zero;
    vcrossu(p, chan, pc);
    
//     printf("%4d p : % .6lf % .6lf % .6lf\n",i,p[0],p[1],p[2]);
//     printf("%4d pc: % .6lf % .6lf % .6lf\n",i,pc[0],pc[1],pc[2]);
    
    dist = vmodule(pc) / vmodule(chan);
    if(sph[i].r[2]<-3.8){
      printf("%3d % .6lf\t\t(dist)\n",i,dist);
    }
    if(dist<1.65){
      sph[i].type=1;
    }
  }
  
  
  
  
}

/*
 * update_sphere_positions(dim, sph, d)
 * Updates positions of spheres for the dimer 'd' with regard to the periodic 
 * boundaries (the simulation box is assumed to be -L/2;L/2)
 */
void update_sphere_positions(DIM3D *dim, SPH *sph, double box_x, int d)
{
  int j;
  int atom0 = dim[d].sph_ind[0];
  int atom1 = dim[d].sph_ind[1];
  double sp_r0[3], sp_r1[3];
  
  for(j=0;j<3;j++){
    // Generate position component for the first atom
    sp_r0[j] = dim[d].R[j] + dim[d].O[j] * dim[d].L / two;
    
    // Generate position component for the second atom
    sp_r1[j] = dim[d].R[j] - dim[d].O[j] * dim[d].L / two;
    
    // Apply paeriodic boundaries
    sp_r0[j] = sp_r0[j] - box_x * round( sp_r0[j]/box_x );
    sp_r1[j] = sp_r1[j] - box_x * round( sp_r1[j]/box_x );
    
    // Assign values to sphere array
    sph[atom0].r[j] = sp_r0[j];
    sph[atom1].r[j] = sp_r1[j];
  }

}

/*
 * bind_spheres_to_dimers(dim, sph, nd)
 * Set links between dimers and spheres data structures (i.e. bind which spheres
 * beong to which dimmers and vice vesra)
 */
void bind_spheres_to_dimers(DIM3D *dim, SPH *sph, int nd)
{
  int i=0, j=0;
  
  while(i<nd){
    // Assign sphere indexes to dimers array
    dim[i].sph_ind[0] = j;
    dim[i].sph_ind[1] = j+1;
//     printf("%3d (%3d %3d)\n",i,dim[i].sph_ind[0], dim[i].sph_ind[1]);
    // Assign dimer index to spheres array
    sph[j  ].dim_ind = i;
    sph[j+1].dim_ind = i;
    
    sph[j  ].d = one;
    sph[j+1].d = one;
//     printf("%3d (%3d)\n",j,sph[j].dim_ind);
//     printf("%3d (%3d)\n",j+1,sph[j+1].dim_ind);
    // Incremend counters
    i++;
    j=j+2;
  }  
}

/*
 * sph_set_fcc(sph, ns, fcc_x)
 * Set fcc structure of ns spheres in a cubic system of fcc_x cells at the
 * edge. The edge length is assumed sqrt(2)
 */
int sph_set_fcc( SPH *sph, int ns, int fcc_x)
{
  int i, x=0, y=0, z=0;
  double cell_edge = sqrt(two);
  double cell_edge_half = cell_edge/two;

  for(i=0; i<ns; i+=4){
    sph[i  ].r[0] = - fcc_x * cell_edge_half + x * cell_edge;
    sph[i  ].r[1] = - fcc_x * cell_edge_half + y * cell_edge;
    sph[i  ].r[2] = - fcc_x * cell_edge_half + z * cell_edge;
    
    sph[i+1].r[0] = - fcc_x * cell_edge_half + x * cell_edge + cell_edge_half;
    sph[i+1].r[1] = - fcc_x * cell_edge_half + y * cell_edge;
    sph[i+1].r[2] = - fcc_x * cell_edge_half + z * cell_edge + cell_edge_half;
    
    sph[i+2].r[0] = - fcc_x * cell_edge_half + x * cell_edge;
    sph[i+2].r[1] = - fcc_x * cell_edge_half + y * cell_edge + cell_edge_half;
    sph[i+2].r[2] = - fcc_x * cell_edge_half + z * cell_edge + cell_edge_half;
    
    sph[i+3].r[0] = - fcc_x * cell_edge_half + x * cell_edge + cell_edge_half;
    sph[i+3].r[1] = - fcc_x * cell_edge_half + y * cell_edge + cell_edge_half;
    sph[i+3].r[2] = - fcc_x * cell_edge_half + z * cell_edge;
    
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
