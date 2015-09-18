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
 * brake_dimers(dim, sph, nd)
 * Check all dimers and brake bonds in case at least one of spheres belong
 * to channel.
 */
int brake_dimers(DIM3D *dim, SPH *sph, int nd)
{
  int i;
  int atom0, atom1;
  int broken_dimers = 0;
  for(i=0; i<nd; i++){
    
    // Check for type set by make_channel()
    if(dim[i].type == 2){
      // Get atom indexes
      atom0 = dim[i].sph_ind[0];
      atom1 = dim[i].sph_ind[1];
      
      // Delete indexes of spheres forming a dimer
//       dim[i].type = 2;
      dim[i].sph_ind[0] = -1;
      dim[i].sph_ind[1] = -1;
      
      // Set the free (non-channel) sphere (if any) to type '3'
      if(sph[atom0].type == 1){
        sph[atom0].type = 3;
      }
      if(sph[atom1].type == 1){
        sph[atom1].type = 3;
      }
      
      // Remove dimer index from the sphere data
      sph[atom0].dim_ind = -1;
      sph[atom1].dim_ind = -1;
      
      // Count broken dimers
      broken_dimers++;
      
//       printf("Dimer %4d has been destroyed\n",i);   
    }  
  }
  return broken_dimers;
}

/*
 * make_channel()
 * Calculates the distance of the sphere 'i' from the line given by vector 'c',
 * with respect to two points on the line: lower left corner and center point
 * of the channel axis (this is required because periodic boundaries may select
 * periodic image that is not the nearest to the channel axis)
 * 
 */
void make_channel(
  DIM3D *dim, SPH *sph, int c[3], double cr, double box_x, int nd)
{
  int i;
  double cd[3] = {one*c[0], one*c[1], one*c[2]};
  double llc[3]={zero,zero,zero}, ccp[3]={zero,zero,zero};
  double llc_r = -box_x/two + one;
  double p1[3], p2[3], pxcd[3], dist0, dist1;
  
  // Find coordinates of the sphere in lower-left-corner of the cube
  for(i=1; i<2*nd; i++){
    if(sph[i].r[0]<llc_r && sph[i].r[1]<llc_r && sph[i].r[2]< llc_r){
      llc[0] = sph[i].r[0];
      llc[1] = sph[i].r[1];
      llc[2] = sph[i].r[2];
    }
  }
  
  // Find the coordinates of the center point on the channel axis
  ccp[0] = llc[0] * (one - cd[0]);
  ccp[1] = llc[1] * (one - cd[1]);
  ccp[2] = llc[2] * (one - cd[2]);
  
  for(i=0; i<2*nd; i++){
    // Determine vector from point 'i' to lowe-left-corner atom of the cube
    p1[0] = llc[0] - sph[i].r[0];
    p1[1] = llc[1] - sph[i].r[1];
    p1[2] = llc[2] - sph[i].r[2];
    
    // Determine vector from point 'i' to channel-center-point of the cube
    p2[0] = ccp[0] - sph[i].r[0];
    p2[1] = ccp[1] - sph[i].r[1];
    p2[2] = ccp[2] - sph[i].r[2];

    // Apply boundary conditions    
    p1[0] = p1[0] - box_x * round( p1[0]/box_x );
    p1[1] = p1[1] - box_x * round( p1[1]/box_x );
    p1[2] = p1[2] - box_x * round( p1[2]/box_x );
    
    p2[0] = p2[0] - box_x * round( p2[0]/box_x );
    p2[1] = p2[1] - box_x * round( p2[1]/box_x );
    p2[2] = p2[2] - box_x * round( p2[2]/box_x );
    
    // Calculate the cross product pxcd = p x cd
    vcrossu(p1, cd, pxcd);
    
    // Calculate the shortest distance of the sphere center from the line
    // with reference to p1 (lower left atom)
    dist0 = vmodule(pxcd) / vmodule(cd);
    
    // Calculate the cross product pxcd = p x cd
    vcrossu(p2, cd, pxcd);
    
    // Calculate the shortest distance of the sphere center from the line
    // with reference to p2 (center point on the channel axis)
    dist1 = vmodule(pxcd) / vmodule(cd);
    
    // If the sphere lies within the channel radius include the former into
    // channel and continue to the next sphere
    if(dist0 < cr || dist1 < cr){
      // Mark sphere as 'channel-sphere' (type '2')
      sph[i].type = 2;
      // Mark dimers that cross the channel as type '2'
      dim[ sph[i].dim_ind ].type = 2;
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

    // Assign dimer index to spheres array
    sph[j  ].dim_ind = i;
    sph[j+1].dim_ind = i;

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
 * memory_clean_spheres(sph, ns)
 * memory_clean_dimers(dim, nd)
 * 
 * Assign default values to data structures' arrays
 */
void memory_clean_spheres(SPH *sph, int ns)
{
  SPH template;
  int i;

  // Default values for initial array of spheres
  template.type = 1;    // regular spheres (forming dimers)
  template.dim_ind = -1;
  template.r[0] = zero;
  template.r[1] = zero;
  template.r[2] = zero;
  template.d = one;
  
  for(i=0; i<12; i++){
    template.ngb[i] = -1;
  }
  
  // Copy default values to the array
  for(i=0; i<ns; i++){
    sph[i] = template;
  }
  
}

void memory_clean_dimers(DIM3D *dim, int nd)
{
  DIM3D template;
  int i;
  
  // Default values for initial array of dimers
  template.type = 1;    // regular dimers
  template.sph_ind[0] = -1;
  template.sph_ind[1] = -1;
  template.R[0] = zero;
  template.R[1] = zero;
  template.R[2] = zero;
  template.O[0] = zero;
  template.O[1] = zero;
  template.O[2] = zero;
  template.L = one;
  
  for(i=0; i<22; i++){
    template.ngb[i][0] = -1;
    template.ngb[i][1] = -1;
  }
  
  // Copy default values to the array
  for(i=0; i<nd; i++){
    dim[i] = template;
  }

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
