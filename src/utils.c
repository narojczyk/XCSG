/*
 * utils.c
 *
 * File defines utility functions used in the setup and control of the program
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
 * ziper(dim,sph,box,sph_ind)
 *
 * Run zipper on the structure, starting from sphere sph_ind, and continue
 * untill another free sphere is found.
 */
int zipper(DIM3D *dim, SPH *sph, double box[3], int nd, int sph_ind, int ms)
{
  int i=0;
  int step=0;
  int sngb_id, sngb_type, rand_ngb, valid_ngb, sother_id;
  int dim_ind, dngb_id, dim_t2_ind;

  if(sph[sph_ind].type == 1){
    // Get the index of dimer for type-1 sphere
    dim_ind = sph[sph_ind].dim_ind;
    // Flag dimer spheres as type-3
    sph[ dim[dim_ind].sph_ind[0] ].type = 3;
    sph[ dim[dim_ind].sph_ind[1] ].type = 3;
    sph[ dim[dim_ind].sph_ind[0] ].dim_ind = -1;
    sph[ dim[dim_ind].sph_ind[1] ].dim_ind = -1;
    // Brake dimer
    dim[dim_ind].type = 2;
    dim[dim_ind].sph_ind[0] = -1;
    dim[dim_ind].sph_ind[1] = -1;
  }else if(sph[sph_ind].type == 2){
    // Exit if the channel sphere selected as starting point
    return -1;
  }

  // Loop until the free spheres meet
  while(1){
    // Select random neighbor of sph_ind
    do{
      valid_ngb = 0;
      // Select neighbor index randomly
      rand_ngb = (int) (u_RNG() * 12);
      // Check not to go outside the naighbor list
      rand_ngb = (rand_ngb < 12) ? rand_ngb : 11;
      // Get neighbor id and type
      sngb_id = sph[sph_ind].ngb[rand_ngb];
      sngb_type = sph[ sngb_id ].type;

      // Select type-1 sphere or, if enough steps passed allow type-3 selection
      if( ((sngb_type == 3) && (step > ms)) || (sngb_type == 1) ) {
        valid_ngb = 1;
      }
      // Security check not to select type-2 spheres EVER!
      if(sngb_type == 2){
        valid_ngb = 0;
      }
    }while(!valid_ngb);

    // Brake neighboring dimer or connect free spheres and brake the cycle
    if(sngb_type == 1){
      // Get dimer index to whom sngb_id belongs to
      dngb_id = sph[sngb_id].dim_ind;

      // Get the index of the other sphere of dimer 'dngb_id'
      sother_id = (dim[dngb_id].sph_ind[0] == sngb_id)
        ? dim[dngb_id].sph_ind[1] : dim[dngb_id].sph_ind[0];

      // Flag 'sother_id' as type-3 and make dimer from the remaining spheres
      sph[sother_id].type = 3;
      sph[sother_id].dim_ind = -1;

      sph[sph_ind].type = 1;
      sph[sph_ind].dim_ind = dngb_id;

      // Ammend the indexes in data structures
      dim[dngb_id].sph_ind[0] = sph_ind;
      dim[dngb_id].sph_ind[1] = sngb_id;

      // Switch to created free sphere
      sph_ind = sother_id;

      // Count steps
      step++;

    }else if(sngb_type == 3){
      // Connect the spheres into dimers with index from the type-2 DIMERS range
      // Find a free (type-2) dimer
      i = 0;
      while(dim[i].type == 1){
        i++;
      }
      dim_t2_ind = i;

      if(dim[dim_t2_ind].type != 1 ){

        // Update sphere data
        sph[sph_ind].type = 1;
        sph[sph_ind].dim_ind = dim_t2_ind;

        sph[sngb_id].type = 1;
        sph[sngb_id].dim_ind = dim_t2_ind;

        // Update dimer data
        dim[dim_t2_ind].type = 1;
        dim[dim_t2_ind].sph_ind[0] = sph_ind;
        dim[dim_t2_ind].sph_ind[1] = sngb_id;
      }
      // Brake the master while loop and finish procedure
      break;
    }
  } // Infinite 'while' loop
  return ++step;
}


/*
 * make_slit(dim.sph,c,nd)
 * Introduces a plane (single atomic layer) into dimer structure. The plane
 * is described by coordinate axis origin and normal vector 'c'. Dimers who's
 * atoms belong to the plane, are flagged for braking.
 * NOTE: periodic boundaries are not taken into account here
 */
void make_slit(DIM3D *dim, SPH *sph, double sth, int c[3], int nd)
{
  int i;
  double cd[3] = {one*c[0], one*c[1], one*c[2]};
  double p[3];
  double dist0;

  // Transform the plane vector to unit vector;
  vnorm(cd);

  // Loop over all spheres in the structure
  for(i=0; i<2*nd; i++){
    // Get the i'th sphere position
    p[0] = sph[i].r[0];
    p[1] = sph[i].r[1];
    p[2] = sph[i].r[2];

    // Calculate the dot product of the normal to the plane located at
    // axes origin and a position vector of i'th sphere
    dist0 = zero - cd[0]*p[0] - cd[1]*p[1] - cd[2]*p[2];

    // To get the distance of a sphere's center from the plane one should
    // (in principle) divide the above by the length of the plane's normal.
    // This is not done here, as the normal is unit-normed.

// printf(" %3d %16.12lf\n",i,dist0);

    // If the sphere lies within the plane, include the sphere into
    // plane and continue to the next sphere
    if(fabs(dist0) < sth ){
      // Mark sphere as 'channel-sphere' (type '2')
      sph[i].type = 2;
      // Mark dimers that cross the channel as type '2'
      dim[ sph[i].dim_ind ].type = 2;
    }
  }
}

/*
 * find_ngb_spheres(sph,ns,box)
 * Build neighburs data for spheres.
 */
int find_ngb_spheres(SPH sph[], int ns, double box[3])
{
  int i,j,cni;
  const double lim = 5e-11;

  // Loop over all spheres
  for(i=0; i<ns; i++){
    // current neighbor index
    cni = 0;

    // Loop over all spheres other than 'i'
    for(j=0; j<ns; j++){
      if(j!=i && fabs(distance(sph[i].r,sph[j].r, box) - one) < lim ){
        sph[i].ngb[cni++] = j;
      }
    }

    // Sanity check
    if(cni != 12){
      fprintf(stderr,
              " [ %s ]: error: wrong cni value (%d != 12)\n", __func__, cni);
      return EXIT_FAILURE;
    }
  }
  return 0;
}

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
  DIM3D *dim, SPH *sph, int c[3], double cr, double box[3], int nd)
{
  int i;
  double cd[3] = {one*c[0], one*c[1], one*c[2]};
  double llc[3]={zero,zero,zero}, ccp[3]={zero,zero,zero};
  double llc_r[3] = {-box[0]/two + one, -box[1]/two + one, -box[2]/two + one};
  double p1[3], p2[3], pxcd[3], dist0, dist1;

  // Find coordinates of the sphere in lower-left-corner of the cube
  for(i=1; i<2*nd; i++){
    if(sph[i].r[0]<llc_r[0] && sph[i].r[1]<llc_r[1] && sph[i].r[2]< llc_r[2]){
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
    p1[0] = p1[0] - box[0] * round( p1[0]/box[0] );
    p1[1] = p1[1] - box[1] * round( p1[1]/box[1] );
    p1[2] = p1[2] - box[2] * round( p1[2]/box[2] );

    p2[0] = p2[0] - box[0] * round( p2[0]/box[0] );
    p2[1] = p2[1] - box[1] * round( p2[1]/box[1] );
    p2[2] = p2[2] - box[2] * round( p2[2]/box[2] );

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
 * update_dimer_parameters(dim,sph,box,d)
 * 
 * Update the orientation and the center of mass coordinates for dimer d,
 * based on the spheres' positions
 */
void update_dimer_parameters(DIM3D *dim, SPH *sph, double box[3], int d)
{
  int atom0 = dim[d].sph_ind[0];
  int atom1 = dim[d].sph_ind[1];
  double O[3] = {zero,zero,zero};
  double R[3] = {zero,zero,zero};
  
  // Calculate current dimer orientation
  vector(sph[atom1].r, sph[atom0].r, O, box);
  
  // Unit norm the calculated vector
  vnorm(O);
  
  // Calculate the center of mass
  R[0] = sph[atom1].r[0] + O[0]/two;
  R[1] = sph[atom1].r[1] + O[1]/two;
  R[2] = sph[atom1].r[2] + O[2]/two;
  
  // Store calculations in the data structure
  dim[d].R[0] = R[0];
  dim[d].R[1] = R[1];
  dim[d].R[2] = R[2];
  dim[d].O[0] = O[0];
  dim[d].O[1] = O[1];
  dim[d].O[2] = O[2];
}
/*
 * update_sphere_positions(dim, sph, d)
 * Updates positions of spheres for the dimer 'd' with regard to the periodic
 * boundaries (the simulation box is assumed to be -L/2;L/2)
 */
void update_sphere_positions(DIM3D *dim, SPH *sph, double box[3], int d)
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
    sp_r0[j] = sp_r0[j] - box[j] * round( sp_r0[j]/box[j] );
    sp_r1[j] = sp_r1[j] - box[j] * round( sp_r1[j]/box[j] );

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
 * sph_set_fcc(sph, ns, fcc)
 * Set fcc structure of ns spheres in a cubic system of fcc cells at the
 * edge. The edge length is assumed sqrt(2)
 */
int sph_set_fcc( SPH *sph, int ns, int fcc[3])
{
  int i, x=0, y=0, z=0;
  double cell_edge = sqrt(two);
  double cell_edge_half = cell_edge/two;

  for(i=0; i<ns; i+=4){
    sph[i  ].r[0] = - fcc[0] * cell_edge_half + x * cell_edge;
    sph[i  ].r[1] = - fcc[1] * cell_edge_half + y * cell_edge;
    sph[i  ].r[2] = - fcc[2] * cell_edge_half + z * cell_edge;

    sph[i+1].r[0] = - fcc[0] * cell_edge_half + x * cell_edge + cell_edge_half;
    sph[i+1].r[1] = - fcc[1] * cell_edge_half + y * cell_edge;
    sph[i+1].r[2] = - fcc[2] * cell_edge_half + z * cell_edge + cell_edge_half;

    sph[i+2].r[0] = - fcc[0] * cell_edge_half + x * cell_edge;
    sph[i+2].r[1] = - fcc[1] * cell_edge_half + y * cell_edge + cell_edge_half;
    sph[i+2].r[2] = - fcc[2] * cell_edge_half + z * cell_edge + cell_edge_half;

    sph[i+3].r[0] = - fcc[0] * cell_edge_half + x * cell_edge + cell_edge_half;
    sph[i+3].r[1] = - fcc[1] * cell_edge_half + y * cell_edge + cell_edge_half;
    sph[i+3].r[2] = - fcc[2] * cell_edge_half + z * cell_edge;

    sph[i  ].d = one;
    sph[i+1].d = one;
    sph[i+2].d = one;
    sph[i+3].d = one;

    // Increment cell in x direction
    x++;
    // Increment cell in y direction
    if( x == fcc[0]){
      x = 0;
      y++;
    };
    // increment cell in z direction
    if( y == fcc[1]){
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
 * u_RNG()
 *
 * Returns a pseudo random number uniformly distributed [0:1]
 */
double u_RNG()
{
#ifdef USE_64BIT_MT19937
  return genrand64_real1();
#endif

#ifdef USE_32BIT_MT19937
  return genrand_real1();
#endif

#ifdef USE_DRAND48
  return drand48();
#endif

  // In case of failure of the above
  return zero;
}

/*
 * init_RNG( s )
 *
 * This function initiates MT19937 rng (either 32 or 64bit, as selected
 * with includes) with a unsigned long int 's' seed.
 */
void init_RNG(unsigned long int s)
{
  #ifdef USE_64BIT_MT19937
    init_genrand64( s );
  #endif

  #ifdef USE_32BIT_MT19937
    init_genrand( s );
  #endif

  #ifdef USE_DRAND48
    srand48( s );
  #endif
}

/* vim: set tw=80 ts=2 sw=2 et: */
