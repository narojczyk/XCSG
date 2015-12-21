
/*
 * structure.c
 *
 * File defines functions to preform operations on the DC structure
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "data.h"
#include "utils.h"
#include "algebra.h"

extern const double zero;
extern const double one;
extern const double two;
extern const double pi;

/*
 * flip_dimers(dim, sph, box, od, d0, d1)
 * 
 * Change configuration of two specified dimers d0, d1 and update orientation
 * table
 */

void flip_dimers(DIM3D *dim, SPH *sph, double box[3], int od[6], int d0, int d1)
{
  int rn_id0, rn_id1;   // array indexes generated randomly
  int a0, a1, b0, b1;   // sphere indexes of first and second dimer resp.
  int o_id0 = check_dimer_direction(dim, d0);
  int o_id1 = check_dimer_direction(dim, d1);
  double dist0, dist1;  // interatomic distances
  
  // Decrement appropriate orientation buckets
  if(o_id0 != -1){
    od[o_id0]--;
  }else{
    fprintf(stderr, "[%s] Illegal orientation index (%d) of dimer %d\n",
      __func__, o_id0, d0);    
  }
  
  if(o_id1 != -1){
    od[o_id1]--;
  }else{
    fprintf(stderr, "[%s] Illegal orientation index (%d) of dimer %d\n",
      __func__, o_id1, d1);    
  }
  
  // Assign corresponding sphere indexes randomly to fliping handler
  rn_id0 = (int) (u_RNG() * 2);
  rn_id0 = ( rn_id0 >= 2 ) ? 2 : rn_id0;
  rn_id1 = fabs(rn_id0 - 1);
  a0 = dim[d0].sph_ind[rn_id0];
  a1 = dim[d0].sph_ind[rn_id1];
    
  rn_id0 = (int) (u_RNG() * 2);
  rn_id0 = ( rn_id0 >= 2 ) ? 2 : rn_id0;
  rn_id1 = fabs(rn_id0 - 1);
  b0 = dim[d1].sph_ind[rn_id0];
  b1 = dim[d1].sph_ind[rn_id1];
  
  // Calculate the interatomic dostances for pairs a0-b0 and a1-b1
  dist0 = distance(sph[a0].r, sph[b0].r, box);  
  dist1 = distance(sph[a1].r, sph[b1].r, box);
  
  // Assign paired spheres to respective dimers if BOTH distances are equal
  // to $\sigm$ (i.e. spheres lie on neighboring lattice sites).
  // In the other case mismatch pairs
  if(fabs(dist0 - one) < 1e-10 && fabs(dist1 - one) < 1e-10){
    // Forming dimer a0-b0 and a1-b1
    dim[d0].sph_ind[0] = a0; 
    dim[d0].sph_ind[1] = b0;
    sph[a0].dim_ind = d0;
    sph[b0].dim_ind = d0;
    
    dim[d1].sph_ind[0] = a1;
    dim[d1].sph_ind[1] = b1;
    sph[a1].dim_ind = d1;
    sph[b1].dim_ind = d1;
  }else{
    // Forming dimer a0-b1 and a1-b0
    dim[d0].sph_ind[0] = a0;
    dim[d0].sph_ind[1] = b1;
    sph[a0].dim_ind = d0;
    sph[b1].dim_ind = d0;
      
    dim[d1].sph_ind[0] = a1;
    dim[d1].sph_ind[1] = b0;
    sph[a1].dim_ind = d1;
    sph[b0].dim_ind = d1;
  }
    
  // Update orientation and the center of mass for fliped dimers
  update_dimer_parameters(dim, sph, box, d0);
  update_dimer_parameters(dim, sph, box, d1);  
  
  // Increment appropriate orientation buckets
  o_id0 = check_dimer_direction(dim, d0);
  o_id1 = check_dimer_direction(dim, d1);  
  
  if(o_id0 != -1){
    od[o_id0]++;
  }else{
    fprintf(stderr," [%s] Illegal orientation index (%d) of fliped dimer %d\n",
      __func__, o_id0, d0);    
  }
  
  if(o_id1 != -1){
    od[o_id1]++;
  }else{
    fprintf(stderr," [%s] Illegal orientation index (%d) of fliped dimer %d\n",
      __func__, o_id1, d1);    
  }
  
  /*        
      printf("After fliping\n");
 printf("dimer %3d <% .5lf % .5lf % .5lf> (%3d %3d)\n",
        d0,dim[d0].O[0],dim[d0].O[1],dim[d0].O[2],
        dim[d0].sph_ind[0],dim[d0].sph_ind[1]
       );
 printf("dimer %3d <% .5lf % .5lf % .5lf> (%3d %3d)\n",
        dngb1,dim[d1].O[0],dim[d1].O[1],dim[d1].O[2],
        dim[d1].sph_ind[0],dim[d1].sph_ind[1]
       );
  */
  
}

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

/* vim: set tw=80 ts=2 sw=2 et: */
