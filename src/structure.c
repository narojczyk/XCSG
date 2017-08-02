
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
#include "structure.h"

extern const double zero;
extern const double one;
extern const double two;

/*
 * find_valid_cluster()
 * 
 * Find a pair of molecules oriented such that they can be flipped to change
 * their orientation
 */

void find_valid_cluster(DIM3D *dim, SPH *sph, double box[3], int nd, 
                        int vclust[2])
{
  int rn_id0, rn_id1;   // array indexes generated randomly
  int a0, a1, b0, b1;   // sphere indexes of first and second dimer resp.
  int d0, nd0, nd1;     // selected dimer index and neighboring dimer ids.
  int candidate0 = 0, candidate1 = 0; // Valid configuration flags
  
  vclust[0] = -1;
  vclust[1] = -1;
  // Do the following untill valid dimer cluster is located
  do{   
    // Select a random type-1 dimer
    do{
      d0 = (int) (u_RNG() * (nd+1));
      d0 = ( d0 >= nd ) ? nd : d0;
    }while(dim[d0].type != 1);
      
    // Get index of spheres for the selected dimer
    a0 = dim[d0].sph_ind[0];
    a1 = dim[d0].sph_ind[1];
      
    // Randomly get the index of neighboring spheres
    rn_id0 = (int) (u_RNG() * 12);
    rn_id0 = ( rn_id0 >= 12 ) ? 11 : rn_id0;
    
    rn_id1 = (int) (u_RNG() * 12);
    rn_id1 = ( rn_id1 >= 12 ) ? 11 : rn_id1;
    
    // Get indexes of respective neighbors of atoms a0,a1
    b0 = sph[a0].ngb[rn_id0];
    b1 = sph[a1].ngb[rn_id1];
    
    // Clear flags
    candidate0 = candidate1 = 0;
    nd0 = nd1 = -1;
      
    // Make sure b0 is not the other atom of dimer d0 and is a type-1 sphere
    if(b0 != a1 && sph[b0].type == 1){
      nd0 = sph[b0].dim_ind;
      candidate0 = check_dimers_configuration(dim, sph, box, d0, nd0);
    }

    // Make sure b1 is not the other atom of dimer d0 and is a type-1 sphere
    if(b1 != a0 && sph[b1].type == 1){
      nd1 = sph[b1].dim_ind;
      candidate1 =  check_dimers_configuration(dim, sph, box, d0, nd1);
    }
  }while(candidate0 == 0 && candidate1 == 0);
  
  if(candidate0 != 0){
    vclust[0] = d0;
    vclust[1] = nd0;
  }else if(candidate1 != 0){
    vclust[0] = d0;
    vclust[1] = nd1;
  }
}

/*
 * check_dimers_configuration(dim,sph,d0,d1)
 * 
 * Check the orientation of d0 with relation to d1 and return the index
 * of possible transformation
 */
int check_dimers_configuration(DIM3D *dim, SPH *sph, double box[3], 
                              int d0, int d1)
{
  // Sphere ID's
  int d0atom0 = dim[d0].sph_ind[0];
  int d0atom1 = dim[d0].sph_ind[1];
  int d1atom0 = dim[d1].sph_ind[0];
  int d1atom1 = dim[d1].sph_ind[1];
  
  // 'O' times 'O', dot product of dimer orientations
  double OOdp = vdotu(dim[d0].O, dim[d1].O);
  
  // Sphere positions
  double d0r0[3] = {sph[d0atom0].r[0], sph[d0atom0].r[1], sph[d0atom0].r[2]};
  double d0r1[3] = {sph[d0atom1].r[0], sph[d0atom1].r[1], sph[d0atom1].r[2]};
  double d1r0[3] = {sph[d1atom0].r[0], sph[d1atom0].r[1], sph[d1atom0].r[2]};
  double d1r1[3] = {sph[d1atom1].r[0], sph[d1atom1].r[1], sph[d1atom1].r[2]};
  
  // Cross-dimers sphere-to-sphere distances
  double L_d0r0_d1r0 = distance(d0r0, d1r0, box);
  double L_d0r0_d1r1 = distance(d0r0, d1r1, box);
  double L_d0r1_d1r0 = distance(d0r1, d1r0, box);
  double L_d0r1_d1r1 = distance(d0r1, d1r1, box);
  
  // The sum of interatomic lengths
  double sumL = L_d0r0_d1r0 + L_d0r0_d1r1 + L_d0r1_d1r0 + L_d0r1_d1r1;

  // Out-of-plane tetrahedron configuration
  if(fabs(sumL - 4e0) < 1e-10 && fabs(OOdp) < 1e-10){
    return 69;
  }
  
  // In-plane square configuration
  if(fabs(sumL - two*(one + sqrt(two)) ) < 1e-10 && fabs(OOdp-one) < 1e-10){
    return 11;
  }
  
  // In-plane dimond configuration
  if(fabs(sumL - 3e0 - sqrt(3e0) ) < 1e-10 && fabs(OOdp-one) < 1e-10){
    return 47;
  }
  
  // If no valid configuration found return zero
  return 0;
}

/*
 * test_dimer_distribution()
 * 
 * Chceck current partition of dimer orientations
 */
void test_dimer_distribution(DIM3D *dim, int od[6], int nd)
{
  int i, j;
  int nd1=0;
  int o_ind;
  int od_local[6]={0,0,0,0,0,0};
  
  for(i=0; i<nd; i++){
    if(dim[i].type == 1){
      o_ind = check_dimer_direction(dim, i);
      if(o_ind != -1){
        od_local[o_ind]++;
        nd1++;
      }
    }
  }

  for(j=0; j<6; j++){
    od[j] = od_local[j];
  }
}

/*
 * display_dimer_distribution()
 * 
 * Display current partition of dimer orientations
 */
void display_dimer_distribution(int od[6])
{
  int j;
  
  for(j=0; j<6; j++){
    fprintf(stdout," %3d  ",od[j]);
  }
  fprintf(stdout,"\n");
}

/*
 * validate_distrib(od, nd, step)
 * 
 * checks the distribution of dimer orientations and brakes the cluster moves
 * if good-enough state is achieved. Returns '1' if the conditions are met and
 * zero otherwise.
 */
int validate_distrib(int od[6], int nd1, int step)
{
  int i;
  int level;
  int exit_code = 1; // Assume the structure is good ;)
  
  // Set exit code criteria depending on the number of dimers in the system
  if( nd1 % 6 == 0 && step < (int) 1e+8 ){
    // Perfect distribution is possible
    level = (int) (((double) nd1) / 6e0);
  }else{
    // Perfect distribution is not possible
    level = (int) (((double) nd1) / 6e0) + 1;
  }
  
  for(i=0; i<6; i++){
    exit_code *= ( od[i] <= level ) ? 1 : 0;
  }
  
  if(exit_code == 1 || step % 1000000 == 0){
    fprintf(stdout," Step %9d distribution :", step);
    for(i=0; i<6; i++){
      fprintf(stdout," %3d  ",od[i]);
    }
    fprintf(stdout,"cond.: <= %d\n",level);    
  }
 
  return exit_code;
}

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
  rn_id0 = ( rn_id0 >= 2 ) ? 1 : rn_id0;
  rn_id1 = fabs(rn_id0 - 1);
  a0 = dim[d0].sph_ind[rn_id0];
  a1 = dim[d0].sph_ind[rn_id1];
    
  rn_id0 = (int) (u_RNG() * 2);
  rn_id0 = ( rn_id0 >= 2 ) ? 1 : rn_id0;
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
}

/*
 * ziper(dim,sph,box,sph_ind)
 *
 * Run zipper on the structure, starting from sphere sph_ind, and continue
 * untill another free sphere is found.
 */
int zipper(DIM3D *dim, SPH *sph, double box[3], int nd, int sph_ind, int ms)
{
  int i=0, nseek=0;
  int step=0;
  int sngb_id, sngb_type, rand_ngb, valid_ngb, sother_id;
  int dim_ind, dngb_id, dim_t2_ind;
  int allow_type3=0;

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
    nseek=0;
    do{
      valid_ngb = 0;
      allow_type3 = 0;
      // Select neighbor index randomly
      rand_ngb = (int) (u_RNG() * 12);
      // Check not to go outside the naighbor list
      rand_ngb = (rand_ngb < 12) ? rand_ngb : 11;
      // Get neighbor id and type
      sngb_id = sph[sph_ind].ngb[rand_ngb];
      sngb_type = sph[ sngb_id ].type;
      // If enough steps passed OR if type-1 sphere is NOT on the list, 
      // allow type-3 selection
      if((sngb_type == 3) && ( (step > ms) || (nseek > 100) )){
        allow_type3 = 1;
      }
      // Select type-1 sphere or, if enough steps passed allow type-3 selection
      if( (allow_type3 == 1) || (sngb_type == 1) ) {
        valid_ngb = 1;
      }
      // Security check not to select type-2 spheres EVER!
      if(sngb_type == 2){
        valid_ngb = 0;
      }     
      // Watchdog: count attempts to select valid neighbor
      nseek++;
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
  
  // Recalculate centers of mass and orientations for type-1 diemrs
  for(i=0; i<nd; i++){
    if(dim[i].type == 1){
      update_dimer_parameters(dim, sph, box, i);
    }
  }
          
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
  DIM3D *dim, SPH *sph, double c[3], double cr, double box[3], double tr[3],
  double csd, int nd)
{
  int i;
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
  
  // Translate channel by specified vector
  for(i=0; i<3; i++){
    llc[i] += tr[i];
  }

  // Find the coordinates of the center point on the channel axis
  ccp[0] = llc[0] * (one - c[0]);
  ccp[1] = llc[1] * (one - c[1]);
  ccp[2] = llc[2] * (one - c[2]);

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
    vcrossu(p1, c, pxcd);

    // Calculate the shortest distance of the sphere center from the line
    // with reference to p1 (lower left atom)
    dist0 = vmodule(pxcd) / vmodule(c);

    // Calculate the cross product pxcd = p x cd
    vcrossu(p2, c, pxcd);

    // Calculate the shortest distance of the sphere center from the line
    // with reference to p2 (center point on the channel axis)
    dist1 = vmodule(pxcd) / vmodule(c);

    // If the sphere lies within the channel radius include the former into
    // channel and continue to the next sphere
    if(dist0 < cr || dist1 < cr){
      // Mark sphere as 'channel-sphere' (type '2')
      sph[i].type = 2;
      // ... and assign the respective channel-sphere-diameter
      sph[i].d = csd;
      // Mark dimers that cross the channel as type '2'
      dim[ sph[i].dim_ind ].type = 2;
    }
  }
}

/*
 * make_slit(dim,sph,c,ns)
 * Introduces a plane (single atomic layer) into dimer structure. The plane
 * is described by coordinate axis origin and normal vector 'c'. Dimers who's
 * atoms belong to the plane, are flagged for braking.
 * NOTE: periodic boundaries are not taken into account here
 */
void make_slit(DIM3D *dim, SPH *sph, double box[3], double thick, double os[3],
               double ssd, double nm[3], int ns){
  int i,j;
  double cd[3] = {one*nm[0], one*nm[1], one*nm[2]};
  double p[7][3];
  double dist;

  // Transform the plane vector to unit vector;
  vnorm(cd);

  // Loop over all spheres in the structure
  for(i=0; i<ns; i++){
    // Get the i'th sphere position relative to the plane inside periodic box
    p[0][0] = sph[i].r[0] - os[0];
    p[0][1] = sph[i].r[1] - os[1];
    p[0][2] = sph[i].r[2] - os[2];
    
    // ... and in the six of the sorrounding images 
    p[1][0] = sph[i].r[0] - os[0] + box[0];
    p[1][1] = sph[i].r[1] - os[1];
    p[1][2] = sph[i].r[2] - os[2];

    p[2][0] = sph[i].r[0] - os[0] - box[0];
    p[2][1] = sph[i].r[1] - os[1];
    p[2][2] = sph[i].r[2] - os[2];
    
    p[3][0] = sph[i].r[0] - os[0];
    p[3][1] = sph[i].r[1] - os[1] + box[1];
    p[3][2] = sph[i].r[2] - os[2];

    p[4][0] = sph[i].r[0] - os[0];
    p[4][1] = sph[i].r[1] - os[1] - box[1];
    p[4][2] = sph[i].r[2] - os[2];

    p[5][0] = sph[i].r[0] - os[0];
    p[5][1] = sph[i].r[1] - os[1];
    p[5][2] = sph[i].r[2] - os[2] + box[2];

    p[6][0] = sph[i].r[0] - os[0];
    p[6][1] = sph[i].r[1] - os[1];
    p[6][2] = sph[i].r[2] - os[2] - box[2];  
    
    for(j=0; j<=6; j++){
      // Calculate the dot product of the normal to the plane located at
      // axes origin and a position vector of i'th sphere
      dist = zero - cd[0]*p[j][0] - cd[1]*p[j][1] - cd[2]*p[j][2];

      // To get the distance of a sphere's center from the plane one should
      // (in principle) divide the above by the length of the plane's normal.
      // This is not done here, as the normal is unit-normed.

      // If the sphere lies within the plane, include the sphere into
      // plane and continue to the next sphere
      if(fabs(dist) < thick ){
        // Mark sphere as 'channel-sphere' (type '2')
        sph[i].type = 2;
        // Mark dimers that cross the channel as type '2'
        dim[ sph[i].dim_ind ].type = 2;
        // Escape the j-loop if sphere is found to lie on the lane
        break;
      }
      
    }
  }
}

/* vim: set tw=80 ts=2 sw=2 et: */
