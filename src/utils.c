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

int draw_ngb_sphere_typeX(SPH *sph, int x, int sph_ind)
{
  int nseek = 0;
  int sngb_id, sngb_type, rand_ngb, valid_ngb;

  // Check if sph_ind is a valid array index
  if(sph_ind == -1){
    return -1;
  }

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

    // Select type-x sphere
    if(sngb_type == x ) {
      valid_ngb = 1;
    }
    // Security check not to select type-2 spheres EVER!
    if(sngb_type == 2){
      valid_ngb = 0;
    }
    // Watchdog: count attempts to select valid neighbor
    nseek++;
    if(nseek > 500){
      fprintf(stderr,
              " [%s] error: Cannot locate type %d neighbour of %d sphere\n",
              __func__, x, sph_ind);
      fprintf(stderr," [%s] terminated by Watchdog\n", __func__);
      return -1;
    }
  }while(!valid_ngb);

  return sngb_id;
}

/*
 * draw_sphere_typeX(sph,x,ns)
 *
 * Randomly select a sphere of type x
 * BEWARE of the infinite loop possibility.
 */
int draw_sphere_typeX(SPH *sph, int x, int ns)
{
  int sp_id = (int) (u_RNG() * ns);

  while(1){
    // Check not to go outside the sphere tab
    sp_id = (sp_id < ns) ? sp_id : ns - 1;
    // Check if the type is ok
    if(sph[sp_id].type == x){
      return sp_id;
    }else{
      // Draw another index in type is not correct.
      sp_id = (int) (u_RNG() * ns);
    }
  }
}

/*
 * find_critical_FS(sph,ns)
 *
 * Find a free type-3 sphere with the lowest count of type-3 neigbours
 */

int find_critical_FS(SPH *sph, int ns)
{
  int i, ic;
  int min_ic = 100;
  int crit_fs = -1;

  for(i=0; i<ns; i++){
    if(sph[i].type == 3){
      // clear type-3 spheres count for i'th sphere
      ic = count_typeX_sp_neighbours(sph, 3, i)  ;

      // remember 'i' index if calculated ic is minimal
      if(ic > 0 && ic < min_ic){
        min_ic = ic;
        crit_fs = i;
      }
    }
  }

  // Return the index of the selected sphere (or error code if none found)
  return crit_fs;
}

/*
 * count_typeX_sp_neighbours(sph,x,id)
 *
 * check the neighbours list of sphere 'id' and count all neighbours of type x
 */
int count_typeX_sp_neighbours(SPH *sph, int x, int id)
{
  int ic = 0, i;

  // Check if id is valid array index
  if(id == -1){
    return -1;
  }

  // count possible neighbours to connect.
  for(i=0; i<12; i++){
    ic += (sph[sph[id].ngb[i]].type == x) ? 1 : 0 ;
  }
  return ic;
}

int count_typeX_dimers(DIM3D *dim, int x, int nd)
{
  int i, count=0;

  for(i=0; i<nd; i++){
    if(dim[i].type == x){
      count++;
    }
  }
  return count;
}

int count_typeX_spheres(SPH *sph, int x, int ns)
{
  int i, count=0;

  for(i=0; i<ns; i++){
    if(sph[i].type == x){
      count++;
    }
  }
  return count;
}

int check_dimer_direction(DIM3D *dim, int i)
{
  int j;
  double fcc_dir[6][3];

  fcc_dir[0][0] =  one;
  fcc_dir[0][1] =  one;
  fcc_dir[0][2] = zero;

  fcc_dir[1][0] = -one;
  fcc_dir[1][1] =  one;
  fcc_dir[1][2] = zero;

  fcc_dir[2][0] =  one;
  fcc_dir[2][1] = zero;
  fcc_dir[2][2] =  one;

  fcc_dir[3][0] = -one;
  fcc_dir[3][1] = zero;
  fcc_dir[3][2] =  one;

  fcc_dir[4][0] = zero;
  fcc_dir[4][1] =  one;
  fcc_dir[4][2] =  one;

  fcc_dir[5][0] = zero;
  fcc_dir[5][1] = -one;
  fcc_dir[5][2] =  one;

  if(dim[i].type == 1){
    for(j=0; j<6; j++){
      if( fabs(fabs(vdotu(fcc_dir[j], dim[i].O)) -sqrt(two)) < 1e-10 ){
        return j;
      }
    }
  }
  fprintf(stderr,"\n [%s] error: unrecognized dimer orientation\n",__func__);
  fprintf(stderr," dimer %d O: %.16le %.16le %.16le\n",
          i, dim[i].O[0], dim[i].O[1], dim[i].O[2]);
  return -1;
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

    // Check for type set by make_channel() or make_slit()
    if(dim[i].type == 2){
      // Get atom indexes
      atom0 = dim[i].sph_ind[0];
      atom1 = dim[i].sph_ind[1];

      // Delete indexes of spheres forming a dimer
      dim[i].sph_ind[0] = -1;
      dim[i].sph_ind[1] = -1;

      // Set the free (non-channel) sphere (if any) to type '3'
      if(atom0 >= 0 && sph[atom0].type == 1){
        sph[atom0].type = 3;
      }
      if(atom1 >= 0 && sph[atom1].type == 1){
        sph[atom1].type = 3;
      }

      // Remove dimer index from the sphere data
      if(atom0 >= 0){
        sph[atom0].dim_ind = -1;
      }
      if(atom1 >= 0){
        sph[atom1].dim_ind = -1;
      }

      // Count broken dimers
      broken_dimers++;
    }
  }
  return broken_dimers;
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
  double com[3]={zero,zero,zero};

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

    // set default type as free spheres
    sph[i  ].type = 3;
    sph[i+1].type = 3;
    sph[i+2].type = 3;
    sph[i+3].type = 3;

    // clear dimer indexes
    sph[i  ].dim_ind = -1;
    sph[i+1].dim_ind = -1;
    sph[i+2].dim_ind = -1;
    sph[i+3].dim_ind = -1;

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
  // Calculate the center of mass of the new structure
  for(i=0; i<ns; i++){
    com[0] += sph[i].r[0];
    com[1] += sph[i].r[1];
    com[2] += sph[i].r[2];

  }
  // Move the center of mass of the new structure to point 0
  for(i=0; i<ns; i++){
    sph[i].r[0] -= com[0]/ns;
    sph[i].r[1] -= com[1]/ns;
    sph[i].r[2] -= com[2]/ns;
  }

  return EXIT_SUCCESS;
}

/*
 * memory_clean_spheres(sph, ns)
 * memory_clean_dimers(dim, nd)
 * memory_clean_slits(sl, nsl)
 * memory_clean_channels(cha, nch)
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
  template.type = 2;    // regular dimers
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

void memory_clean_slits(SLI *sli, int nsl)
{
  SLI template;
  int i;

  // Default values for initial array of spheres
  template.os[0] = zero;
  template.os[1] = zero;
  template.os[2] = zero;
  template.nm[0] = zero;
  template.nm[1] = zero;
  template.nm[2] = zero;
  template.thickness = zero;
  template.sph_d = one;

  // Copy default values to the array
  for(i=0; i<nsl; i++){
    sli[i] = template;
  }
}

void memory_clean_channels(CHA *cha, int nch)
{
  CHA template;
  int i;

  // Default values for initial array of spheres
  template.os[0] = zero;
  template.os[1] = zero;
  template.os[2] = zero;
  template.nm[0] = zero;
  template.nm[1] = zero;
  template.nm[2] = zero;
  template.radius = zero;
  template.sph_d = one;

  // Copy default values to the array
  for(i=0; i<nch; i++){
    cha[i] = template;
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

void display_stats(int vd, int bd, int is, int fs, int ns)
{
    fprintf(stdout, " Dimers valid      (if any): %4d %6.2lf %%\n",
            vd, (1e2*vd)/(2e0*ns) );
    fprintf(stdout, " Dimers broken     (if any): %4d %6.2lf %%\n",
            bd, (1e2*bd)/(2e0*ns) );
    fprintf(stdout, " Dimers spheres    (if any): %4d %6.2lf %%\n",
            2 * (vd), (2e2*(vd))/(1e0*ns) );
    fprintf(stdout, " Inclusion spheres (if any): %4d %6.2lf %%\n",
            is, (1e2*(is))/(1e0*ns));
    fprintf(stdout, " Free spheres      (if any): %4d %6.2lf %%\n",
            fs, (1e2*fs)/(1e0*ns) );
}


/* vim: set tw=80 ts=2 sw=2 et: */
