/*
 * utils.c
 *
 * File defines utility functions used in the setup and control of the program
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "config.h"

#ifdef PRNG_64BIT_MT19937
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

extern const int TYPE_INVALID;
extern const int TYPE_SPHERE;
extern const int TYPE_SPHERE_DIMER;
extern const int TYPE_DIMER;
extern const int TYPE_INCLUSION_SPHERE;

const char *fcc = "fcc";

int str_validate(const char *src, const char *tgt){
  int len_src = strlen(src);
  int len_tgt = strlen(tgt);

  if(len_src == len_tgt && !strncmp(src, tgt, len_tgt)){
    return 0;
  }
  return 1;
}

int count_particles_by_type(MODEL *md, SPH *sph, DIM3D *dim){
  int i;
  int tp_sph = 0, tp_inc_sph = 0;
  int tp_sph_dim = 0, tp_dim = 0;

  const char *fmt_data_corruption =
    " [%s] ERR: data corruption found in sphere and dimer arrays\n";
  for(i=0; i<md->Nsph; i++){
    if(sph[i].type == TYPE_SPHERE) tp_sph++;
    if(sph[i].type == TYPE_SPHERE_DIMER) tp_sph_dim++;
    if(sph[i].type == TYPE_INCLUSION_SPHERE) tp_inc_sph++;
  }

  for(i=0; i<md->Ndim; i++){
    if(dim[i].type == TYPE_DIMER) tp_dim++;
  }

  // Check for data corruption between tables;
  if(tp_sph_dim != 2*tp_dim){
    fprintf(stderr, fmt_data_corruption, __func__);
    return EXIT_FAILURE;
  }

  md->mtrx_sph = tp_sph;
  md->mtrx_dim = tp_dim;
  md->incl_sph = tp_inc_sph;
  return EXIT_SUCCESS;
}


int container_dimensions(MODEL *md, CONFIG cf){
  if(!str_validate(cf.symmetry, fcc)){
    for(int i=0; i<3; i++){
      md->box[i] = cf.cells[i] * sqrt(two);
    }
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

int number_of_spheres(CONFIG cf){
  // Variant for f.c.c. structure
  if(!str_validate(cf.symmetry, fcc)){
    return 4 * cf.cells[0] * cf.cells[1] * cf.cells[2];
  }
  return EXIT_FAILURE;
}

void bouble_sort_double(double *array, int s, int ascending){
  int sorted=1;
  int i,j;
  double tmp;

  if( ascending == 0){
    // Sort in descending order
    // Before starting, check if the array is sorted.
    for(i=0; i<s-1; i++){
      if(array[i] < array[i+1]){
        sorted = 0;
        break;
      }
    }

    // Start sorting if required
    if(sorted != 1){
      for(j=0;j<s;j++){
        for(i=0;i<s-1;i++){
          if(array[i] < array[i+1]){
            tmp = array[i+1];
            array[i+1] = array[i];
            array[i] = tmp;
          }
        }
      }
    }
  }else{
    // Sort in ascending order
    // Before starting, check if the array is sorted.
    for(i=1; i<s; i++){
      if(array[i-1] > array[i]){
        sorted = 0;
        break;
      }
    }

    // Start sorting if required
    if(sorted != 1){
      for(j=0;j<s;j++){
        for(i=1;i<s;i++){
          if(array[i-1] > array[i]){
            tmp = array[i-1];
            array[i-1] = array[i];
            array[i] = tmp;
          }
        }
      }
    }
  }

}

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
      // Draw another index if type is not correct.
      sp_id = (int) (u_RNG() * ns);
    }
  }
}

/*
 * find_critical_sphere(sph,ns)
 * Find an id of sphere with the lowest count of TYPE_SPHERE neigbours
 */

int find_critical_sphere(SPH *sph, int ns)
{
  int i, ngb_count;
  int min_ngb_count = 100000;
  int crit_fs = -1;

  for(i=0; i<ns; i++){
    if(sph[i].type == TYPE_SPHERE){
      // Count TYPE_SPHERE neighbours for i'th sphere
      ngb_count = count_typeX_sp_neighbours(sph, TYPE_SPHERE, i, ns);

      // remember 'i' index if calculated ngb_count is minimal
      if(ngb_count > 0 && ngb_count < min_ngb_count){
        min_ngb_count = ngb_count;
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
int count_typeX_sp_neighbours(SPH *sph, int x, int id, int nsph)
{
  int ic = 0, i;

  // Check if id is valid array index
  if(id < 0 && id >= nsph){
    return -1;
  }

  // count possible neighbours to connect.
  for(i=0; i<12; i++){
    if(sph[sph[id].ngb[i]].type == x){
      ic++;
    }
  }
  return ic;
}

int check_dimer_direction(DIM3D *dim, int i)
{
  int j;
  double fcc_dir[6][3];
  const char *fmt_unknown_orientation =
    " [%s] ERR: unrecognized dimer orientation\n";
  const char *fmt_error_details =
    " dimer ID %d O: %.16le %.16le %.16le (no match)\n";

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

  if(dim[i].type == TYPE_DIMER){
    for(j=0; j<6; j++){
      if( fabs(fabs(vdotu(fcc_dir[j], dim[i].O)) -sqrt(two)) < 1e-10 ){
        return j;
      }
    }
  }
  fprintf(stderr, fmt_unknown_orientation, __func__);
  fprintf(stderr, fmt_error_details, i, dim[i].O[0], dim[i].O[1], dim[i].O[2]);
  return -1;
}

/*
 * find_ngb_spheres(sph,ns,box)
 * Build neighburs data for spheres.
 */
int find_ngb_spheres(SPH sph[], int ns, double box[3]){
  //NOTE: this algorithm is not efficient
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
  return EXIT_SUCCESS;
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
    // NOTE update this comment when the above changes.
    if(dim[i].type == TYPE_INVALID){
      // Get atom indexes
      atom0 = dim[i].sph_ind[0];
      atom1 = dim[i].sph_ind[1];

      // Delete indexes of spheres forming a dimer
      dim[i].sph_ind[0] = -1;
      dim[i].sph_ind[1] = -1;

      // NOTE: is this still relevant ? ##########################
      // Set the free (non-channel) sphere (if any) to type '3'
      if(atom0 >= 0 && sph[atom0].type == 1){
        sph[atom0].type = TYPE_SPHERE;
      }
      if(atom1 >= 0 && sph[atom1].type == 1){
        sph[atom1].type = TYPE_SPHERE;
      }
      // ##########################################################

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
int update_dimer_parameters(MODEL md, DIM3D *dim, SPH *sph, int d){
  int atom0 = dim[d].sph_ind[0];
  int atom1 = dim[d].sph_ind[1];
  double O[3] = {zero,zero,zero};
  double R[3] = {zero,zero,zero};
  const char *fmt_data_corruption =
    " [%s] ERR: for dimer ID %d passed, invalid sphere indices found %d %d\n";

  if((atom0 > -1 && atom0 < md.Nsph) && (atom1 > -1 && atom1 < md.Nsph)){
    // Calculate current dimer orientation
    vector(sph[atom1].r, sph[atom0].r, O, md.box);

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

    return EXIT_SUCCESS;
  }else{
    fprintf(stderr, fmt_data_corruption, __func__, d, atom0, atom1);
    return EXIT_FAILURE;
  }
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
  template.type = TYPE_INVALID;
  template.dim_ind = TYPE_INVALID;
  template.r[0] = zero;
  template.r[1] = zero;
  template.r[2] = zero;
  template.lattice_ind[0] = TYPE_INVALID;
  template.lattice_ind[1] = TYPE_INVALID;
  template.lattice_ind[2] = TYPE_INVALID;
  template.d = one;

  for(i=0; i<12; i++){
    template.ngb[i] = TYPE_INVALID;
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
  template.type = TYPE_INVALID;
  template.sph_ind[0] = TYPE_INVALID;
  template.sph_ind[1] = TYPE_INVALID;
  template.R[0] = zero;
  template.R[1] = zero;
  template.R[2] = zero;
  template.O[0] = zero;
  template.O[1] = zero;
  template.O[2] = zero;
  template.L = one;

  for(i=0; i<22; i++){
    template.ngb[i][0] = TYPE_INVALID;
    template.ngb[i][1] = TYPE_INVALID;
  }

  // Copy default values to the array
  for(i=0; i<nd; i++){
    dim[i] = template;
  }

}

void memory_clean_inclusion(INC *inc, int n)
{
  INC template;
  int i;

  // Default values for an initial (undefined) inclusion
  template.os[0] = zero;
  template.os[1] = zero;
  template.os[2] = zero;
  template.nm[0] = zero;
  template.nm[1] = zero;
  template.nm[2] = zero;
  template.radius = zero;
  template.thickness = zero;
  template.sph_d = one;

  // Copy default values to the array
  for(i=0; i<n; i++){
    inc[i] = template;
  }
}

/*
 * u_RNG()
 *
 * Returns a pseudo random number uniformly distributed [0:1]
 */
double u_RNG()
{
#ifdef PRNG_64BIT_MT19937
  return genrand64_real1();
#endif

#ifdef PRNG_32BIT_MT19937
  return genrand_real1();
#endif

#ifdef PRNG_DRAND48
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
  #ifdef PRNG_64BIT_MT19937
    init_genrand64( s );
  #endif

  #ifdef PRNG_32BIT_MT19937
    init_genrand( s );
  #endif

  #ifdef PRNG_DRAND48
    srand48( s );
  #endif
}

/* vim: set tw=80 ts=2 sw=2 et: */
