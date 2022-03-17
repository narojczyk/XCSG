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

extern const char *fcc;

static int find_free_ngb_slot(SPH *sph);
static void bouble_sort_double(double *array, int s, int ascending);

/* # SEC ############## GENERAL STRUCTURE ################################### */

/*
 * container_dimensions()
 * Calculate box dimensions for the specified number of unit cells of
 * given symmetry.
 */
int container_dimensions(MODEL *md, CONFIG cf){
  if(!str_validate(cf.symmetry, fcc)){
    for(int i=0; i<3; i++){
      md->box[i] = cf.cells[i] * sqrt(two);
    }
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}

/*
 * number_of_spheres()
 * Returns the number of spheres for the specified number of unit cells of
 * given symmetry.
 */
int number_of_spheres(CONFIG cf){
  // Variant for f.c.c. structure
  if(!str_validate(cf.symmetry, fcc)){
    return 4 * cf.cells[0] * cf.cells[1] * cf.cells[2];
  }
  return EXIT_FAILURE;
}

/* # SEC ############## PARTICLES - SPHERES ################################# */

/*
 * count_typeX_sp_neighbours(sph, type_x, id)
 *
 * check the neighbours list of sphere 'id' and count all neighbours of type_x
 */
int count_typeX_sp_neighbours(SPH *sph, int type_x, int id, int ns)
{
  extern const int ngb_list_size;
  int ic = 0, i, i_ngb;

  // Check if id is valid array index
  if(id < 0 && id >= ns){
    return -1;
  }

  // count possible neighbours to connect.
  for(i=0; i<ngb_list_size; i++){
    i_ngb = sph[id].ngb[i];
    if(sph[i_ngb].type == type_x && sph[i_ngb].type_locked == 0){
      ic++;
    }
  }
  return ic;
}

/*
 * draw_sphere_typeX(sph,x,ns)
 *
 * Randomly select a sphere of type x
 */
int draw_sphere_typeX(SPH *sph, int x, int ns)
{
  extern const char *fmt_watchdog_activated;
  int sp_id = (int) (u_RNG() * ns);
  int nseek = 0, nmax = 1000 * ns;
  const char *fmt_no_valid_sphere = " [%s] ERR: Cannot locate type %d sphere\n";

  while(1){
    // Check not to go outside the sphere tab
    sp_id = (sp_id < ns) ? sp_id : ns - 1;
    // Check if the type is ok
    if(sph[sp_id].type == x && sph[sp_id].type_locked == 0){
      return sp_id;
    }else{
      // Draw another index if type is not correct.
      sp_id = (int) (u_RNG() * ns);
    }
    // Watchdog: count attempts to select valid neighbour
    nseek++;
    if(nseek > nmax){
      fprintf(stderr, fmt_no_valid_sphere, __func__, x);
      fprintf(stderr, fmt_watchdog_activated, __func__);
      return -1;
    }
  }
}

/*
 * draw_ngb_sphere_typeX()
 * Randomly select neighbour of sphere 'sph_ind'. The neighbour must be of
 * type 'x'.
 */
int draw_ngb_sphere_typeX(SPH *sph, int x, int sph_ind)
{
  extern const int TYPE_SPHERE_DIMER;
  extern const int ngb_list_size;
  extern const char *fmt_watchdog_activated;
  const char *fmt_no_valid_neighbour =
    " [%s] ERR: Cannot locate type %d neighbour of %d sphere\n";

  int nseek = 0, nmax = 500;
  int sngb_id, sngb_type, sngb_type_lock, rand_ngb, valid_ngb;

  // Check if sph_ind is a valid array index
  if(sph_ind == -1){
    return -1;
  }

  // Select random neighbour of sph_ind
  do{
    valid_ngb = 0;
    // Select neighbour index randomly
    rand_ngb = get_random_range12();

    // Get neighbour id and type
    sngb_id = sph[sph_ind].ngb[rand_ngb];
    sngb_type = sph[ sngb_id ].type;
    sngb_type_lock = sph[ sngb_id ].type_locked;

    // Select type-x sphere
    if(sngb_type == x && sngb_type_lock == 0) {
      valid_ngb = 1;
    }
    // Security check not to select type-2 spheres EVER!
    if(sngb_type == TYPE_SPHERE_DIMER){
      valid_ngb = 0;
    }
    // Watchdog: count attempts to select valid neighbour
    nseek++;
    if(nseek > nmax){
      fprintf(stderr, fmt_no_valid_neighbour, __func__, x, sph_ind);
      fprintf(stderr, fmt_watchdog_activated, __func__);
      return -1;
    }
  }while(!valid_ngb);

  return sngb_id;
}

/*
 * find_critical_sphere(sph, type_x, ns)
 * Find an id of sphere of type_x, with the lowest count of type_x neigbours
 */

int find_critical_sphere(SPH *sph, int type_x, int ns)
{
  int i, ngb_count;
  int min_ngb_count = 12;
  int crit_fs = -1;

  for(i=0; i<ns; i++){
    if(sph[i].type == type_x && sph[i].type_locked == 0){
      // Count type_x neighbours for i'th sphere
      ngb_count = count_typeX_sp_neighbours(sph, type_x, i, ns);

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
 * find_ngb_spheres(sph,ns,box)
 * Build neighburs data for spheres.
 */
int find_ngb_spheres(SPH sph[], int ns, double box[3]){
  const char *fmt_too_small_ngb_list =
    " [%s] ERR: neighbour list for %d is nof full\n";
  const char *fmt_no_free_ngb_slot =
    " [%s] ERR: no free slots for neighbour id of sph %d\n";
  int i, j, slot_id, ngb_slot_id;
  const double lim = 5e-11;

  // Loop over all spheres
  for(i=0; i<ns-1; i++){
    // get free neighbour index
    slot_id = find_free_ngb_slot(&sph[i]);

    // Loop over all spheres other than 'i'
    for(j=i+1; j<ns; j++){
      if(fabs(distance(sph[i].r,sph[j].r, box) - one) < lim ){
        // Save ngb index for i'th sphere
        sph[i].ngb[slot_id++] = j;

        // Find first free slot in the ngb tab of j'th sphere
        ngb_slot_id = find_free_ngb_slot(&sph[j]);
        if(ngb_slot_id == TYPE_INVALID){
          fprintf(stderr, fmt_no_free_ngb_slot, __func__, j);
          return EXIT_FAILURE;
        }

        // And save the index of i'th sphere as neighbour of 'j'
        sph[j].ngb[ngb_slot_id] = i;
      }
    } // j-loop end

    // Sanity check
    if(find_free_ngb_slot(&sph[i]) != TYPE_INVALID){
      fprintf(stderr, fmt_too_small_ngb_list, __func__, i);
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

static int find_free_ngb_slot(SPH *sph){
  extern const int ngb_list_size;

  for(int i=0; i<ngb_list_size; i++){
    if(sph->ngb[i] == TYPE_INVALID){
      return i;
    }
  }

  return TYPE_INVALID;
}

/*
 * sph_assign_lattice_indexes( sph, ns)
 *
 * Based on posisions at closepacking assign lattice indexes (x,y,z) to all
 * spheres
 * TODO: Maybe move to utils or misc section
 */

int sph_assign_lattice_indexes( SPH *sph, int ns)
{
  const char *fmt_sph_number = " Sphere number %d\n";
  const char *fmt_unassigned_latt_idx =
    " [%s] ERR: Unasigned lattice index for the direction %d\n";
  const char *fmt_index_out_of_range = " [%s] ERR: index out of range\n";
  int tmp_size = 300;
  int dir, i=0, j=0, c, present;
  double uniq_coords[tmp_size];
  double c_value;

  for(dir=0; dir<3; dir++){
    c=0;
    present=0;

    // Assign the coordinate of the first sphere into the array
    uniq_coords[c++] = sph[0].r[dir];

    // Counting from the next sphere, check whether its coordinate
    // matches any of previously found and stored in the array
    for(i=1; i<ns; i++){
      // flag to signal that value was previously found
      present=0;
      // get the coordinate of the next sphere
      c_value = sph[i].r[dir];
      // sweep the array from the start to the current number of found
      // elements and try to match the current coordinate
      for(j=0; j<c; j++){
        if( fabs(c_value/uniq_coords[j] -1) < 1e-10 ){
          // set flag and brake the loop if value already found in the array
          present=1;
          break;
        }
      }
      // if the flag is not set after the sweep add the current coordinate
      // to the array and increment the array counter
      if(present==0){
        if(c < tmp_size){
          uniq_coords[c++] = c_value;
        }else{
          fprintf(stderr, fmt_index_out_of_range, __func__);
          return EXIT_FAILURE;
        }
      }
    }

    // Make sure the array is sorted ascending
    bouble_sort_double(uniq_coords, c, 1);

    // For all spheres compare their coordinats ...
    for(i=0; i<ns; i++){
      // ... with consecutive values from the array ...
      for(j=0; j<c; j++){
        if( fabs(sph[i].r[dir]/uniq_coords[j] -1) < 1e-10 ){
          // ... and assing respective index value when found
          sph[i].lattice_ind[dir] = j+1;
          break;
        }
      }
      // For security check if the index does not remain unassigned
      if(sph[i].lattice_ind[dir] == -1){
        fprintf(stderr, fmt_unassigned_latt_idx, __func__, dir);
        fprintf(stderr, fmt_sph_number, i);
        return EXIT_FAILURE;
      }
    }
  }
  return EXIT_SUCCESS;
}

static void bouble_sort_double(double *array, int s, int ascending){
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

/* # SEC ############## PARTICLES - DIMERS ################################## */

/*
 * check_dimer_direction()
 * Return id 'j' of one of six possible dimer orientations in f.c.c. symmetry
 */
int check_dimer_direction(DIM3D *dim, int i)
{
  extern const int TYPE_DIMER;
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

/* # SEC ############## MISCELLANEOUS ####################################### */

/*
 * count_particles_by_type()
 * Calculate all the different particle types currently in the structure
 */
int count_particles_by_type(MODEL *md, SPH *sph, DIM3D *dim){
  extern const int TYPE_SPHERE;
  extern const int TYPE_SPHERE_DIMER;
  extern const int TYPE_DIMER;
  extern const int TYPE_INCLUSION_SPHERE;
  extern const int TYPE_INCLUSION_SPHERE_DIMER;
  extern const int TYPE_INCLUSION_DIMER;
  int i;
  int tp_sph = 0, tp_inc_sph = 0, tp_inc_sph_dim = 0;
  int tp_sph_dim = 0, tp_dim = 0, tp_inc_dim = 0;
  const char *fmt_data_corruption =
    " [%s] ERR: data corruption found in sphere and dimer arrays\n";

  for(i=0; i<md->Nsph; i++){
    if(sph[i].type == TYPE_SPHERE) tp_sph++;
    if(sph[i].type == TYPE_SPHERE_DIMER) tp_sph_dim++;
    if(sph[i].type == TYPE_INCLUSION_SPHERE) tp_inc_sph++;
    if(sph[i].type == TYPE_INCLUSION_SPHERE_DIMER) tp_inc_sph_dim++;
  }

  for(i=0; i<md->Ndim; i++){
    if(dim[i].type == TYPE_DIMER) tp_dim++;
    if(dim[i].type == TYPE_INCLUSION_DIMER) tp_inc_dim++;
  }

  // Check for data corruption between tables;
  if(tp_sph_dim != 2*tp_dim || tp_inc_sph_dim != 2*tp_inc_dim){
    fprintf(stderr, fmt_data_corruption, __func__);
    return EXIT_FAILURE;
  }

  md->mtrx_sph = tp_sph;
  md->mtrx_dim = tp_dim;
  md->incl_sph = tp_inc_sph;
  md->incl_dim = tp_inc_dim;
  return EXIT_SUCCESS;
}

/*
 * str_validate(s, t)
 * Compare two strings s and t, and check if they match
 */
int str_validate(const char *src, const char *tgt){
  int len_src = strlen(src);
  int len_tgt = strlen(tgt);

  if(len_src == len_tgt && !strncmp(src, tgt, len_tgt)){
    return 0;
  }
  return 1;
}


/* # SEC ############## P. R. N. G. ######################################### */

/*
 * get_random_range12()
 *
 * Draw random integer number in the range <0:11>
 */
int get_random_range12(){
  int rn = (int) (u_RNG() * 12);
  return (rn < 12) ? rn : 11;
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

/* # SEC ############## MEMORY CLEANING AND INITIATING ###################### */

/*
 * memory_clean_spheres(sph, ns)
 * memory_clean_dimers(dim, nd)
 * memory_clean_inclusion(inc, n)
 *
 * Assign default values to data structures' arrays
 */
void memory_clean_spheres(SPH *sph, int ns)
{
  extern const int ngb_list_size;
  SPH template;
  int i;

  // Default values for initial array of spheres
  template.type = TYPE_INVALID;
  template.type_locked = 0;
  template.dim_ind = TYPE_INVALID;
  template.r[0] = zero;
  template.r[1] = zero;
  template.r[2] = zero;
  template.lattice_ind[0] = TYPE_INVALID;
  template.lattice_ind[1] = TYPE_INVALID;
  template.lattice_ind[2] = TYPE_INVALID;
  template.d = one;

  for(i=0; i<ngb_list_size; i++){
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

/* vim: set tw=80 ts=2 sw=2 et: */
