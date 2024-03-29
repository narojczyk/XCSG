
/*
 * structure.c
 *
 * File defines functions to preform operations on the DC structure
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "data.h"
#include "utils.h"
#include "algebra.h"
#include "structure.h"

extern const double zero;
extern const double one;
extern const double two;

extern const int TYPE_MATRIX_BASE;
extern const int TYPE_INCLUSION_BASE;
extern const int TYPE_INVALID;
extern const int TYPE_SPHERE;
extern const int TYPE_SPHERE_DIMER;
extern const int TYPE_DIMER;
extern const int TYPE_INCLUSION_SPHERE;
extern const int TYPE_INCLUSION_SPHERE_DIMER;

extern const char *fmt_internal_call_failed;

static int set_fcc(SPH *sph, int ns, int cells[3]);
static int chk_dim_config(PARTICLES pts, double box[3], int d0, int d1);
static int validate_distrib(int od[6], int nd1, int step);
static int zipper(MODEL md, PARTICLES pts, int s_id, int workload, int TYPE);

static double shortest_distance_to_line(MODEL md, POINT r, POINT a, POINT b);
static double p_AB_distance(POINT p, POINT a, POINT b);

static void display_dimer_distribution(int od[6]); //(->i/o ?)
static void flip_dimers(MODEL md, PARTICLES pts, int od[6], int d_ids[2]);
static void find_valid_cluster(MODEL md, PARTICLES pts, int vclust[2]);
static void make_dimer(PARTICLES pts, MODEL md, int s1, int s2, int type_tgt);
static void update_dimer_type(DIM3D *dim, SPH *sph, int type, double sph_d);

/* # SEC ############## SYMMETRY OF THE STRUCTURE ########################### */

/*
 * set_structure()
 * Call particular function based on the selected symmetry.
 */
int set_structure(CONFIG cf, MODEL md, SPH *sph){
  extern const char *fcc;
  const char *fmt_generate_structure = " Generating initial %s structure\n";

  if(!str_validate(cf.symmetry, fcc)){
    fprintf(stdout, fmt_generate_structure, fcc);
    return set_fcc(sph, md.Nsph, cf.cells);
  }

  return EXIT_FAILURE;
}

/*
 * set_fcc(sph, ns, cells)
 * Set fcc structure of ns spheres in a cubic system of fcc cells at the
 * edge. The edge length is assumed sqrt(2)
 */
static int set_fcc(SPH *sph, int ns, int cells[3]){
  int i,j,k,l, x=0, y=0, z=0, last_index;
  // Lattice constant
  double a = sqrt(two);
  double half_a = a/two;
  double sys_half_size[3] = {cells[0]*half_a, cells[1]*half_a, cells[2]*half_a};
  double com[3]={zero,zero,zero};
  const char *fmt_index_out_of_range = " [%s] ERR: index out of range\n";
  const char *fmt_index_value = " index %d > MODEL.Nsph\n";

  for(i=0, j=i+1, k=i+2, l=i+3; i<ns; i+=4, j+=4, k+=4, l+=4){
    sph[i].r[0] = -sys_half_size[0] + x * a;
    sph[i].r[1] = -sys_half_size[1] + y * a;
    sph[i].r[2] = -sys_half_size[2] + z * a;

    sph[j].r[0] = -sys_half_size[0] + x * a + half_a;
    sph[j].r[1] = -sys_half_size[1] + y * a;
    sph[j].r[2] = -sys_half_size[2] + z * a + half_a;

    sph[k].r[0] = -sys_half_size[0] + x * a;
    sph[k].r[1] = -sys_half_size[1] + y * a + half_a;
    sph[k].r[2] = -sys_half_size[2] + z * a + half_a;

    sph[l].r[0] = -sys_half_size[0] + x * a + half_a;
    sph[l].r[1] = -sys_half_size[1] + y * a + half_a;
    sph[l].r[2] = -sys_half_size[2] + z * a;

    sph[i].d = one;
    sph[j].d = one;
    sph[k].d = one;
    sph[l].d = one;

    // set default type as free spheres
    sph[i].type = TYPE_SPHERE;
    sph[j].type = TYPE_SPHERE;
    sph[k].type = TYPE_SPHERE;
    sph[l].type = TYPE_SPHERE;

    // clear dimer indexes
    sph[i].dim_ind = TYPE_INVALID;
    sph[j].dim_ind = TYPE_INVALID;
    sph[k].dim_ind = TYPE_INVALID;
    sph[l].dim_ind = TYPE_INVALID;

    // Increment cell in x direction
    x++;
    // Increment cell in y direction
    if( x == cells[0]){
      x = 0;
      y++;
    };
    // increment cell in z direction
    if( y == cells[1]){
      y = 0;
      z++;
    };
  }
  // Security check
  last_index = (x+1)*(y+1)*(z+1)*4;
  if( last_index > ns ){
    fprintf(stderr, fmt_index_out_of_range, __func__);
    fprintf(stderr, fmt_index_value, last_index);
    return EXIT_FAILURE;
  }
  // Calculate the center of mass of the new structure
  for(i=0; i<ns; i++){
    com[0] += sph[i].r[0];
    com[1] += sph[i].r[1];
    com[2] += sph[i].r[2];
  }
  com[0] /= ns;
  com[1] /= ns;
  com[2] /= ns;

  // Move the center of mass of the new structure to point 0
  for(i=0; i<ns; i++){
    sph[i].r[0] -= com[0];
    sph[i].r[1] -= com[1];
    sph[i].r[2] -= com[2];
  }

  return EXIT_SUCCESS;
}

/*
 * find_smallest_coordinates()
 *
 * Look for sphere with position closest to the corner of the periodic box
 * with the lowest coordinates
 */
void find_smallest_coordinates(SPH *sph, MODEL *md){
  int i;
  double lcs[3] = {zero, zero, zero};
  double ref[3] = {-md->box[0]/two, -md->box[1]/two, -md->box[2]/two};
  double dist, min_dist = md->box[0]; // any relatively large number will do

  for(i=0; i < md->Nsph; i++){
    // Calculate distance (no p.b.) of a sphere from a reference point
    dist = distance_absolute(ref, sph[i].r);
    if(dist < min_dist){
      // remember coordinates if they are closer to ref than the previous one
      lcs[0] = sph[i].r[0];
      lcs[1] = sph[i].r[1];
      lcs[2] = sph[i].r[2];
      // remember the new minimal distance
      min_dist = dist;
    }
  }

  // Store the coordinates in the MODEL structure
  for(i=0; i<3; i++){
    md->lcs[i] = lcs[i];
  }
}

/* # SEC ############## INCLUSIONS ########################################## */

/*
 * make_cluster()
 * Calculates the distance of the sphere 'i' from the centre of inclusion
 * cluster.
 *
 */
void make_cluster(MODEL md, INC *inc, PARTICLES pts){
  SPH *sph = pts.spheres, *s = NULL;
  DIM3D *dim = pts.dimers;
  int i,k;
  int type_lock = (inc->tgt_Nmer == 1) ? 1 : 0;
  double centre[3] = {
    md.lcs[0] + inc->os[0], md.lcs[1] + inc->os[1], md.lcs[2] + inc->os[2]};
  double radius = inc->radius, sph_diameter = inc->sph_d;

  // Test all sphere positions against channel's axis through 'a' and 'b'
  for(i = 0; i < md.Nsph; i++){
    s = &sph[i];
    // Skip spheres that already belong to any inclusion
    if((*s).type < TYPE_INCLUSION_BASE){
      // Check if the distance from the sphere position 'r' to
      // the cluster centre is less than inclusion radius
      if(distance((*s).r, centre, md.box) < radius){
        // Based on the sphere type
        if((*s).type == TYPE_SPHERE){
          // ... mark regular spheres as 'inclusion-spheres'
          (*s).type = TYPE_INCLUSION_SPHERE;
          // ... and assign the respective inclusion-sphere-diameter
          (*s).d = sph_diameter;
          // set type lock flag for this inclusion
          (*s).type_locked = type_lock;
        }else if ((*s).type == TYPE_SPHERE_DIMER){
          // ... mark regular dimers as inclusion dimers only if their
          // centres of mass are inside the channel
          k = (*s).dim_ind;
          if(distance(dim[k].R, centre, md.box) < radius){
            update_dimer_type(&dim[k], sph, TYPE_INCLUSION_SPHERE_DIMER,
                              inc->sph_d);
          }
        }
      }
    } // end if sph.type <TYPE_INCLUSION_BASE condition
  }
}

/*
 * make_channel()
 * Calculates the distance of the sphere 'i' from the line given by vector 'c',
 * with respect to two points on the line: lower left corner and center point
 * of the channel axis (this is required because periodic boundaries may select
 * periodic image that is not the nearest to the channel axis)
 *
 */
void make_channel(MODEL md, INC *inc, PARTICLES pts){
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;
  POINT a, b, r;
  int i,k;
  int type_lock = (inc->tgt_Nmer == 1) ? 1 : 0;
  double n[3] = {inc->nm[0], inc->nm[1], inc->nm[2]};
  double L = 0.85 * md.box[0]; // an arbitrary number (box/2 < L < box)

  // Determine two points on an axis parallel to inc->nm
  for(i=0; i<3; i++){
    // (i) coordinates of the sphere closest to ref (the vertex of the cube
    // with the lowest coordinates), off setted by specified inc->os
    a.c[i] = md.lcs[i] + inc->os[i];
    // (ii) find the coordinates of the channel's second point
    // by translating 'a' by 'L' in the direction of normal vector
    b.c[i] = a.c[i] + L * n[i];
  }

  // Test all sphere positions against channel's axis through 'a' and 'b'
  for(i=0; i<md.Nsph; i++){
    // Skip spheres that already belong to any inclusion
    if(sph[i].type < TYPE_INCLUSION_BASE){
      // Check if the shortest distance from the sphere position 'r' to
      // the channel's axis is less than inclusion radius
      r.c[0] = sph[i].r[0];
      r.c[1] = sph[i].r[1];
      r.c[2] = sph[i].r[2];
      if(shortest_distance_to_line(md, r, a, b) < inc->radius){
        // Based on the sphere type
        if(sph[i].type == TYPE_SPHERE){
          // ... mark regular spheres as 'inclusion-spheres'
          sph[i].type = TYPE_INCLUSION_SPHERE;
          // ... and assign the respective inclusion-sphere-diameter
          sph[i].d = inc->sph_d;
          // set type lock flag for this inclusion
          sph[i].type_locked = type_lock;
        }else if (sph[i].type == TYPE_SPHERE_DIMER){
          // ... mark regular dimers as inclusion dimers only if their
          // centres of mass are inside the channel
          k = sph[i].dim_ind;
          r.c[0] = dim[k].R[0];
          r.c[1] = dim[k].R[1];
          r.c[2] = dim[k].R[2];
          if(shortest_distance_to_line(md, r, a, b) < inc->radius){
            update_dimer_type(&dim[k], sph, TYPE_INCLUSION_SPHERE_DIMER,
                              inc->sph_d);
          }
        }
      }
    } // end if sph.type <TYPE_INCLUSION_BASE condition
  } // end for loop
}

/*
 * make_slit(dim,sph,c,ns)
 * Introduces a plane (single atomic layer) into dimer structure. The plane
 * is described by coordinate axis origin and normal vector 'c'. Dimers who's
 * atoms belong to the plane, are flagged for braking.
 * NOTE: periodic boundaries are not taken into account here
 */
void make_slit(MODEL md, INC *inc, PARTICLES pts){
  SPH *sph = pts.spheres, *s = NULL;
  DIM3D *dim = pts.dimers, *d = NULL;
  POINT p[7];
  int i, j;
  int type_lock = (inc->tgt_Nmer == 1) ? 1 : 0;
  double n[3] = {inc->nm[0], inc->nm[1], inc->nm[2]};
  double dist;

  // Transform the plane vector to unit vector;
  vnorm(n);

  // Loop over all spheres in the structure
  for(i=0; i<md.Nsph; i++){
    s = &sph[i];

    // Skip spheres that already belong to any inclusion
    if((*s).type < TYPE_INCLUSION_BASE){
      // Get the i'th sphere position relative to the plane inside periodic md.box
      p[0].c[0] = (*s).r[0] - inc->os[0];
      p[0].c[1] = (*s).r[1] - inc->os[1];
      p[0].c[2] = (*s).r[2] - inc->os[2];

      // ... and in the six of the surrounding images
      for(j=1; j<7; j++){
        p[j] = p[0];
      }
      p[1].c[0] += md.box[0];
      p[2].c[1] += md.box[1];
      p[3].c[2] += md.box[2];
      p[4].c[0] -= md.box[0];
      p[5].c[1] -= md.box[1];
      p[6].c[2] -= md.box[2];

      for(j=0; j<=6; j++){
        // Calculate the dot product of the normal to the plane located at
        // axes origin and a position vector of i'th sphere
        dist = zero - n[0]*p[j].c[0] - n[1]*p[j].c[1] - n[2]*p[j].c[2];

        // To get the distance of a sphere's center from the plane one should
        // (in principle) divide the above by the length of the plane's normal.
        // This is not done here, as the normal is unit-normed.

        // If the sphere lies within the plane, include the sphere into
        // plane and continue to the next sphere
        if(fabs(dist) < inc->thickness){
          // Based on the sphere type
          if( (*s).type == TYPE_SPHERE ){
            // ... mark regular spheres as 'inclusion-spheres'
            (*s).type = TYPE_INCLUSION_SPHERE;
            // ... and assign the respective inclusion-sphere-diameter
            (*s).d = inc->sph_d;
            // set type lock flag for this inclusion
            (*s).type_locked = type_lock;
          }else if ( (*s).type == TYPE_SPHERE_DIMER ){
            // ... mark regular dimers as inclusion dimers only if their
            // centres of mass are inside the slit
            d = &dim[(*s).dim_ind];
            if(fabs(zero
                - n[0] * (*d).R[0] - n[1] * (*d).R[1] - n[2] * (*d).R[2])
                  < inc->thickness){
              update_dimer_type(d, sph, TYPE_INCLUSION_SPHERE_DIMER,
                                inc->sph_d);
            }
          }

          // Escape the j-loop if sphere is found to lie on the plane
          break;
        }
      }
    } // end if sph.type <TYPE_INCLUSION_BASE condition
  } // end for loop
}

/*
 * shortest_distance_to_line()
 *
 * Calculates the shortest distance from point r to the line described by
 * wersor c and two points on the axis lc and cp
 */
static double shortest_distance_to_line(MODEL md, POINT r, POINT a, POINT b){
  POINT p[6]; // periodic images of r
  int i;
  double min_distance, distance;

  // Get coordinates of 'r' in six adjacent periodic boxes
  for(i=0; i<6; i++){
    p[i] = r;
  }
  p[0].c[0] += md.box[0];
  p[1].c[1] += md.box[1];
  p[2].c[2] += md.box[2];
  p[3].c[0] -= md.box[0];
  p[4].c[1] -= md.box[1];
  p[5].c[2] -= md.box[2];

  // Get distance for r inside the box
  min_distance = p_AB_distance(r, a, b);

  // Get distance for r in respective images and store the shortest one
  for(i=0; i<6; i++){
    distance = p_AB_distance(p[i], a, b);
    if(distance < min_distance){
      min_distance = distance;
    }
  }

  return min_distance;
}

/*
 * p_AB_distance()
 *
 * Returns the shortest distance between point p and a line given by points
 * a and b, in 3D
 * Extention of the example for 2D case:
 * The area of the parallelogram spanned by points A,B (on the line), and P is
 * |(B-A)x(C-A)| = |(B.x−A.x)(P.y−A.y)−(B.y−A.y)(P.x−A.x)|. If we divide this by
 * the length Sqrt[(B-A)^2] = Sqrt[(B.x−A.x)^2 + (B.y−A.y)^2] of its base,
 * we obtain ist height. Final formula:
 * |(B.x−A.x)(P.y−A.y)−(B.y−A.y)(P.x−A.x) | / Sqrt[(B.x−A.x)^2 + (B.y−A.y)^2]
 */
static double p_AB_distance(POINT p, POINT a, POINT b){
  double ab_xSQ = (a.c[0] - b.c[0]) * (a.c[0] - b.c[0]);
  double ab_ySQ = (a.c[1] - b.c[1]) * (a.c[1] - b.c[1]);
  double ab_zSQ = (a.c[2] - b.c[2]) * (a.c[2] - b.c[2]);

  double exp_1 =  a.c[1]*b.c[0] - a.c[0]*b.c[1] - a.c[1]*p.c[0]
                  + b.c[1]*p.c[0] + a.c[0]*p.c[1] - b.c[0]*p.c[1];
  double exp_2 =  a.c[2]*b.c[0] - a.c[0]*b.c[2] - a.c[2]*p.c[0]
                  + b.c[2]*p.c[0] + a.c[0]*p.c[2] - b.c[0]*p.c[2];
  double exp_3 =  a.c[2]*b.c[1] - a.c[1]*b.c[2] - a.c[2]*p.c[1]
                  + b.c[2]*p.c[1] + a.c[1]*p.c[2] - b.c[1]*p.c[2];

  double denominator = exp_1 * exp_1 + exp_2 * exp_2 + exp_3 * exp_3;

  return one / sqrt( (ab_xSQ + ab_ySQ + ab_zSQ) / denominator );
}

/*
 * update_dimer_type()
 * Changes the dimer and its spheres to the selected 'type' and assign 'sph_d'
 * as sphere diameters
 */
static void update_dimer_type(DIM3D *dim, SPH *sph, int type, double sph_d){
  SPH *sph0 = &sph[dim->sph_ind[0]];
  SPH *sph1 = &sph[dim->sph_ind[1]];

  dim->type = type;
  sph0->type = type;
  sph1->type = type;

  sph0->d = sph_d;
  sph1->d = sph_d;
}

/* # SEC ############## INSERTING DIMERS #################################### */

/*
 * introduce_random_dimers()
 * Randomly (where possible) connect neighbouring spheres of type_src into
 * dimers. In critical cases start with spheres with fewest possible
 * connections.
 * NOTE: make the function return success or failure code.
 * BUG: infinite loop possible.
 */
void introduce_random_dimers(PARTICLES pts, MODEL *md, int type_src,
                             int type_tgt){
  SPH *sph = pts.spheres;
  int s_id=-1, s_ngb_id=-1, s_ngb_qty=0;
  int nsph = md->Nsph;
  const char *fmt_mk_dimers_at_random =
    "\n Randomly join type-%d spheres into dimers\n";

  fprintf(stdout, fmt_mk_dimers_at_random, type_src);

  do{
    // Find a type_src sphere with the lowes count of chanses to form a dimer
    s_id = find_critical_sphere(sph, type_src, nsph);
    s_ngb_qty = count_typeX_sp_neighbours(sph, type_src, s_id, nsph);

    // If the possibilities are high enough, select sphere randomly
    if(s_ngb_qty >= 5){
      s_id = draw_sphere_typeX(sph, type_src, nsph);
      s_ngb_qty = count_typeX_sp_neighbours(sph, type_src, s_id, nsph);
    }

    // Randomly select neighbour of s_id only when valid index found
    if(s_id >= 0 && s_id < nsph){
      s_ngb_id = draw_ngb_sphere_typeX(sph, type_src, s_id);
    }
    // Create a valid dimer from the pair of spheres
    if( s_id != -1 && s_ngb_id != -1){
      make_dimer(pts, *md, s_id, s_ngb_id, type_tgt);
      // Update free spheres count
      md->mtrx_sph -= 2;
    }
  }while(s_ngb_qty > 0 && md->mtrx_sph > 0);

  if(count_particles_by_type(md, pts) == EXIT_FAILURE){
    fprintf(stderr, fmt_internal_call_failed, __func__);
  }
}

/*
 * introduce_dimers_by_zipper()
 * Using zipper mechanism eliminate any remaining spheres. The TYPE variable
 * is used to restrict zipper to run on matrix (TYPE_MATRIX_BASE) or inclusion
 * (TYPE_INCLUSION_BASE) particles.
 */
int introduce_dimers_by_zipper(PARTICLES pts, MODEL *md, int TYPE){
  extern const int lazy;
  SPH *sph = pts.spheres;
  // local redefinition of TYPE* variables basen on whether matix or inclusion
  // particles are to be considered
  const int TYPE_sphere = TYPE + TYPE_SPHERE;

  int i, zipper_steps_passed, zipper_runs = 0, zip_init_sph;
  const char *fmt_ziper_run =
    " Zipper from sphere ID %5d completed after %7d steps\n";
  const char *fmt_required_zipper_runs =
    " Running zipper %d times to create dimers\n";
  const char *fmt_mk_dimers_with_zipper =
    "\n Connecting remaining spheres with zipper\n";
  const char *fmt_illegal_type =
    " [%s] ERR: Not recognized value of TYPE (%d)\n";
  const char *fmt_no_valid_startpoint =
    " No valid starting point to initiate zipper\n";

  if(TYPE != TYPE_MATRIX_BASE && TYPE != TYPE_INCLUSION_BASE){
    fprintf(stderr, fmt_mk_dimers_with_zipper, __func__);
    fprintf(stderr, fmt_illegal_type, __func__, TYPE);
    return EXIT_FAILURE;
  }

  // Calculate the number of zipper runs to perform
  // NOTE: This condition only holds while there are two base types
  // of particles: matrix and inclusion
  zipper_runs = (TYPE == TYPE_MATRIX_BASE)
    ? (md->mtrx_sph - (md->mtrx_sph & 1))/2
    : (md->incl_sph - (md->incl_sph & 1))/2;

  if(zipper_runs != 0){
    fprintf(stdout, fmt_mk_dimers_with_zipper);
    fprintf(stdout, fmt_required_zipper_runs, zipper_runs);
  }else{
    return EXIT_SUCCESS;
  }

  while(zipper_runs != 0){
    // Select a valid type-sphere object as the starting point for zipper
    zip_init_sph = -1;
    for(i=0; i<md->Nsph; i++){
      if(sph[i].type == TYPE_sphere && sph[i].type_locked == 0){
        zip_init_sph = i;
        break;
      }
    }

    // Do not start zipper when no valid staring point found
    if(zip_init_sph == -1){
      fprintf(stdout, fmt_no_valid_startpoint);
      break;
    }

    // Start zipper from selected sphere
    fprintf(stdout,"  (%d)", zipper_runs);
    zipper_steps_passed = zipper( *md, pts, zip_init_sph, lazy, TYPE);
    // Terminate if zipper failed
    if(zipper_steps_passed == EXIT_FAILURE){
      fprintf(stderr, fmt_internal_call_failed, __func__);
      return EXIT_FAILURE;
    }
    fprintf(stdout, fmt_ziper_run, zip_init_sph, zipper_steps_passed);
    zipper_runs--;
  }

  if(count_particles_by_type(md, pts) == EXIT_FAILURE){
    fprintf(stderr, fmt_internal_call_failed, __func__);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

/*
 * refine_dimer_distribution()
 * Reorganize molecules to get best possible DC distribution of dimers
 * NOTE: This procedure will most likely fail on sets with low amount of dimers
 * or groups of molecules with high special confinements. Thus, TYPE of
 * particles to run zipper on has been hardcoded as matrix type only.
 */
int refine_dimer_distribution(PARTICLES pts, MODEL md, int od[6]){
  extern const int diligent;
  extern const int run_zipper;
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;

  int zip_init_sph = -1, flip_count = 0, flipable_dimers[2] = {-1, -1};

  const char *fmt_ziper_run =
    " Zipper from sphere ID %5d completed after %7d steps\n";
  const char *fmt_dimer_distr_header = " %-28s %5s %5s %5s %5s %5s %5s\n";

  flip_count = 0;
  fprintf(stdout,fmt_dimer_distr_header, "Initial distribution:",
          "[110]", "[T10]", "[101]", "[T01]", "[011]", "[0T1]");
  if(!validate_distrib(od, md.mtrx_dim, flip_count)){
    do{
      // Find a valid dimer configuration to flip orientations
      find_valid_cluster(md, pts, flipable_dimers);
      // Flip dimers
      if(flipable_dimers[0] != -1 && flipable_dimers[1] != -1){
        flip_dimers(md, pts, od, flipable_dimers);
      }

      // Run zipper every once in a while (zipper length = (const) diligent)
      if(flip_count % run_zipper == 0){
        do{
          // Select random sphere (forming any dimer) to start from
          zip_init_sph = (int) (u_RNG() * md.Nsph);
          zip_init_sph = (zip_init_sph < md.Nsph ? zip_init_sph : md.Nsph - 1);
        }while(sph[zip_init_sph].type != TYPE_SPHERE_DIMER);

        fprintf(stdout, fmt_ziper_run, zip_init_sph,
                zipper(md, pts, zip_init_sph, diligent, TYPE_MATRIX_BASE));

        // Display distribution after zipper
        if(test_dimer_distribution(dim,  od, md.Ndim) == EXIT_FAILURE){
          fprintf(stderr, fmt_internal_call_failed, __func__);
          return EXIT_FAILURE;
        }
        fprintf(stdout," %-27s :", "Zipper distribution");
        display_dimer_distribution(od);
      }

      flip_count++;
    }while(validate_distrib(od, md.mtrx_dim, flip_count) == 0);
  }
  return EXIT_SUCCESS;
}

/* # SUB-SEC ########## INSERTING DIMERS - PARTICLE FUNCTIONS ############### */

/*
 * ziper(dim,sph,box,sph_ind)
 *
 * Run zipper on the structure, starting from sphere sph_ind, and continue
 * untill another free sphere is found.
 */
static int zipper(MODEL md, PARTICLES pts, int s_id, int workload, int TYPE)
{
  extern const int diligent;
  extern const char *fmt_watchdog_activated;
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;
  const char *fmt_ziper_failed =
    " [%s] WRN: Unable to close a zipper run."
    " Probably no available spheres left\n";
  const char *fmt_brake_zipper_at = " [%s] Braking cycle at sphere id %d\n";
  const char *fmt_invalid_id = " [%s] ERR: invalid ID (%d) at input\n";
  const int TYPE_BASE = TYPE;
  const int TYPE_BASE_FORBIDDEN = fabs(TYPE - TYPE_INCLUSION_BASE);
  // local redefinition of TYPE* variables basen on whether matix or inclusion
  // particles are to be considered
  const int TYPE_sphere_dimer     = TYPE_BASE + TYPE_SPHERE_DIMER;
  const int TYPE_sphere           = TYPE_BASE + TYPE_SPHERE;
  const int TYPE_dimer            = TYPE_BASE + TYPE_DIMER;
  const int TYPE_forbidden_sphere = TYPE_BASE_FORBIDDEN + TYPE_SPHERE;

  const unsigned short nseek_limit = 100;
  int i=0, nseek=0, step=0, s0, s1;
  int sngb_id, sngb_type, sngb_type_lock, rand_ngb, valid_ngb, next_s_id;
  int d_id, dngb_id;
  int allow_tp_sph=0;

  // Check input parameters to avoid reading from outside the array
  if(s_id < 0 || s_id >= md.Nsph){
    fprintf(stderr, fmt_invalid_id, __func__, s_id);
    return EXIT_FAILURE;
  }

  if(sph[s_id].type == TYPE_sphere_dimer){
    // Get the index of dimer for type-sphere-dimer object and id's of both
    // spheres NOTE: one of them is s_id, but who cares.
    d_id = sph[s_id].dim_ind;
    s0 = dim[d_id].sph_ind[0];
    s1 = dim[d_id].sph_ind[1];
    // Flag its spheres as type-spheres and brake dimer
    sph[s0].type = TYPE_sphere;
    sph[s1].type = TYPE_sphere;
    sph[s0].dim_ind = TYPE_INVALID;
    sph[s1].dim_ind = TYPE_INVALID;

    dim[d_id].type = TYPE_INVALID;
    dim[d_id].sph_ind[0] = TYPE_INVALID;
    dim[d_id].sph_ind[1] = TYPE_INVALID;
  }else if(sph[s_id].type == TYPE_forbidden_sphere){
    // Exit if the inclusion or matrix sphere (depending on the current TYPE),
    // has been selected as starting point
    return -1;
  }

  // Loop until the two spheres meet
  while(1){
    // Master Watchdog on zipper steps, to brake infinite loop when zipper ran
    // on single free sphere in the system.
    if(step > 10 * diligent){
      fprintf(stdout, fmt_ziper_failed , __func__);
      fprintf(stdout, fmt_brake_zipper_at, __func__, s_id);
      fprintf(stderr, fmt_watchdog_activated, __func__);
      break;
    }

    // Select random neighbor of s_id
    nseek=0;
    do{
      valid_ngb = 0;
      // Select neighbor index randomly
      rand_ngb = get_random_range12();

      // Get neighbor id and type
      sngb_id = sph[s_id].ngb[rand_ngb];
      sngb_type = sph[sngb_id].type;
      sngb_type_lock = sph[sngb_id].type_locked;

      // If enough steps passed OR if type-sphere-dimer is NOT on the list,
      // allow type-sphere selection
      allow_tp_sph = (
        (sngb_type == TYPE_sphere) && (sngb_type_lock == 0) &&
        ((step > workload) || (nseek > nseek_limit)) ) ? 1 : 0;

      // Select type-sphere-dimer or, if enough steps passed allow regular
      // sphere selection
      valid_ngb = (allow_tp_sph || (sngb_type == TYPE_sphere_dimer)) ? 1 : 0;

      // Security check not to select
      // type-(inclusion|matrix)-sphere spheres (depending on selected TYPE)
      valid_ngb = (sngb_type == TYPE_forbidden_sphere) ? 0 : valid_ngb;

      // Watchdog: count attempts to select valid neighbor
      nseek++;
    }while(!valid_ngb);

    // Brake neighboring dimer or connect free spheres and brake the cycle
    if(sngb_type == TYPE_sphere_dimer){
      // Get dimer index to whom sngb_id belongs to
      dngb_id = sph[sngb_id].dim_ind;

      // Get the index of the other sphere of dimer 'dngb_id'
      next_s_id = (dim[dngb_id].sph_ind[0] == sngb_id)
        ? dim[dngb_id].sph_ind[1] : dim[dngb_id].sph_ind[0];

      // Flag 'next_s_id' as type-sphere ...
      sph[next_s_id].type = TYPE_sphere;
      sph[next_s_id].dim_ind = TYPE_INVALID;

      // ...  and make dimer from the remaining two spheres
      sph[s_id].type = TYPE_sphere_dimer;
      sph[s_id].dim_ind = dngb_id;

      // Ammend the indices in data structures
      dim[dngb_id].sph_ind[0] = s_id;
      dim[dngb_id].sph_ind[1] = sngb_id;

      // Switch to created free sphere
      s_id = next_s_id;

      // Count steps
      step++;

    }else if(sngb_type == TYPE_sphere){
      // Make final dimer
      make_dimer(pts, md, s_id, sngb_id, TYPE_dimer);
      // Brake the master while loop and finish procedure
      break;
    }
  } // Infinite 'while' loop

  // Recalculate centers of mass and orientations for type-diemrs
  for(i=0; i<md.Ndim; i++){
    if(dim[i].type == TYPE_dimer){
      if(update_dimer_parameters(md, pts, i) == EXIT_FAILURE){
        fprintf(stderr, fmt_internal_call_failed, __func__);
      }
    }
  }
  return step;
}

/*
 * make_dimer(dim, sph, s1, s2, type_x)
 * Make dimer from s1,s2 spheres, mark it as type_x and update all relevant
 * table data
 */
static void make_dimer(PARTICLES pts, MODEL md, int s1, int s2, int type_tgt){
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;
  const char *fmt_no_free_dimer_slots =
    " [%s] ERR: No free dimers slots found in array.\n";
  const char *fmt_failed_to_connect =
    " Cannot make dimer from %d and %d spheres\n";
  const char *fmt_data_corruption =
    " [%s] ERR: Inconsistent data in dimer structure.\n"
    "  * recived sphere IDs: %d %d\n";
  int i = 0;
  int valid_ids = 1;

  // Find first free slot in dimer array
  while(dim[i].type != TYPE_INVALID){
    i++;
    // Brake the loop if no free dimers are found - this should not happen
    if(i == md.Ndim){
      i = -1;
      fprintf(stderr, fmt_no_free_dimer_slots, __func__);
      fprintf(stderr, fmt_failed_to_connect, s1, s2);
      break;
    }
  }

  // Validate obtained parameters
  if(s1 < 0 || s2 < 0 || s1 >= md.Nsph || s2 >= md.Nsph){
    valid_ids = 0;
  }

  if(i >= 0 && dim[i].type == TYPE_INVALID && valid_ids == 1){
    // Update sphere data
    // NOTE: This should be an inverse of type_tgt (egz. type_sphere_tgt),
    // but as long as they are kept equal in global settings, it is OK.
    sph[s1].type = type_tgt;
    sph[s1].dim_ind = i;

    sph[s2].type = type_tgt;
    sph[s2].dim_ind = i;

    // Update dimer data
    dim[i].type = type_tgt;
    dim[i].sph_ind[0] = s1;
    dim[i].sph_ind[1] = s2;
    update_dimer_parameters(md, pts, i);
  }else{
    fprintf(stderr, fmt_data_corruption, __func__, s1, s2);
  }
}

/*
 * find_valid_cluster()
 *
 * Find a pair of molecules oriented such that they can be flipped to change
 * their orientation
 */

static void find_valid_cluster(MODEL md, PARTICLES pts, int vclust[2]){
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;
  int rn_id0, rn_id1;   // array indexes generated randomly
  int a0, a1, b0, b1;   // sphere indexes of first and second dimer resp.
  int d0, nd0, nd1;     // selected dimer index and neighboring dimer ids.
  int candidate0 = 0, candidate1 = 0; // Valid configuration flags

  vclust[0] = -1;
  vclust[1] = -1;
  // Do the following untill valid dimer cluster is located
  do{
    // Select a random dimer of type-dimer
    do{ // Be aware of infinite loop
      d0 = (int) (u_RNG() * (md.Ndim+1));
      d0 = (d0 >= md.Ndim) ? md.Ndim : d0;
    }while(dim[d0].type != TYPE_DIMER);

    // Get index of spheres for the selected dimer
    a0 = dim[d0].sph_ind[0];
    a1 = dim[d0].sph_ind[1];

    // Randomly get the index of neighboring spheres
    rn_id0 = get_random_range12();
    rn_id1 = get_random_range12();

    // Get indexes of respective neighbors of atoms a0,a1
    b0 = sph[a0].ngb[rn_id0];
    b1 = sph[a1].ngb[rn_id1];

    // Clear flags
    candidate0 = candidate1 = 0;
    nd0 = nd1 = -1;

    // Make sure b0 is not the other atom of dimer d0 and it is a dimer sphere
    if(b0 != a1 && sph[b0].type == TYPE_SPHERE_DIMER){
      nd0 = sph[b0].dim_ind;
      candidate0 = chk_dim_config(pts, md.box, d0, nd0);
    }

    // Make sure b1 is not the other atom of dimer d0 and it is a dimer sphere
    if(b1 != a0 && sph[b1].type == TYPE_SPHERE_DIMER){
      nd1 = sph[b1].dim_ind;
      candidate1 =  chk_dim_config(pts, md.box, d0, nd1);
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

/* # SUB-SEC ########## INSERTING DIMERS - UTILITIES ######################## */

/*
 * chk_dim_config(dim,sph,d0,d1)
 *
 * Check the orientation of d0 with relation to d1 and return the index
 * of possible transformation
 */
static int chk_dim_config(PARTICLES pts, double box[3], int d0, int d1){
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;

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
int test_dimer_distribution(DIM3D *dim, int od[6], int Ndim)
{
  int i, j;
  int o_ind;
  int od_local[6]={0,0,0,0,0,0};

  for(i=0; i<Ndim; i++){
    if(dim[i].type == TYPE_DIMER){
      o_ind = check_dimer_direction(dim, i);
      if(o_ind != -1){
        od_local[o_ind]++;
      }else{
        fprintf(stderr, fmt_internal_call_failed, __func__);
        return EXIT_FAILURE;
      }
    }
  }

  for(j=0; j<6; j++){
    od[j] = od_local[j];
  }
  return EXIT_SUCCESS;
}

/*
 * display_dimer_distribution()
 *
 * Display current partition of dimer orientations
 */
static void display_dimer_distribution(int od[6])
{
  int j;

  for(j=0; j<6; j++){
    fprintf(stdout," %3d  ",od[j]);
  }
  fprintf(stdout,"\n");
}

/*
 * validate_distrib(od, nd, flip)
 *
 * checks the distribution of dimer orientations and brakes the cluster moves
 * if good-enough state is achieved. Returns '1' if the conditions are met and
 * zero otherwise.
 */
static int validate_distrib(int od[6], int Ndim, int flip)
{
  int i, round_up=0;
  int dimer_per_dir;
  int exit_code = 1; // Assume the structure is good ;)
  extern const int display_interval;
  extern const int min_flips;
  extern const int max_flips;

  // Set exit code criteria depending on the number of dimers in the system
  if(Ndim % 6 != 0 && flip < max_flips) round_up = 1;
  dimer_per_dir = ((int) (((double) Ndim) / 6e0)) + round_up;

  for(i=0; i<6; i++){
    exit_code *= ( od[i] <= dimer_per_dir ) ? 1 : 0;
  }

  if(exit_code == 1 || flip % display_interval == 0){
    fprintf(stdout," Step %9d distribution :", flip);
    for(i=0; i<6; i++){
      fprintf(stdout," %3d  ",od[i]);
    }
    fprintf(stdout,"cond.: <= %d\n",dimer_per_dir);
  }

  // Force a minimum number of flips to be made
  if(flip < min_flips) exit_code = 0 ;

  return exit_code;
}

/*
 * flip_dimers(dim, sph, box, od, d0, d1)
 *
 * Change configuration of two specified dimers d0, d1 and update orientation
 * table
 */
static void flip_dimers(MODEL md, PARTICLES pts, int od[6], int d_ids[2])
{
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;
  int d0 = d_ids[0];
  int d1 = d_ids[1];
  int rn_id0, rn_id1;   // array indexes generated randomly
  int a0, a1, b0, b1;   // sphere indexes of first and second dimer resp.
  int o_id0 = check_dimer_direction(dim, d0);
  int o_id1 = check_dimer_direction(dim, d1);
  double dist0, dist1;  // interatomic distances
  const char *fmt_illegal_orientation_before =
    " [%s] ERR: Illegal orientation index (%d) of dimer %d before flip\n";
  const char *fmt_illegal_orientation_after =
    " [%s] ERR: Illegal orientation index (%d) of dimer %d after flip\n";

  // Decrement appropriate orientation buckets
  if(o_id0 != -1){
    od[o_id0]--;
  }else{
    fprintf(stderr, fmt_illegal_orientation_before, __func__, o_id0, d0);
  }

  if(o_id1 != -1){
    od[o_id1]--;
  }else{
    fprintf(stderr, fmt_illegal_orientation_before, __func__, o_id1, d1);
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
  dist0 = distance(sph[a0].r, sph[b0].r, md.box);
  dist1 = distance(sph[a1].r, sph[b1].r, md.box);

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
  update_dimer_parameters(md, pts, d0);
  update_dimer_parameters(md, pts, d1);

  // Increment appropriate orientation buckets
  o_id0 = check_dimer_direction(dim, d0);
  o_id1 = check_dimer_direction(dim, d1);

  if(o_id0 != -1){
    od[o_id0]++;
  }else{
    fprintf(stderr, fmt_illegal_orientation_after, __func__, o_id0, d0);
  }

  if(o_id1 != -1){
    od[o_id1]++;
  }else{
    fprintf(stderr, fmt_illegal_orientation_after, __func__, o_id1, d1);
  }
}

/* vim: set tw=80 ts=2 sw=2 et: */
