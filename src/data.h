#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;         // Sphere type
  int ngb[12];      // Neighbors list
//   TODO: possibly change the following to pointer to DIM array element
  int dim_ind;      // Index of dimer the sphere belongs to
  double r[3];      // Position vector
  double d;         // Sphere diameter
} SPH;

typedef struct
{
  int type;         // Dimer type
//   TODO: possibly change the following to array of pointers to SPH array elem.
  int sph_ind[2];   // Indexes of spheres forming dimer
  int ngb[22][2];   // Neighbors list
  double R[3];      // Dimer position
  double O[3];      // Dimer orientation (unit vector)
  double L;         // Dimer length (dist. between spheres' centers)
} DIM3D;

typedef struct
{
  double offset[3]; // Channel offset from the base atom in the structure
  double normal[3]; // Vector parallel to the channel axis
  double radius;    // The radius of the channel
  double sph_d;     // Diameter of spheres in the channel
} CHA;

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
