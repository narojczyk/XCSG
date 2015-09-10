#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;         // Sphere type
  int ngb[12];      // Neighbors list
  double r[3];      // Position vector
  double d;         // Sphere diameter
} SPH;

typedef struct
{
  int type;         // Dimer type
//   TODO: possibly change the following to array of pointers to SPH array
  int sph_ind[2];   // Indexes of spheres forming dimer
  int ngb[22][2];   // Neighbors list
  double R[3];      // Dimer position
  double O[3];      // Dimer orientation (unit vector)
  double L;         // Dimer length (dist. between spheres' centers)
} DIM3D;

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
