#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;         // Sphere type
  int ngb[12];      // Neighbors list
//   TODO: possibly change the following to pointer to DIM array element
  int dim_ind;      // Index of dimer the sphere belongs to
  int lattice_ind[3]; //sphere lattice coordinates
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
  double os[3];     // Channel offset from the base atom in the structure
  double nm[3];     // Vector parallel to the channel axis
  double radius;    // The radius of the channel
  double sph_d;     // Diameter of spheres in the channel
} CHA;

typedef struct
{
  double os[3];     // Channel offset from the base atom in the structure
  double nm[3];     // Vector orthogonal to the plane
  double thickness; // The thickness of the plane
  double sph_d;     // Diameter of spheres in the plane
} SLI;

typedef struct
{
  unsigned long int seed; // Seed for MT19937 p.r.n.g.
  int fcc_cells[3];       // Number of f.c.c. cells on edge of the box
  int first;              // Index of first structures to generate
  int last;               // Index of last structure to generate
  int mk_channel;         // Bolean [0|1] flag to enable channel inclusion
  int mk_slit;            // Bolean [0|1] flag to enable slit inclusion
  int mk_dimers;          // Bolean [0|1] flag to connect free spheres to dimers
  int num_channels;       // Number of channels described in cfg_channels
  int num_slits;          // Number of slits described in cfg_slits
  char cfg_channels[41];  // File name for the channels' parameters
  char cfg_slits[41];     // File name for the slits' parameters
  char symmetry[8];       // Symmetry of the generated structure
} CONFIG;

typedef struct
{
  int Nmon;     // Number of monomers forming system (based on structure) //old: Ns
  int Ndim;     // Maximum possible number of dimers (based on Nmon) // old: Nd

} MODEL;

typedef struct
{

} PTINC;

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
