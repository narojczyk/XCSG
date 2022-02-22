#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;         // Sphere type
  int ngb[12];      // Neighbors list
  int dim_ind;      // Index of dimer the sphere belongs to
  int lattice_ind[3]; //sphere lattice coordinates
  double r[3];      // Position vector
  double d;         // Sphere diameter
} SPH;

typedef struct
{
  int type;         // Dimer type
  int sph_ind[2];   // Indexes of spheres forming dimer
  int ngb[22][2];   // Neighbors list
  double R[3];      // Dimer position
  double O[3];      // Dimer orientation (unit vector)
  double L;         // Dimer length (dist. between spheres' centers)
} DIM3D;

typedef struct INC
{
  int tgt_Nmer;     /* Targeted n-mer particle to be formed by spheres
                       of this inclusion */
  double os[3];     // Inclusion offset from the base atom in the structure
  double nm[3];     // Normal vector of the inclusion
  double radius;    // Inclusion radius (for spherically symmetric inclusions)
  double thickness; // Inclusion thickness (in direction orthogonal to 'nm')
  double sph_d;     // Diameter of spheres in the inclusion
} INC;

typedef struct
{
  unsigned long int seed; // Seed for MT19937 p.r.n.g.
  int cells[3];           // Number of f.c.c. cells on edge of the box
  int first;              // Index of first structures to generate
  int last;               // Index of last structure to generate
  int mk_channel;         // Bolean [0|1] flag to enable channel inclusion
  int mk_slit;            // Bolean [0|1] flag to enable slit inclusion
  int mk_dimers;          // Bolean [0|1] flag to connect free spheres to dimers
  int mk_inc_dimers;      /* Bolean [0|1] flag to connect inclusion spheres
                             to dimers */
  int num_channels;       // Number of channels described in cfg_channels
  int num_slits;          // Number of slits described in cfg_slits
  char cfg_channels[41];  // File name for the channels' parameters
  char cfg_slits[41];     // File name for the slits' parameters
  char symmetry[8];       // Symmetry of the generated structure
} CONFIG;

typedef struct
{
  int Nsph;       // Number of monomers forming system (based on structure)
  int Ndim;       // Maximum possible number of dimers (based on Nsph)
  int mtrx_sph;   // Current number of single spheres in matrix
  int mtrx_dim;   // Current number of dimers in matrix (not in inclusion)
  int incl_sph;   // Current number of spheres forming inclusion
  int incl_dim;   // Current number of dimers forming inclusion
  double box[3];  // Size of the container for the structure  //old: box_edge
} MODEL;

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
