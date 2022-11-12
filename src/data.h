#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;         // Sphere type
  int type_locked;  // Boolean [0|1] flag to enable channel inclusion
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

typedef struct{
  SPH *spheres;     // Pointer to the array with spheres
  DIM3D *dimers;    // Pointer to the array with dimers
} PARTICLES;

typedef struct INC
{
  struct INC *next; // Pointer to next instance of INC (eg., if part of a set)
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
  int mk_clusters;        // Boolean [0|1] flag to enable cluster inclusion
  int mk_channel;         // Boolean [0|1] flag to enable channel inclusion
  int mk_slit;            // Boolean [0|1] flag to enable slit inclusion
  int mk_dimers;          // Boolean [0|1] flag to connect free spheres to dim.
  int mk_inc_dimers;      /* Boolean [0|1] flag to connect inclusion spheres
                             to dimers */
  int num_clusters;       // Number of cluster incl. described in cfg_clusters
  int num_channels;       // Number of channels described in cfg_channels
  int num_slits;          // Number of slits described in cfg_slits
  int rough_inclusions;   /* Boolean [0|1] flag to allow molecules to extend
                             beyond the inclusion boundary */
  double dimer_length;    /* Target length of dimers (0e0 for auto, based on
                             atomic diameters)*/
  char cfg_clusters[41];  // File name for the cluster inclusions' parameters
  char cfg_channels[41];  // File name for the channels' parameters
  char cfg_slits[41];     // File name for the slits' parameters
  char symmetry[8];       // Symmetry of the generated structure
} CONFIG;

typedef struct
{
  int Nsph;       // Number of monomers forming system (based on structure)
  int Ndim;       // Maximum possible number of dimers (based on Nsph)
  int mtrx_sph;   // Current number of single spheres in matrix
  int mtrx_dim;   // Current number of dimers in matrix
  int incl_sph;   // Current number of spheres forming inclusion
  int incl_dim;   // Current number of dimers forming inclusion
  double box[3];  // Size of the container for the structure
  double lcs[3];  // Position of a sphere with smallest coordinates
} MODEL;

typedef struct{
  double c[3];    // Coordinates of a point in 3D space (x,y,z -> 0,1,2)
} POINT;

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
