#ifndef _UTILS_H
#define _UTILS_H

int find_ngb_spheres(SPH sph[], int ns, double box[3]);
int brake_dimers(DIM3D *dim, SPH *sph, int nd);
int sph_set_fcc( SPH *sph, int ns, int fcc[3]);

void memory_clean_spheres(SPH *sph, int nd);
void memory_clean_dimers(DIM3D *dim, int nd);
void make_channel(
  DIM3D *dim, SPH *sph, int c[3], double cr, double box[3], int nd);
void make_slit(DIM3D *dim, SPH *sph, double sth, int c[3], int nd);
void update_sphere_positions(DIM3D *dim, SPH *sph, double box[3], int d);
void bind_spheres_to_dimers(DIM3D *dim, SPH *sph, int nd);
void init_MT19937(unsigned long int s);

#endif
