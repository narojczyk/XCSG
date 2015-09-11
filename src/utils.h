#ifndef _UTILS_H
#define _UTILS_H

int sph_set_fcc( SPH *sph, int ns, int fcc_x);

void update_sphere_positions(DIM3D *dim, SPH *sph, double box_x, int d);
void bind_spheres_to_dimers(DIM3D *dim, SPH *sph, int nd);
void init_MT19937(unsigned long int s);

#endif
