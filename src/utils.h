#ifndef _UTILS_H
#define _UTILS_H

int find_ngb_spheres(SPH sph[], int ns, double box[3]);
int brake_dimers(DIM3D *dim, SPH *sph, int nd);
int sph_set_fcc( SPH *sph, int ns, int fcc[3]);

int check_dimers_configuration(DIM3D *dim, SPH *sph, double box[3], 
                              int d0, int d1);
int check_dimer_direction(DIM3D *dim, int i);

void memory_clean_spheres(SPH *sph, int nd);
void memory_clean_dimers(DIM3D *dim, int nd);
void update_sphere_positions(DIM3D *dim, SPH *sph, double box[3], int d);
void update_dimer_parameters(DIM3D *dim, SPH *sph, double box[3], int d);
void bind_spheres_to_dimers(DIM3D *dim, SPH *sph, int nd);
void init_MT19937(unsigned long int s);
void init_RNG(unsigned long int s);
void dimer_distribution(DIM3D *dim, int od[6], int nd);
void find_valid_cluster(DIM3D *dim, SPH *sph, double box[3], int nd, 
                        int vclust[2]);

double u_RNG();

#endif
