#ifndef _UTILS_H
#define _UTILS_H

int count_particles_by_type(MODEL *md, SPH *sph, DIM3D *dim);
int container_dimensions(MODEL *md, CONFIG cf);
int number_of_spheres(CONFIG cf);
int find_ngb_spheres(SPH sph[], int ns, double box[3]);
int brake_dimers(DIM3D *dim, SPH *sph, int nd);
int check_dimer_direction(DIM3D *dim, int i);
int count_typeX_spheres(SPH *sph, int x, int ns);
int count_typeX_dimers(DIM3D *dim, int x, int nd);
int count_typeX_sp_neighbours(SPH *sph, int x, int id);
int find_critical_sphere(SPH *sph, int ns);
int draw_sphere_typeX(SPH *sph, int x, int ns);
int draw_ngb_sphere_typeX(SPH *sph, int x, int sph_ind);

void memory_clean_spheres(SPH *sph, int nd);
void memory_clean_dimers(DIM3D *dim, int nd);
void memory_clean_slits(SLI *sli, int nsl);
void memory_clean_channels(CHA *cha, int nch);
void update_sphere_positions(DIM3D *dim, SPH *sph, double box[3], int d);
int update_dimer_parameters(MODEL md, DIM3D *dim, SPH *sph, int d);
void bind_spheres_to_dimers(DIM3D *dim, SPH *sph, int nd);
void init_MT19937(unsigned long int s);
void init_RNG(unsigned long int s);
void bouble_sort_double(double *array, int s, int ascending);

double u_RNG();

#endif
