#ifndef _UTILS_H
#define _UTILS_H

int str_validate(const char *src, const char *tgt);
int count_particles_by_type(MODEL *md, SPH *sph, DIM3D *dim);
int container_dimensions(MODEL *md, CONFIG cf);
int number_of_spheres(CONFIG cf);
int find_ngb_spheres(SPH sph[], int ns, double box[3]);
int check_dimer_direction(DIM3D *dim, int i);
int count_typeX_sp_neighbours(SPH *sph, int x, int id, int ns);
int find_critical_sphere(SPH *sph, int type_x, int ns);
int draw_sphere_typeX(SPH *sph, int x, int ns);
int draw_ngb_sphere_typeX(SPH *sph, int x, int sph_ind);
int update_dimer_parameters(MODEL md, DIM3D *dim, SPH *sph, int d);

void memory_clean_spheres(SPH *sph, int nd);
void memory_clean_dimers(DIM3D *dim, int nd);
void memory_clean_inclusion(INC *inc, int n);

void bouble_sort_double(double *array, int s, int ascending);
void init_RNG(unsigned long int s);

double u_RNG();

#endif
