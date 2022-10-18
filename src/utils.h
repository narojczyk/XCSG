#ifndef _UTILS_H
#define _UTILS_H

double u_RNG();

int check_dimer_direction(DIM3D *dim, int i);
int container_dimensions(MODEL *md, CONFIG cf);
int count_inclusions(INC *inclusions);
int count_particles_by_type(MODEL *md, PARTICLES pts);
int count_typeX_sp_neighbours(SPH *sph, int x, int id, int ns);
int draw_ngb_sphere_typeX(SPH *sph, int x, int sph_ind);
int draw_sphere_typeX(SPH *sph, int x, int ns);
int find_critical_sphere(SPH *sph, int type_x, int ns);
int find_ngb_spheres(SPH sph[], int ns, double box[3]);
int get_random_range12();
int number_of_spheres(CONFIG cf);
int sph_assign_lattice_indexes(SPH *sph, int ns);
int str_validate(const char *src, const char *tgt);
int update_dimer_parameters(MODEL md, PARTICLES pts, int d);

void adjust_dimer_length(DIM3D *dim, SPH *sph);
void init_RNG(unsigned long int s);
void memory_clean_spheres(SPH *sph, int nd);
void memory_clean_dimers(DIM3D *dim, int nd);
INC template_inclusion();

#endif
