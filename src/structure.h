#ifndef _STRUCTURE_H
#define _STRUCTURE_H

int set_structure(CONFIG cf, MODEL md, SPH *sph);
int set_fcc(SPH *sph, int ns, int cells[3]);
int validate_distrib(int od[6], int nd1, int step);
int zipper(MODEL md, DIM3D *dim, SPH *sph, int s_id, int workload);
int check_dimers_configuration(DIM3D *dim, SPH *sph, double box[3], 
                              int d0, int d1);
int sph_assign_lattice_indexes( SPH *sph, int ns);

void make_dimer(DIM3D *dim, SPH *sph, MODEL md, int s1, int s2);
void flip_dimers(MODEL md, DIM3D *dim, SPH *sph, int od[6], int d_ids[2]);
void make_channel(MODEL md, DIM3D *dim, SPH *sph, INC *inc);
void make_slit(MODEL md, DIM3D *dim, SPH *sph, INC *inc);
void find_valid_cluster(MODEL md, DIM3D *dim, SPH *sph, int vclust[2]);
int test_dimer_distribution(DIM3D *dim, int od[6], int nd);
void display_dimer_distribution(int od[6]);

#endif
