#ifndef _STRUCTURE_H
#define _STRUCTURE_H


int refine_dimer_distribution(DIM3D *dim, SPH *sph, MODEL md, int od[6]);
int set_structure(CONFIG cf, MODEL md, SPH *sph);
int sph_assign_lattice_indexes( SPH *sph, int ns);
int test_dimer_distribution(DIM3D *dim, int od[6], int nd);

void introduce_random_dimers(DIM3D *dim, SPH *sph, MODEL md);
void introduce_dimers_by_zipper(DIM3D *dim, SPH *sph, MODEL md);
void make_channel(MODEL md, SPH *sph, INC *inc);
void make_slit(MODEL md, SPH *sph, INC *inc);


#endif
