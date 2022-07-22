#ifndef _STRUCTURE_H
#define _STRUCTURE_H


int introduce_dimers_by_zipper(DIM3D *dim, SPH *sph, MODEL *md, int TYPE);
int refine_dimer_distribution(DIM3D *dim, SPH *sph, MODEL md, int od[6]);
int set_structure(CONFIG cf, MODEL md, SPH *sph);
int test_dimer_distribution(DIM3D *dim, int od[6], int nd);

void find_smallest_coordinates(SPH *sph, MODEL *md);
void introduce_random_dimers(DIM3D *dim, SPH *sph, MODEL *md, int type_src,
                             int type_tgt);
void make_channel(MODEL md, INC *inc, SPH *sph, DIM3D *dim);
void make_slit(MODEL md, INC *inc, SPH *sph, DIM3D *dim);


#endif
