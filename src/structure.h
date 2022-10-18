#ifndef _STRUCTURE_H
#define _STRUCTURE_H


int introduce_dimers_by_zipper(PARTICLES pts, MODEL *md, int TYPE);
int refine_dimer_distribution(PARTICLES pts, MODEL md, int od[6]);
int set_structure(CONFIG cf, MODEL md, SPH *sph);
int test_dimer_distribution(DIM3D *dim, int od[6], int nd);

void find_smallest_coordinates(SPH *sph, MODEL *md);
void introduce_random_dimers(PARTICLES pts, MODEL *md, int type_src,
                             int type_tgt);
void make_channel(MODEL md, INC *inc, PARTICLES pts);
void make_slit(MODEL md, INC *inc, PARTICLES pts);


#endif
