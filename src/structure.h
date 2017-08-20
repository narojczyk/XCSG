#ifndef _STRUCTURE_H
#define _STRUCTURE_H

int validate_distrib(int od[6], int nd1, int step);
int zipper(DIM3D *dim, SPH *sph, double box[3], int nd, int sph_ind, int ms);
int check_dimers_configuration(DIM3D *dim, SPH *sph, double box[3], 
                              int d0, int d1);

void make_dimer(DIM3D *dim, SPH *sph, double box[3], int s1, int s2, int nd);
void flip_dimers(DIM3D *dim, SPH *sph, double box[3], int od[6], int d0, 
                 int d1);
void make_channel(DIM3D *dim, SPH *sph, double c[3], double cr, double box[3], 
                  double tr[3], double csd, int ns);
void make_slit(DIM3D *dim, SPH *sph, double box[3], double thick, double os[3],
               double ssd, double nm[3], int ns);
void find_valid_cluster(DIM3D *dim, SPH *sph, double box[3], int nd, 
                        int vclust[2]);
void test_dimer_distribution(DIM3D *dim, int od[6], int nd);
void display_dimer_distribution(int od[6]);
#endif
