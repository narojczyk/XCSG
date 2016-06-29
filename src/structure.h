#ifndef _STRUCTURE_H
#define _STRUCTURE_H

int validate_distrib(int od[6], int nd1, int step);
int zipper(DIM3D *dim, SPH *sph, double box[3], int nd, int sph_ind, int ms);

void flip_dimers(DIM3D *dim, SPH *sph, double box[3], int od[6], int d0, 
                 int d1);
void make_channel(
  DIM3D *dim, SPH *sph, double c[3], double cr, double box[3], double tr[3],
  double csd, int nd);
void make_slit(DIM3D *dim, SPH *sph, double sth, int c[3], int nd);

#endif
