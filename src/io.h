#ifndef _IO_H
#define _IO_H

int export_dimers(FILE *file, DIM3D *dim, int nd);
int export_spheres(FILE *file, SPH *sph, int ns);
int export_to_GLviewer(SPH *sph, double box_x, int sn, int ns);

int load_dcsgen(FILE *file, DIM3D *dim, double box_x, int nd);

#endif