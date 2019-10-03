#ifndef _IO_H
#define _IO_H

int exp_str_data(DIM3D *dim, SPH *sph, double box[3], int ns, int nd, int strn);

int export_dimers(FILE *file, DIM3D *dim, int nd);
int export_spheres(FILE *file, SPH *sph, int ns);
int export_to_GLviewer(DIM3D *dim, SPH *sph, double box[3], int strn, int ns, 
                       int nd);

#endif
