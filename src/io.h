#ifndef _IO_H
#define _IO_H

int export_dimers(FILE *file, DIM3D *dim, int nd);
int export_spheres(FILE *file, SPH *sph, int ns);
int export_to_GLviewer(DIM3D *dim, SPH *sph, double box[3], int sn, int ns, 
                       int nd);

int load_dcsgen(FILE *file, DIM3D *dim, double box[3], int nd);

#endif