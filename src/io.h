#ifndef _IO_H
#define _IO_H

int export_spheres(FILE *file, SPH *sp_tab, int ns);
int export_to_GLviewer(SPH *sp_tab, double box_x, int ns);

int load_dcsgen(FILE *file, DIM3D *dim, double box_x, int Nmax);

#endif