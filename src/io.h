#ifndef _IO_H
#define _IO_H

void display_configuration_summary(CONFIG cfg, MODEL md, SLI *slits,
                                   CHA *channels);
void display_stats(MODEL md);

int exp_str_data(CONFIG cf, MODEL md, DIM3D *dim, SPH *sph,
                 int ns3, int nd2, int strn);
int export_dimers(FILE *file, DIM3D *dim, int nd);
int export_spheres(FILE *file, SPH *sph, int ns);
int povray_export_spheres(FILE *file, SPH *sph, int ns);
int export_to_GLviewer(MODEL md, DIM3D *dim, SPH *sph, int strn);
int legacy_GLexport_dimer_type_converter(int type);
int legacy_GLexport_sphere_type_converter(int type);


#endif
