#ifndef _IO_H
#define _IO_H

void display_configuration_summary(CONFIG cfg, MODEL md, INC *slits,
                                   INC *channels);
void display_stats(MODEL md, CONFIG cfg);
int exp_str_data(CONFIG cf, MODEL md, DIM3D *dim, SPH *sph, int strn);

// These functions can close program on error
FILE* open_to_read(const char *file);
FILE* open_to_write(const char *file);

#endif
