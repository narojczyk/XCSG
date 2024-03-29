#ifndef _IO_H
#define _IO_H
#include "data.h"

void show_cfg_summary(CONFIG cfg, MODEL md, INC *slits, INC *channels,
                      INC *clusters);
void show_particle_stats(MODEL md, CONFIG cfg);
int export_structure(CONFIG cf, MODEL md, PARTICLES pts, int strn);

// These functions may return NULL
FILE* try_open_to_read(const char *file);

// These functions can close program on error
FILE* open_to_read(const char *file);
FILE* open_to_write(const char *file);

#endif
