#ifndef _TERMINATORS_H
#define _TERMINATORS_H

/*
 * Declarations of terminator functions (ran when exit signal is received)
 */

void free_dimers  (int stat, void *ptr);
void free_spheres (int stat, void *ptr);
void free_slits   (int stat, void *ptr);
void free_channels(int stat, void *ptr);

#endif
/* vim: set tw=80 ts=2 sw=2 et: */