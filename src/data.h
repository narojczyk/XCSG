#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;         // Sphere type
  int ngb[12];      // Neighbors list
  double r[3];      // Position vector
  double d;         // Sphere diameter
} SPH;


#endif
/* vim: set tw=80 ts=2 sw=2 et: */
