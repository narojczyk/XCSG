#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct 
{
  int type;       // Walker type [1|2] (searcher|target)
  int ngb[12];
  double r[3];   // Initial position vector (for <R^2>)
  double d;
} SPH;


#endif
/* vim: set tw=80 ts=2 sw=2 et: */
