#ifndef _DATASTRUCTURES_H
#define _DATASTRUCTURES_H

typedef struct Walker_tag
{
//   struct Walker_tag *point_to;
  int status;     // Walker status [0|1] (non-active/active)
  int type;       // Walker type [1|2] (searcher|target)
  int u;          // exponent of the power-law distribution
  int steps;      // The number of steps generated for the walker
  int hits;       // The number of encounters
  int step_completed; // BOOL flag [0|1] that signals completed step
  double r0[2];   // Initial position vector (for <R^2>)
  double r[2];    // Position vector
  double d[2];    // unit vector with the movement direction
  double l0;      // step length scaling factor
  double l;       // step length
  double tot_L;   // total distance traveled
  double R;       // Radius
  double v;       // velocity
  double t;       // time left to travel distance 'l' at velocity 'v'
} Walker;


#endif
/* vim: set tw=80 ts=2 sw=2 et: */
