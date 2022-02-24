/*
 * algebra.c
 *
 * Definitions of commonly used matrix operations.
 *
 * We are dealing only with matrices of rank 3, thus dimension parameters is
 * NOT explicitly specified in any of functions below.
 */
#include <string.h>
#include <math.h>
#include "algebra.h"


extern const double zero;
extern const double one;

/*
 * vdotu(v,u)
 *
 * Returns the dot product of two vectors v*u
 */
double vdotu(double v[3], double u[3])
{
  return v[0]*u[0] + v[1]*u[1] + v[2]*u[2];
}

/*
 * vector(p,q,v, box)
 * Calculates the vector v in 3D real box,taking into account the periodic
 * boundary conditions. v is a vector from point p to point q.
 */
void vector(double p[3], double q[3], double pq[3], double box[3])
{
  int i;

  for (i=0; i<3; i++){
    // Calculat distance between p and q
    pq[i] = q[i] - p[i];

    // Apply periodic boundaries
    pq[i] = pq[i] - box[i] * round( pq[i]/box[i] );
  }
}

/*
 * distance(p, q, box)
 *
 * Computes the distance between points p and q in real cubic box with
 * edge lenght equal to box
  * Takes into account periodic boundary conditions (returns the distance of
 * q with reference to p -- not that the order makes any difference :D )
 * NOTE: Works ONLY for cubic boxes centered around the axes origin.
 */
double distance(double p[3], double q[3], double box[3])
{
  int i;
  double pq[3];

  for (i=0; i<3; i++){
    // Calculat absolute distance between p and q
    pq[i] = q[i] - p[i];

    // Apply periodic boundaries
    pq[i] = pq[i] - box[i] * round( pq[i]/box[i] );
  }

  return sqrt(pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2]);
}

/*
 * distance_absolute(p, q)
 *
 * Computes the absolute distance between points p and q in real space
 */
double distance_absolute(double p[3], double q[3])
{
  int i;
  double pq[3];

  for (i=0; i<3; i++){
    // Calculat absolute distance between p and q
    pq[i] = q[i] - p[i];
  }

  return sqrt(pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2]);
}

/*
 *
 * vmodule(v)
 *
 * Returns the length of the vecotr v
 */
double vmodule(double v[3])
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/*
 * vnorm(v)
 *
 * Norms 'v' to unit vector
 */
void vnorm( double v[3] )
{
  double vl = sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
  if( vl > 0e0 ){
    v[0] = v[0]/vl;
    v[1] = v[1]/vl;
    v[2] = v[2]/vl;
  }else{
    v[0] = zero;
    v[1] = zero;
    v[2] = zero;
  }
}

/*
 * vcrossu(v,u,w)
 *
 * Calculates the cross product of v and u vectors and stores it in w.
 * It calculates v x u.
 */
void vcrossu(double v[3], double u[3], double w[3])
{
  w[0] =   v[1]*u[2] - v[2]*u[1];
  w[1] = -(v[0]*u[2] - v[2]*u[0]);
  w[2] =   v[0]*u[1] - v[1]*u[0];
}

/* vim: set tw=80 ts=2 sw=2 et: */
