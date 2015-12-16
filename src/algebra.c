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

/*
 * vdotM(v, m, w)
 *
 * vector 'v=(v0,v1,v2)' and matrix 'm(3x3)' product. The result is stored in
 * 'w=(w0,w1,w2)'
 * |w0|                |m00 m01 m02|
 * |w0| = |v0,v1,v2| . |m10 m11 m12|
 * |w2|                |m20 m21 m22|
 */
void vdotM(double v[3], double m[3][3], double w[3])
{
  w[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  w[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  w[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/*
 * matrix_product(m, n, r)
 *
 * Multiplies matrices m and n storing the result in matrix r.
 */
void matrix_product(double m[3][3], double n[3][3], double r[3][3])
{
  int i, j, k;

  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){
      r[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        r[i][j] += m[i][k] * n[k][j];
    }
  }
}


/*
 * scalar_matrix_product(m, c)
 *
 * Multiplies every element of matrix m by a constant c.
 */
void scalar_matrix_product(double m[3][3], double c)
{
  int i, j;

  for (i = 0; i < 3; i++){
    for (j = 0; j < 3; j++){
      m[i][j] *= c;
    }
  }
}


/*
 * invert_matrix(i, m)
 *
 * Calculetes an invert matrix of m and places the result in i.
 */
void invert_matrix(double im[3][3], double m[3][3])
{
  double inv_detM = 1.0 / det_matrix(m);

  // macierz dopelnien algebraicznych
  // 1st row
  im[0][0] =  (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
  im[0][1] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
  im[0][2] =  (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  // 2nd row
  im[1][0] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]);
  im[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
  im[1][2] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]);

  // 3rd row
  im[2][0] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]);
  im[2][1] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]);
  im[2][2] =  (m[0][0] * m[1][1] - m[0][1] * m[1][0]);

  // macierz dolaczona (transpozycja macierzy dopelnien alg.)
  transpose_matrix(im);

  // macierz odwrotna (macierz dolaczona / wyznacznik mac. podst.)
  scalar_matrix_product(im, inv_detM);
}


/*
 * transpose_matrix(m)
 *
 * Transposes in-place a given matrix m.
 */
void transpose_matrix(double m[3][3])
{
  double a[3][3];
  int i, j;

  memcpy(a, m, 3 * sizeof *m);
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      m[j][i] = a[i][j];
}

/*
 * det_matrix(m)
 *
 * Retruns the determinant of matrix m.
 */
double det_matrix(double m[3][3])
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) +
    m[0][1] * (m[1][2] * m[2][0] - m[1][0] * m[2][2]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/* vim: set tw=80 ts=2 sw=2 et: */
