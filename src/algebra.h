/*
 * algebra.h
 *
 * Header file containing function prototypes of commonly used matrix
 * operations.
 */
#ifndef _ALGEBRA_H
#define _ALGEBRA_H

// Matrix functions
double det_matrix(double m[3][3]);
void transpose_matrix(double m[3][3]);
void invert_matrix(double inv_m[3][3], double m[3][3]);
void scalar_matrix_product(double m[3][3], double c);
void matrix_product(double m[3][3], double n[3][3], double r[3][3]);

void vnorm(double v[3]);
void vdotM(double v[3], double m[3][3], double w[3]);
void vcrossu(double v[3], double u[3], double w[3]);
void vector(double p[3], double q[3], double pq[3], double box[3]);
double distance(double p[3], double q[3], double box[3]);
double vmodule(double v[3]);

//void box_norm(double m[3][3], double norm[3][3]);
#endif
