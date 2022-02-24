/*
 * algebra.h
 *
 * Header file containing function prototypes of commonly used math operations.
 */

#ifndef _ALGEBRA_H
#define _ALGEBRA_H

void vnorm(double v[3]);
void vcrossu(double v[3], double u[3], double w[3]);
void vector(double p[3], double q[3], double pq[3], double box[3]);
double distance(double p[3], double q[3], double box[3]);
double distance_absolute(double p[3], double q[3]);
double vmodule(double v[3]);
double vdotu(double v[3], double u[3]);

#endif
