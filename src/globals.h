#ifndef _GLOBALS_H
#define _GLOBALS_H

/*
 * Declarations of global variables
 */

// constants
const double zero = (double) 0;
const double one = (double) 1;
const double two = (double) 2;
const double pi = 4e0 * atan(1e0);

//config file variables
int i_edge_fcc_N = 0;           // Number of f.c.c. cells on edge of the box
int i_chanel[3] = {0, 0, 0};    // Nano-chanel direction, all zeros - no chanel
double i_chanel_R = 0e0;       // Nano-chanel radius in sigma units


#endif

/* vim: set tw=80 ts=2 sw=2 et: */
