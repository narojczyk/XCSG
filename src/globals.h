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

// config file variables
unsigned long int i_seed;       // Seed for MT19937 r.n.g.
int i_edge_fcc_N = 0;           // Number of f.c.c. cells on edge of the box
int i_normal[3] = {0, 0, 0};    // Normal vector describing channel or slit
int i_iDCfrom = -1;             // Index of input DC structure (if <0,no input)
int i_iDCto = -1;               // Index of last DC structure
int i_make_channel = 0;         // Bolean flag to enable nano-channel [0|1]
int i_make_slit = 0;            // Bolean flag to enable nano-slit [0|1]
double i_channel_R = 0e0;       // Nano-channel radius in sigma units

// general purpose global variables
char *prog_name;


#endif

/* vim: set tw=80 ts=2 sw=2 et: */
