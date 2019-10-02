#ifndef _GLOBALS_H
#define _GLOBALS_H

/*
 * Declarations of global variables
 */

// constants
const double zero = (double) 0;
const double one = (double) 1;
const double two = (double) 2;

// config file variables
unsigned long int i_seed;       // Seed for MT19937 r.n.g.
int i_edge_fcc_N[3] = {0, 0, 0};// Number of f.c.c. cells on edge of the box
int i_load_DC_from_file = 1;    // OBSOLETE TODO: remove this while cleaning ini file structure
int i_iDCfrom = -1;             // Index of input DC structure (if <0,no input)
int i_iDCto = -1;               // Index of last DC structure
int i_make_channel = 0;         // Bolean flag to enable nano-channel [0|1]
int i_make_slit = 0;            // Bolean flag to enable nano-slit [0|1]
int i_fs_connect = 0;           // Bolean flag to connect free spheres [0|1]
int i_n_channels = 0;           // Number of channels described in i_Fchannels
int i_n_slits = 0;              // Number of slits described in i_Fslits
char i_Fchannels[41];           // File name for the channels' parameters
char i_Fslits[41];              // File name for the slits' parameters

// general purpose global variables
char *prog_name;

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
