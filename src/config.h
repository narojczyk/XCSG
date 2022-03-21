/*
 * config.h
 *
 * Header file containing macrodefinitions for code pre-processing
 */

#ifndef _CONFIG_H
#define _CONFIG_H

// Uncoment the following to save the system's evolution to trajectory files
#define DATA_VISGL_OUTPUT

// Uncoment the following to enable debugging output
// #define DEBUG_MODE

// Uncoment the following to enable verbose output
// #define VERBOSE_MDOE

/* 
 * Pseudo random number generator section
 * 
 * NOTE: Pay attention that ONLY ONE of the following three options is active
 * Any TWO of the following #define's should be ALWAYS COMMENTED OUT
 */
// Uncoment the following to use the standard dran48 r.n.g.
// #define PRNG_DRAND48

// Uncoment the following to use the 32bit implementation of MT 19937 r.n.g.
#define PRNG_32BIT_MT19937

// Uncoment the following to use the 64bit implementation of MT 19937 r.n.g.
// #define PRNG_64BIT_MT19937

#endif
