/*
 * config.h
 *
 * Header file containing macrodefinitions for code pre-processing
 */

#ifndef _CONFIG_H
#define _CONFIG_H

// Uncoment the following to save the system's evolution to trajectory files
// #define DATA_VISGL_OUTPUT

// Uncoment the following to export mean square displacement for particles
// #define PUSH_MSD

// Uncoment the following to use the 64bit implementation of MT 19937 r.n.g.
// #define USE_64BIT_MT19937

// Uncoment the following to enable permanent kill mode (target are not
// recreated at random positions after hit)
// #define PERMANENT_KILL

// Compile program with capture-the-flag-mode, i.e. targets must reach specified
// area. This mode assumes permanentKill=1
// #define CAPTURE_THE_FLAG

// Uncoment the following to enable debugging output
// #define DEBUG_MODE

#endif
