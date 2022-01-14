#ifndef _GLOBALS_H
#define _GLOBALS_H

/*
 * Declarations of global variables
 */

// constants
const double zero = (double) 0;
const double one = (double) 1;
const double two = (double) 2;

// particle types identifiers
const int TYPE_INVALID = -1;
const int TYPE_SPHERE = 1;
const int TYPE_SPHERE_DIMER = 2;
const int TYPE_DIMER = TYPE_SPHERE_DIMER;
const int TYPE_INCLUSION_SPHERE = 100 + TYPE_SPHERE;
const int TYPE_INCLUSION_DIMER = 100 + TYPE_DIMER;

// zipper work length
const int lazy = 10;
const int diligent = 1000000;

// DC structure
const int run_zipper = 1000000;
const int display_interval = 1000000;
const int max_flips = 100000000;
//  * optimise DC structure if this % of dimers is present
const double minimal_dimer_qty = 2e-1;

// general purpose global variables
char *prog_name;

// message strings
// *** errors
const char *fmt_internal_call_failed = " [%s] ERR: internall call failed\n";
const char *fmt_open_failed = " [%s] ERR: cannot open file: %s\n";
const char *fmt_writting_failed = "  [%s] ERR: writing %s data failed\n";

// *** warnings

// *** notifications
const char *fmt_write_notify = " Writting %-6s data to file %s\n";

// *** additional debuging messages
const char *fmt_dbg_opening_file = " [%s] DBG opening file: %s\n";

// *** I/O formats

#endif
/* vim: set tw=80 ts=2 sw=2 et: */
