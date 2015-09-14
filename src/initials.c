/*
 * initials.c
 *
 * File defines the function and it auxilaries for parsing program's command
 * line arguments. It utilizes GNU getopt framework.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "config.h"
#include "checksum.h"
#include "initials.h"

extern char *prog_name;

extern unsigned long int i_seed;
extern int i_edge_fcc_N;
extern int i_channel[3];
extern int i_iDCfrom;
extern int i_iDCto;
extern double i_channel_R;

static struct option long_opts[] = {
  {"config",   required_argument, NULL, 'c'},
  {"help",     no_argument,       NULL, 'h'},
  {"version",  no_argument,       NULL, 'v'},
  {"template", no_argument,       NULL, 't'},
};


/*
 * print_greetings()
 *
 * Print author information and program info
 */
void print_greetings()
{
  fputs("\n\
  Generate a D.C. phase of dimers on the f.c.c. lattice with nano-channels\n\
  Institute of Molecular Physics, Polish Academy of Sciences\n\n", stdout);
}


/*
 * print_usage()
 *
 * Display description of command line arguments
 */
void print_usage(int status)
{
  if(status != EXIT_SUCCESS){
    fprintf(stderr, "  Try '%s --help' for more information.\n", prog_name);
  }else{
    fprintf(stdout, "  Usage: %s [OPTIONS] [CONFIG_FILE]\n\n", prog_name);
    fputs("\
    Option list:\n\
      -h, --help     print this information\n\
      -v, --version  print version and build details\n\
      -t, --template generate template config file\n", stdout);
  }
  exit(status);
}

/*
 * print_version()
 *
 * Display build info and checksums of sources used
 */
void print_version(int status)
{
  print_greetings();

  fprintf(stdout, "  Created by:\t%s\n\t\t%s\n", author0, author1);
  fprintf(stdout, "  Build by:\t%s\n", builder);
  fprintf(stdout, "  Build date:\t%s\n\n", build);

  fprintf(stdout, "  Program build-in features:\n");
  fprintf(stdout, "  * use 64bit MT19937 random num. gen.: ");
  #ifdef USE_64BIT_MT19937
    fprintf(stdout, "Yes\n");
  #else
    fprintf(stdout, "No\n");
  #endif
    
  fprintf(stdout, "  * verbose output for debugging      : ");
  #ifdef DEBUG_MODE
    fprintf(stdout, "Yes\n");
  #else
    fprintf(stdout, "No\n");
  #endif
    fprintf(stdout, "\n");

  fprintf(stdout, "  Version control:\t");
  fprintf(stdout, "  sources checksum's (SHA1)\n");
  fprintf(stdout, "  %s  %s\n", fccdcgen_c_SHA1, "fccdcgen.c");
  #ifdef USE_64BIT_MT19937
    fprintf(stdout, "  %s  %s\n", mt19937_64_h_SHA1, "mt19937_64.h");
  #else
    fprintf(stdout, "  %s  %s\n", mt19937_h_SHA1, "mt19937.h");
  #endif
  fprintf(stdout, "  %s  %s\n", config_h_SHA1, "config.h");
  fprintf(stdout, "  %s  %s\n", data_h_SHA1, "data.h");
  fprintf(stdout, "  %s  %s\n", globals_h_SHA1, "globals.h");
  fprintf(stdout, "  %s  %s\n", algebra_h_SHA1, "algebra.h");
  fprintf(stdout, "  %s  %s\n", algebra_c_SHA1, "algebra.c");
  fprintf(stdout, "  %s  %s\n", initials_h_SHA1, "initials.h");
  fprintf(stdout, "  %s  %s\n", initials_c_SHA1, "initials.c");
  fprintf(stdout, "  %s  %s\n", io_h_SHA1, "io.h");
  fprintf(stdout, "  %s  %s\n", io_c_SHA1, "io.c");
  fprintf(stdout, "  %s  %s\n", utils_h_SHA1, "utils.h");
  fprintf(stdout, "  %s  %s\n", utils_c_SHA1, "utils.c");

  exit(status);
}

/*
 * generate_template_config(status)
 *
 * Function to generate a template config file for the program.
 * Functins terminates with exit code from 'status'.
 */
void generate_template_config(int status)
{

  FILE  *f;
  char *template = strcat(prog_name,".ini.sed");

  // Opening file
  if ((f = fopen(template,"w")) == NULL ) {
    fprintf(stderr,
      "  [%s]: error: failed to open file %s\n", __func__, template);
    exit(1);
  }

  fprintf(f, "MT19937 seed            : LUINT\n");
  fprintf(f, "Number of edge fcc cells: INT\n");
  fprintf(f, "Nano-channel direction   : INT_h INT_k INT_l\n");
  fprintf(f, "Chanel radius [sigma]   : DOUBLE\n");
  fprintf(f, "Load initial DC struct. : INT INT\n");

  if(fclose(f)==0) {
    fprintf(stdout,"  Template config file written to:\n%s\n",template);
    exit(status);
  } else {
    fprintf(stderr,
      "  [%s]: error: failed to write to: %s\n", __func__, template);
    exit(EXIT_FAILURE);
  }
}

/*
 * parse_config(file)
 *
 * Reads configuration option using file descriptor provided.
 */
void parse_config(FILE *file)
{
  fscanf(file, "%*26c %lu\n", &i_seed);
  fscanf(file, "%*26c %d\n", &i_edge_fcc_N);
  fscanf(file, "%*26c %d %d %d\n", &i_channel[0], &i_channel[1], &i_channel[2]);
  fscanf(file, "%*26c %lf\n",&i_channel_R);
  fscanf(file, "%*26c %d %d\n", &i_iDCfrom, &i_iDCto);

}

/*
 * parse_options(argc, argv)
 *
 * Process the optional arguments, using custom error messages.
 */
void parse_options(int argc, char *argv[], char**f)
{
  char optc;
  opterr = 0;
  while ((optc = getopt_long(argc, argv, "hvt", long_opts, &optind)) != -1) {
    switch (optc) {
      case 'h':
        print_usage(EXIT_SUCCESS);
        break;
      case 'v':
        print_version(EXIT_SUCCESS);
        break;
      case 't':
        generate_template_config(EXIT_SUCCESS);
        break;
      default:
        fprintf(stderr, "  [%s]: Unknown option '-%c'\n", prog_name, optopt);
        print_usage(EXIT_FAILURE);
        break;
    }
  }

  /*
   * After processing all options 'optind' contains the index of the next
   * command line argument.  Thus, if it is equal to the total number of
   * arguments, we are missing mandatory config file specification.  Complain
   * and exit in such a case, otherwise try to process it.
   */
  if(optind == argc) {
    fprintf(stderr, "  [%s]: error: missing config file specification\n",
            prog_name);
    print_usage(EXIT_FAILURE);
  }
  *f = argv[optind];
}

/* vim: set tw=80 ts=2 sw=2 et: */
