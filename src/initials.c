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
#include "data.h"
#include "checksum.h"
#include "initials.h"

extern char *prog_name;

extern unsigned long int i_seed;
extern int i_edge_fcc_N[3];
extern int i_ch_layout[3];
extern int i_iDCfrom;
extern int i_iDCto;
extern int i_make_channel;
extern int i_make_slit;
extern int i_fs_connect;
extern int i_n_channels; 
extern int i_n_slits; 
extern char i_Fchannels[41];
extern char i_Fslits[41];

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
  Generate a D.C. phase of dimers (with inclusions) on the f.c.c. lattice\n\
  Institute of Molecular Physics, Polish Academy of Sciences\n\n", stdout);
}

/*
 * print_help()
 *
 * Display description of command line arguments
 */
void print_help(int status)
{
  if(status != EXIT_SUCCESS){
    fprintf(stderr, "  Try '%s --help' for more information.\n", prog_name);
  }else{
    fprintf(stdout, "  Usage: %s [OPTIONS] [CONFIG_FILE]\n\n", prog_name);
    fputs("\
    Option list:\n\
      -h, --help     print this information\n\
      -i, --info     print detailed usage description\n\
      -t, --template generate template config file\n\
      -v, --version  print version and build details\n", stdout);
  }
  exit(status);
}

/*
 * print_info()
 *
 * Display information about the program usage
 */
void print_info(int status)
{
  fprintf(stdout, "  %s usage description\n\n", prog_name);
  
  fprintf(stdout, "  will be provided here eventually :)\n");
  
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

//   fprintf(stdout, "  Created by:\t%s\n\t\t%s\n", author0, author1);
  fprintf(stdout, "  Created by:\t%s\n", author0);
  fprintf(stdout, "  Build by:\t%s\n", builder);
  fprintf(stdout, "  Build date:\t%s\n", build);
  fprintf(stdout, "  Build host:\t%s\n\n", buildAt);
  
  fprintf(stdout, "  Try '-i' or '--info' for usage details\n\n");

  fprintf(stdout, "  Program build-in features:\n");
    
  fprintf(stdout, "  * verbose output for debugging      : ");
  #ifdef DEBUG_MODE
    fprintf(stdout, "Yes\n");
  #else
    fprintf(stdout, "No\n");
  #endif
  
    fprintf(stdout, "  * use 64bit MT19937 random num. gen.: ");
  #ifdef USE_64BIT_MT19937
    fprintf(stdout, "Yes\n");
  #else
    fprintf(stdout, "No\n");
  #endif
    
    fprintf(stdout, "  * use 32bit MT19937 random num. gen.: ");
  #ifdef USE_32BIT_MT19937
    fprintf(stdout, "Yes\n");
  #else
    fprintf(stdout, "No\n");
  #endif
    
    fprintf(stdout, "  * use drand48 random num. gen.      : ");
  #ifdef USE_DRAND48
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
  #endif    
  #ifdef USE_32BIT_MT19937
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
  fprintf(stdout, "  %s  %s\n", structure_h_SHA1, "structure.h");
  fprintf(stdout, "  %s  %s\n", structure_c_SHA1, "structure.c");
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
  char *template_ch = "channels.dat.sed";
  char *template_sl = "slits.dat.sed";

  // Opening file
  if ((f = fopen(template,"w")) == NULL ) {
    fprintf(stderr,
      "  [%s]: error: failed to open file %s\n", __func__, template);
    exit(1);
  }

  fprintf(f, "RNG seed                : LUINT\n");
  fprintf(f, "Number of edge fcc cells: INT_x INT_y INT_z\n");
//   fprintf(f, "Slit normal vector      : INT_h INT_k INT_l\n");
//   fprintf(f, "Chanel radius [sigma]   : DOUBLE\n");
  fprintf(f, "Load initial DC struct. : INT INT\n");
  fprintf(f, "Make nano-channel (bool): INT\n");
  fprintf(f, "Make nano-slit    (bool): INT\n");
//   fprintf(f, "Slit thickness [sigma]  : DOUBLE\n");
  fprintf(f, "Free sph. connect (bool): INT\n");
  fprintf(f, "Number of channels      : INT\n");
  fprintf(f, "Channels desc. file name: STRING\n");
  fprintf(f, "Number of slits         : INT\n");
  fprintf(f, "Slits descrip. file name: STRING\n");
  
  if(fclose(f)==0) {
    fprintf(stdout," Main  program template config file written to: %s\n",
            template);
  } else {
    fprintf(stderr,
      "  [%s]: error: failed to write to: %s\n", __func__, template);
    exit(EXIT_FAILURE);
  }
  
  // Opening file
  if ((f = fopen(template_ch,"w")) == NULL ) {
    fprintf(stderr,
      "  [%s]: error: failed to open file %s\n", __func__, template_ch);
    exit(1);
  }
  
  fprintf(f, "#axis offset (3F); axis wersor (3F); channel radius (1F);");
  fprintf(f, " ch. sph. diameter (1F). The 1st line is skipped.\n");
  fprintf(f, " DOUBLE_ox DOUBLE_oy DOUBLE_oz ");
  fprintf(f, " DOUBLE_cx DOUBLE_cy DOUBLE_cz ");
  fprintf(f, " DOUBLE_cr DOUBLE_sd\n");
  
  if(fclose(f)==0) {
    fprintf(stdout," Channel desc. template config file written to: %s\n",
            template_ch);
  } else {
    fprintf(stderr,
      "  [%s]: error: failed to write to: %s\n", __func__, template_ch);
    exit(EXIT_FAILURE);
  }
  
  // Opening file
  if ((f = fopen(template_sl,"w")) == NULL ) {
    fprintf(stderr,
      "  [%s]: error: failed to open file %s\n", __func__, template_sl);
    exit(1);
  }
                
  fprintf(f, "#plane offset (3F); plane normal (3F); pl. thickness (1F);");
  fprintf(f, " plane sph. diameter (1F). 1st line is skipped.\n");
  fprintf(f, " DOUBLE_ox DOUBLE_oy DOUBLE_oz ");
  fprintf(f, " DOUBLE_cx DOUBLE_cy DOUBLE_cz ");
  fprintf(f, " DOUBLE_cr DOUBLE_sd\n");
  
  if(fclose(f)==0) {
    fprintf(stdout," Slit descrip. template config file written to: %s\n",
            template_sl);
    exit(status);
  } else {
    fprintf(stderr,
      "  [%s]: error: failed to write to: %s\n", __func__, template_sl);
    exit(EXIT_FAILURE);
  }
  
}

/*
 * parse_channels(file, ch_tab)
 * 
 * Reads channels' configuration from 'file' and stores data for each channel
 * respectively.
 */
int parse_channels(FILE *file, CHA ch_tab[])
{
  int i=0;
  double c_off[3] = {0e0, 0e0, 0e0};
  double c_nor[3] = {0e0, 0e0, 0e0};
  double c_r = 0e0;
  double s_d = 0e0;
  
  // Skip first line in the file
  fscanf(file, "%*[^\n]\n", NULL);

  // Read data from the file
  for(i=0; i<i_n_channels; i++){
    if(fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf\n", 
           &c_off[0], &c_off[1], &c_off[2],
           &c_nor[0], &c_nor[1], &c_nor[2], &c_r, &s_d) == EOF){
      fprintf(stderr,"  [%s]: faile ended unexpectedly\n", __func__);
      fprintf(stderr,"  %d data lines read; %d data lines expected\n", 
              i, i_n_channels);
      return EXIT_FAILURE;      
    }else{
    
    // Store offset for the channel 'i'
    ch_tab[i].os[0] = c_off[0];
    ch_tab[i].os[1] = c_off[1];
    ch_tab[i].os[2] = c_off[2];
    
    // Store normal vector for the channel 'i'
    ch_tab[i].nm[0] = c_nor[0];
    ch_tab[i].nm[1] = c_nor[1];
    ch_tab[i].nm[2] = c_nor[2];
    
    // Store radius for the channel 'i'
    ch_tab[i].radius = c_r;
    
    // Store diameter of spheres in channel 'i'
    ch_tab[i].sph_d = s_d;
    }
  }
  return EXIT_SUCCESS;
}

/*
 * parse_slits(file, sl_tab)
 * 
 * Reads slits' configuration from 'file' and stores data for each slit
 * respectively.
 */
int parse_slits(FILE *file, SLI sl_tab[])
{
  int i=0;
  double sl_off[3] = {0e0, 0e0, 0e0};
  double sl_nor[3] = {0e0, 0e0, 0e0};
  double sl_th = 0e0;
  double s_d = 0e0;
  
  // Skip first line in the file
  fscanf(file, "%*[^\n]\n", NULL);

  // Read data from the file
  for(i=0; i<i_n_slits; i++){
    if(fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf\n", 
           &sl_off[0], &sl_off[1], &sl_off[2],
           &sl_nor[0], &sl_nor[1], &sl_nor[2], &sl_th, &s_d) == EOF){
      fprintf(stderr,"  [%s]: faile ended unexpectedly\n", __func__);
      fprintf(stderr,"  %d data lines read; %d data lines expected\n", 
              i, i_n_slits);
      return EXIT_FAILURE;      
    }else{
    
    // Store offset for the slit 'i'
    sl_tab[i].os[0] = sl_off[0];
    sl_tab[i].os[1] = sl_off[1];
    sl_tab[i].os[2] = sl_off[2];
    
    // Store normal vector for the slit 'i'
    sl_tab[i].nm[0] = sl_nor[0];
    sl_tab[i].nm[1] = sl_nor[1];
    sl_tab[i].nm[2] = sl_nor[2];
    
    // Store radius for the slit 'i'
    sl_tab[i].thickness = sl_th;
    
    // Store diameter of spheres in slit 'i'
    sl_tab[i].sph_d = s_d;
    }
  }
  return EXIT_SUCCESS;
}

/*
 * parse_config(file)
 *
 * Reads configuration option using file descriptor provided.
 */
void parse_config(FILE *file)
{
  fscanf(file, "%*26c %lu\n",      &i_seed);
  fscanf(file, "%*26c %d %d %d\n", &i_edge_fcc_N[0], &i_edge_fcc_N[1], 
         &i_edge_fcc_N[2]);
  fscanf(file, "%*26c %d %d\n",    &i_iDCfrom, &i_iDCto);
  fscanf(file, "%*26c %d\n",       &i_make_channel);
  fscanf(file, "%*26c %d\n",       &i_make_slit);
  fscanf(file, "%*26c %d\n",       &i_fs_connect);
  fscanf(file, "%*26c %d\n",       &i_n_channels);
  fscanf(file, "%*26c %s\n",        i_Fchannels);
  fscanf(file, "%*26c %d\n",       &i_n_slits);
  fscanf(file, "%*26c %s\n",        i_Fslits);
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
  while ((optc = getopt_long(argc, argv, "hvit", long_opts, &optind)) != -1) {
    switch (optc) {
      case 'h':
        print_help(EXIT_SUCCESS);
        break;
      case 'i':
        print_info(EXIT_SUCCESS);
        break;        
      case 't':
        generate_template_config(EXIT_SUCCESS);
        break;        
      case 'v':
        print_version(EXIT_SUCCESS);
        break;
      default:
        fprintf(stderr, "  [%s]: Unknown option '-%c'\n", prog_name, optopt);
        print_help(EXIT_FAILURE);
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
    print_help(EXIT_FAILURE);
  }
  *f = argv[optind];
}

/* vim: set tw=80 ts=2 sw=2 et: */
