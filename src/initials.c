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
#include "io.h"

extern char *prog_name;

static void gen_template_confirmation(int ec, const char *desc,
  const char *file);

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
  Generate monomer, multimer or mixed strutures with optional inclusions.\n\
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
      -t, --template generate template config files\n\
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
  const char *fmt_version = "  %s version %s\n";
  const char *fmt_comit_id = "  git comit ID: %s\n";
  const char *fmt_comit_date = "  committed on: %s\n\n";
  const char *fmt_sha1 = "  %s  %s\n";
  const char *fmt_greet   = "  %-12s\t%s\n";
  const char *fmt_bullets = "  * %-34s: ";
  const char *fmt_category= "\n  %s:\n";


  #ifdef DATA_VISGL_OUTPUT
    const char *datavisGL = "Yes";
  #else
    const char *datavisGL = "No";
  #endif

  #ifdef PRNG_64BIT_MT19937
    const char *prng_type = "64bit MT19937";
  #endif

  #ifdef PRNG_32BIT_MT19937
    const char *prng_type = "32bit MT19937";
  #endif

  #ifdef PRNG_DRAND48
    const char *prng_type = "Drand48";
  #endif

  #ifdef VERBOSE_MDOE
    const char *verbose_mode = "Yes";
  #else
    const char *verbose_mode = "No";
  #endif

  #ifdef DEBUG_MODE
    const char *debug_mode = "Yes";
  #else
    const char *debug_mode = "No";
  #endif

  fprintf(stdout, "  X-mer Crystal Structure Generator (XCSG)\n");
  fprintf(stdout, fmt_version,  prog_name, code_version);
  fprintf(stdout, fmt_comit_id,  code_comit_id);
  fprintf(stdout, fmt_comit_date,  code_comit_date);
  print_greetings();

  fprintf(stdout, fmt_greet, "Created by:", author0);
  fprintf(stdout, fmt_greet, "Build by", builder);
  fprintf(stdout, fmt_greet, "Build date:", build);
  fprintf(stdout, fmt_greet, "Build host:", buildAt);

  fprintf(stdout, "  Try './%s -i' or '--info' for usage details\n", prog_name);

  fprintf(stdout, fmt_category, "Program build-in features");

  fprintf(stdout, fmt_bullets, "Data_visGL output files");
  fprintf(stdout, "%s\n", datavisGL);

  fprintf(stdout, fmt_bullets, "PRGN");
  fprintf(stdout, "%s\n", prng_type);

  fprintf(stdout, fmt_bullets, "Verbose mode");
  fprintf(stdout, "%s\n", verbose_mode);

  fprintf(stdout, fmt_bullets, "Enable debugging code");
  fprintf(stdout, "%s\n", debug_mode);

  fprintf(stdout, "\n");

  fprintf(stdout, "  Version control:\t  sources checksum's (SHA1)\n");
  fprintf(stdout, fmt_sha1, xcsg_c_SHA1,         "xcsg.c");
  fprintf(stdout, fmt_sha1, config_h_SHA1,      "config.h");
  fprintf(stdout, fmt_sha1, data_h_SHA1,        "data.h");
  fprintf(stdout, fmt_sha1, globals_h_SHA1,     "globals.h");
#ifdef PRNG_64BIT_MT19937
  fprintf(stdout, fmt_sha1, mt19937_64_h_SHA1,  "mt19937_64.h");
#endif
#ifdef PRNG_32BIT_MT19937
  fprintf(stdout, fmt_sha1, mt19937_h_SHA1,     "mt19937.h");
#endif
  fprintf(stdout, fmt_sha1, algebra_h_SHA1,     "algebra.h");
  fprintf(stdout, fmt_sha1, algebra_c_SHA1,     "algebra.c");
  fprintf(stdout, fmt_sha1, initials_h_SHA1,    "initials.h");
  fprintf(stdout, fmt_sha1, initials_c_SHA1,    "initials.c");
  fprintf(stdout, fmt_sha1, io_h_SHA1,          "io.h");
  fprintf(stdout, fmt_sha1, io_c_SHA1,          "io.c");
  fprintf(stdout, fmt_sha1, structure_h_SHA1,   "structure.h");
  fprintf(stdout, fmt_sha1, structure_c_SHA1,   "structure.c");
  fprintf(stdout, fmt_sha1, terminators_h_SHA1, "terminators.h");
  fprintf(stdout, fmt_sha1, terminators_c_SHA1, "terminators.c");
  fprintf(stdout, fmt_sha1, utils_h_SHA1,       "utils.h");
  fprintf(stdout, fmt_sha1, utils_c_SHA1,       "utils.c");
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
  extern const int length_max;
  FILE *f = NULL;
  const char *template = "%s.%s.template";
  char template_cfg[length_max], template_cha[length_max],
    template_sli[length_max];
  int fname_length[] = {
    snprintf(template_cfg, length_max -1, template, prog_name, "ini"),
    snprintf(template_cha, length_max -1, template, prog_name, "channels.dat"),
    snprintf(template_sli, length_max -1, template, prog_name, "slits.dat")};
  const int fname_qty = (int) sizeof(fname_length) / (int) sizeof(int);

  // Exit with error in case of array overflow and when binary name length
  // maches the generated template name (function would try to overwrite the
  // binary)
  for(int i=0; i<fname_qty; i++){
    if(fname_length[i] >= length_max || fname_length[i] == strlen(prog_name)+2){
      fprintf(stderr,
              " [%s] ERR: File name for config template (%d) is too short\n"
              " or the same as the binary. Shorten the name of the binary"
              " or increse the value\n of 'length_max' variable\n",
              __func__, i);
      exit(EXIT_FAILURE);
    }
  }
  // Opening file
  f = open_to_write(template_cfg);

  fprintf(f, "RNG seed                 : LUINT\n");
  fprintf(f, "Number of edge fcc cells : INT_x INT_y INT_z\n");
  fprintf(f, "Generated symmetry       : STRING\n");
  fprintf(f, "Structures (start end)   : INT INT\n");
  fprintf(f, "Make nano-channel (bool) : INT\n");
  fprintf(f, "Make nano-slit    (bool) : INT\n");
  fprintf(f, "Insert dimers DC  (bool) : INT\n");
  fprintf(f, "Number of channels       : INT\n");
  fprintf(f, "Channels desc. file name : STRING\n");
  fprintf(f, "Number of slits          : INT\n");
  fprintf(f, "Slits descrip. file name : STRING\n");

  gen_template_confirmation(fclose(f), "Main  program", template_cfg);

  // Opening file
  f = open_to_write(template_cha);

  fprintf(f, "#axis offset (3F); axis wersor (3F); channel radius (1F);");
  fprintf(f, " ch. sph. diameter (1F). The 1st line is skipped.\n");
  fprintf(f, " DOUBLE_ox DOUBLE_oy DOUBLE_oz ");
  fprintf(f, " DOUBLE_cx DOUBLE_cy DOUBLE_cz ");
  fprintf(f, " DOUBLE_cr DOUBLE_sd\n");

  gen_template_confirmation(fclose(f), "Channel desc.", template_cha);

  // Opening file
  f = open_to_write(template_sli);

  fprintf(f, "#plane offset (3F); plane normal (3F); pl. thickness (1F);");
  fprintf(f, " plane sph. diameter (1F). 1st line is skipped.\n");
  fprintf(f, " DOUBLE_ox DOUBLE_oy DOUBLE_oz ");
  fprintf(f, " DOUBLE_cx DOUBLE_cy DOUBLE_cz ");
  fprintf(f, " DOUBLE_cr DOUBLE_sd\n");

  gen_template_confirmation(fclose(f), "Slit descrip.", template_sli);
  exit(status);
}

static void gen_template_confirmation(int ec, const char *desc,
  const char *file){
  extern const char *fmt_writting_failed;
  if(ec == 0){
    fprintf(stdout," %s template config file written to: %s\n", desc, file);
  }else{
    fprintf(stderr, fmt_writting_failed, __func__, file);
    exit(EXIT_FAILURE);
  }
}

/*
 * parse_channels(file, ch_tab)
 *
 * Reads channels' configuration from 'file' and stores data for each channel
 * respectively.
 */
int parse_channels(FILE *file, CHA ch_tab[], int num_channels)
{
  extern const char *fmt_IO_8f;
  extern const char *fmt_sudden_eof;
  extern const char *fmt_null_ptr;
  int i=0;
  double c_off[3] = {0e0, 0e0, 0e0};
  double c_nor[3] = {0e0, 0e0, 0e0};
  double c_r = 0e0;
  double s_d = 0e0;
  unsigned char trash;

  if(file != NULL){
    // Skip the header line in file
    // https://stackoverflow.com/questions/2799612/how-to-skip-the-first-line-when-fscanning-a-txt-file
    do{
      trash = fgetc(file);
    }while (trash != '\n');

    // Read data from the file
    for(i=0; i<num_channels; i++){
      if(fscanf(file, fmt_IO_8f,
            &c_off[0], &c_off[1], &c_off[2],
            &c_nor[0], &c_nor[1], &c_nor[2], &c_r, &s_d) == EOF){
        fprintf(stderr, fmt_sudden_eof, __func__, i, num_channels);
        fclose(file);
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
    fclose(file);
    return EXIT_SUCCESS;
  }else{
    fprintf(stderr, fmt_null_ptr, __func__);
    return EXIT_FAILURE;
  }
}

/*
 * parse_slits(file, sl_tab)
 *
 * Reads slits' configuration from 'file' and stores data for each slit
 * respectively.
 */
int parse_slits(FILE *file, SLI sl_tab[], int num_slits)
{
  extern const char *fmt_IO_8f;
  extern const char *fmt_sudden_eof;
  extern const char *fmt_null_ptr;
  int i=0;
  double sl_off[3] = {0e0, 0e0, 0e0};
  double sl_nor[3] = {0e0, 0e0, 0e0};
  double sl_th = 0e0;
  double s_d = 0e0;
  unsigned char trash;

  if(file != NULL){
    // Skip the header line in file
    // https://stackoverflow.com/questions/2799612/how-to-skip-the-first-line-when-fscanning-a-txt-file
    do{
      trash = fgetc(file);
    }while (trash != '\n');

    // Read data from the file
    for(i=0; i<num_slits; i++){
      if(fscanf(file, fmt_IO_8f,
            &sl_off[0], &sl_off[1], &sl_off[2],
            &sl_nor[0], &sl_nor[1], &sl_nor[2], &sl_th, &s_d) == EOF){
        fprintf(stderr, fmt_sudden_eof, __func__, i, num_slits);
        fclose(file);
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
    fclose(file);
    return EXIT_SUCCESS;
  }else{
    fprintf(stderr, fmt_null_ptr, __func__);
    return EXIT_FAILURE;
  }
}

/*
 * parse_config(file)
 *
 * Reads configuration option using file descriptor provided.
 */
int parse_config(FILE *file, CONFIG *cfg)
{
  extern const char *fmt_null_ptr;
  const char *fmt_bad_range_values =
    "  [%s] ERR: Incorrect parameters for structure indexes\n";
  const char *fmt_lu  = "%*26c %lu\n";
  const char *fmt_ddd = "%*26c %d %d %d\n";
  const char *fmt_dd  = "%*26c %d %d\n";
  const char *fmt_d   = "%*26c %d\n";
  const char *fmt_s   = "%*26c %s\n";
  int exit_code = EXIT_SUCCESS;

  if(file != NULL){
    fscanf(file, fmt_lu, &cfg->seed);
    fscanf(file, fmt_ddd,
          &cfg->cells[0], &cfg->cells[1], &cfg->cells[2]);
    fscanf(file, fmt_s,   cfg->symmetry);
    fscanf(file, fmt_dd, &cfg->first, &cfg->last);
    fscanf(file, fmt_d,  &cfg->mk_channel);
    fscanf(file, fmt_d,  &cfg->mk_slit);
    fscanf(file, fmt_d,  &cfg->mk_dimers);
    fscanf(file, fmt_d,  &cfg->num_channels);
    fscanf(file, fmt_s,   cfg->cfg_channels);
    fscanf(file, fmt_d,  &cfg->num_slits);
    fscanf(file, fmt_s,   cfg->cfg_slits);

    // Parameters sanity check
    if(cfg->last < cfg->first || cfg->first < 0 || cfg->last < 0){
      fprintf(stderr, fmt_bad_range_values, __func__);
      exit_code = EXIT_FAILURE;
    }

    exit_code = config_inclusions_sanity_check(
      cfg->mk_channel, cfg->num_channels, "num_channels");

    exit_code = config_inclusions_sanity_check(
      cfg->mk_slit, cfg->num_slits, "num_slits");

    fclose(file);
  }else{
    fprintf(stderr, fmt_null_ptr, __func__);
    return EXIT_FAILURE;
  }
  return exit_code;
}

/*
 * config_inclusions_sanity_check()
 *
 * Prints error messages for wrong inlcusion parameters
 */
int config_inclusions_sanity_check(
  int mk_inclusion, int num_inclusions, const char *msg)
{
  const char *fmt_switch_sanity_check =
    "  [%s] ERR: %s must be greater than 0\n";
  if(mk_inclusion != 0 && num_inclusions < 1){
    fprintf(stderr, fmt_switch_sanity_check, msg, __func__);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/*
 * parse_options(argc, argv)
 *
 * Process the optional arguments, using custom error messages.
 */
void parse_options(int argc, char *argv[], char**f)
{
  char optc;
  const char *fmt_missing_config =
    "  [%s] ERR: missing config file specification\n";
  const char *fmt_unknown_option = "  [%s]: Unknown option '-%c'\n";

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
        fprintf(stderr, fmt_unknown_option, prog_name, optopt);
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
    fprintf(stderr, fmt_missing_config, prog_name);
    print_help(EXIT_FAILURE);
  }
  *f = argv[optind];
}

/* vim: set tw=80 ts=2 sw=2 et: */
