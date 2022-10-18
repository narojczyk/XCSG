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
#include "utils.h"

extern char *prog_name;
extern const char* valid_config_version;

extern const double zero;

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
  const char *fmt_usage_info =
"%s is controlled by one main config file  and two more (optional) files  with\n"
"the description of inclusions. Running the program with -t option will generate\n"
"all possible config file templates. Details of the setup options:\n"
"\n\tMain configuration file\n"
"* line 1: version compliance signature (current sig. written by default),\n"
"* line 2: random number generator seed (long integer),\n"
"* line 3: number of unit cells in each direction (three integers),\n"
"* line 4: symmetry of the output structure (string) - only 'fcc' at the  moment,\n"
"* line 5: number of the first and last of generated  structures (two  integers),\n"
"* line 6,7: channel and  layer  inclusion, respectively.  Integer string integer\n"
"  [0|1] flag  to insert  inclusion,  name  of  the file  with inclusion details,\n"
"  integer value with the number of inclusions in said file,\n"
"* line 8: Integer  [0|1] flag  whether to join  particles  outside the inclusion\n"
"  into dimers and arrange them into a disordered (degenerate crystal) phase,\n"
"* line 9: 'Rough'  inclusions - concerns  inclusions  inserted into DC phase and\n"
"  filled with molecules.  This flag [0|1] has no effect in  monomer systems, but\n"
"  in molecular systems, it toggles whether the molecules that cross the boundary\n"
"  of the inclusion are broken (0) or added to inclusion as a whole (1).\n"
"Values must be entered after the colon ':', first 26 characters  in each line is\n"
"omitted.\n"
"\n\tInclusion configuration files\n"
"The first line is a header and is omited. The generated template contains header\n"
"line describing each field(s). They are: inclusion offset (three floats), inclu-\n"
"sion orientation vector (three floats), inclusion size (one float),  diameter of\n"
"spheres forming this inclusion, and an integer value indicating  targeted parti-\n"
"cles, that this inclusion should be  filled with (e.g. 1 - spheres, 2 - dimers).\n"
"By default, an inclusion is pinned to the atom with the smallest coordinates (at\n"
"one of the corners of the system). Inclusion offset is a vector  translating the\n"
"point of the axis (channel) or the position of the  layer in  the  crystal.  The\n"
"orientation vector is  the vector along  the axis (channel) or the normal to the\n"
"plane (layer).\n"
"\n\tOutput files\n"
"The program outputs structure in at least two files:\n"
"\ta) s3d0_summary_XXXXX.ini\n"
"\tb) s3d1_monomers_XXXXX.csv\n"
"\tc) s3d2_dimers_XXXXX.csv (if and dimers are present)\n"
"In a) the statistics of different types of particles and parameters  of the  box\n"
"matrix are given. In b) spheres data is saved. The columns  represent  (resp.):\n"
"sphere id, sphere type (depending if it forms a molecule or belongs to an inclu-\n"
"sion), sphere  position (three  floats), sphere diameter, twelve integer indica-\n"
"ting id of the nearest neighbours,and three integers being indices of the sphere\n"
"position in the crystalline  structure.  The positions  correspond to  the close\n"
"packing of monodisperse spheres of the diameter equal to one. In c)  the columns\n"
"represent: id of the dimer, id of its type, position  of  the  centre of  dimers\n"
"mass (three floats), dimer orientation (unit vector - three floats),  and  dimer\n"
"length along with the two integer id's of spheres that create the  dimer.  Other\n"
"files that contain structure data in formats that can be easily fed to graphical\n"
"programs are also generated. At the moment this includes PoV-RAY scripts as well\n"
"as the not-generally available OpenGL viewer used by the authors of this code.\n";


  fprintf(stdout, "\tUsing %s\n\n", prog_name);

  fprintf(stdout, fmt_usage_info, prog_name);

  exit(status);
}

/*
 * print_version()
 *
 * Display build info and checksums of sources used
 */
void print_version(int status)
{
  const char *fmt_version = "\n  %s version %s\n";
  const char *fmt_commit_id = "  git commit ID: %s on branch: %s\n";
  const char *fmt_commit_date = "  committed on: %s\n";
  const char *fmt_cc_vendor = "  build with: %s\n";
  const char *fmt_cc_version = "  CC version: %s\n";
  const char *fmt_cc_flags = "  CC flags: %s\n\n";
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
  fprintf(stdout, fmt_version,     prog_name, code_version);
  fprintf(stdout, fmt_commit_id,   code_commit_id, code_branch_name);
  fprintf(stdout, fmt_commit_date, code_commit_date);
  fprintf(stdout, fmt_cc_vendor,   cc_vendor);
  fprintf(stdout, fmt_cc_version,  cc_version);
  fprintf(stdout, fmt_cc_flags,    cc_flags);

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
  fprintf(stdout, fmt_sha1, xcsg_c_SHA1,        "xcsg.c");
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

  fprintf(f, "Version compliance sig.  : %s\n", valid_config_version);
  fprintf(f, "RNG seed                 : LUINT\n");
  fprintf(f, "Number of edge fcc cells : INT_x INT_y INT_z\n");
  fprintf(f, "Generated symmetry       : STRING\n");
  fprintf(f, "Structures (start end)   : INT INT\n");
  fprintf(f, "Channels (file name)     : STRING\n");
  fprintf(f, "Slits    (file name)     : STRING\n");
  fprintf(f, "Rough inclusions (bool)  : INT\n");
  fprintf(f, "Insert dimers DC (bool)  : INT\n");
  fprintf(f, "Dimer length at sigma    : INT\n");

  gen_template_confirmation(fclose(f), "Main  program", template_cfg);

  // Opening file
  f = open_to_write(template_cha);

  fprintf(f, "#axis offset (3F); axis wersor (3F); channel radius (1F);");
  fprintf(f, " ch. sph. diameter (1F); n-mer to form (1I).");
  fprintf(f, " The 1st line is skipped.\n");
  fprintf(f, " DOUBLE_ox DOUBLE_oy DOUBLE_oz ");
  fprintf(f, " DOUBLE_cx DOUBLE_cy DOUBLE_cz ");
  fprintf(f, " DOUBLE_cr DOUBLE_sd INT_n-mer\n");

  gen_template_confirmation(fclose(f), "Channel desc.", template_cha);

  // Opening file
  f = open_to_write(template_sli);

  fprintf(f, "#plane offset (3F); plane normal (3F); pl. thickness (1F);");
  fprintf(f, " plane sph. diameter (1F); n-mer to form (1I).");
  fprintf(f, " The 1st line is skipped.\n");
  fprintf(f, " DOUBLE_ox DOUBLE_oy DOUBLE_oz ");
  fprintf(f, " DOUBLE_cx DOUBLE_cy DOUBLE_cz ");
  fprintf(f, " DOUBLE_cr DOUBLE_sd INT_n-mer\n");

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
 * parse_inclusions(file, inc)
 *
 * Reads inclusions' configuration from 'file' and stores data for each array
 * element of 'inc'.
 */
INC* parse_inclusions(FILE *file)
{
  INC template = template_inclusion();
  INC *inclusions_array = NULL;
  INC *list_of_inclusions = NULL;
  INC *current = NULL, *last = NULL;
  extern const char *fmt_null_ptr;
  const char *fmt_IO_8f1d = "%lf %lf %lf %lf %lf %lf %lf %lf %d\n";
  int i=0, incl_counter=0, n_mer;
  double off[3] = {zero, zero, zero};
  double nor[3] = {zero, zero, zero};
  double size = zero;
  double sd = zero;
  unsigned char trash;

  if(file != NULL){
    // Skip the header line in file
    // https://stackoverflow.com/questions/2799612/how-to-skip-the-first-line-when-fscanning-a-txt-file
    do{
      trash = fgetc(file);
    }while (trash != '\n');

    while(fscanf(file, fmt_IO_8f1d, &off[0], &off[1], &off[2],
      &nor[0], &nor[1], &nor[2], &size, &sd, &n_mer) != EOF){

        current = malloc(sizeof(INC));
        (*current) = template;

        // Store offset for the inclusion
        (*current).os[0] = off[0];
        (*current).os[1] = off[1];
        (*current).os[2] = off[2];

        // Store normal vector for the inclusion
        (*current).nm[0] = nor[0];
        (*current).nm[1] = nor[1];
        (*current).nm[2] = nor[2];

        // Store size of the inclusion for all the possible cases
        (*current).thickness = size;
        (*current).radius = size;

        // Store the diameter of spheres in inclusion
        (*current).sph_d = sd;

        // Store the target n-mer size to be formed
        (*current).tgt_Nmer = n_mer;

        if(list_of_inclusions != NULL){
          // On the consecutive iterations update the 'next' field of
          // the previous entry and advance the end-pointer to the current entry
          (*last).next = current;
          last = current;
        }else{
          // On first iteration initiate the list of inclusions to
          // the current, newly generated entry
          list_of_inclusions = current;
          last = current;
        }
        // Count the number of list elements
        incl_counter++;
      }

      // If any inclusions loaded from file
      if(incl_counter){
        // allocate an array of the corresponding size
        inclusions_array = malloc(incl_counter * sizeof(INC));
        // point at the first inclusion on the linked list
        current = list_of_inclusions;
        // For all elelemts on the linked list
        while(current != NULL){
          // copy the inclusion into the array
          inclusions_array[i] = (*current);
          inclusions_array[i].next = &(inclusions_array[i+1]);
          i++;
          // advance the pointer to the next list element and free
          // the memory for the old one
          last = current;
          current = (*current).next;
          free(last);
        }
        // NULL-out the 'next'-pointer of the last array element
        inclusions_array[incl_counter - 1].next = NULL;
      }

      // Will return NULL if no inclusions had been read
      return inclusions_array;
  }

  // Complain of no vaild file pointer received on input
  fprintf(stderr, fmt_null_ptr, __func__);
  return NULL;
}

/*
 * parse_config(file)
 *
 * Reads configuration option using file descriptor provided.
 */
int parse_config(FILE *file, CONFIG *cfg)
{
  extern const char *fmt_null_ptr;
  extern const int config_version_length;
  extern const char* older_config_version;
  const char *fmt_bad_range_values =
    "  [%s] ERR: Incorrect parameters for structure indexes\n";
  const char *fmt_unknown_config_format =
    "  [%s] ERR: Config version not recognized.\n";
//   const char *fmt_legacy_config =
//     "  [%s] WARN: Config version not recognized."
//     " Attempting to parse in legacy format\n";
  const char *fmt_lu  = "%*26c %lu\n";
  const char *fmt_ddd = "%*26c %d %d %d\n";
  const char *fmt_dsd = "%*26c %d %s %d\n";
  const char *fmt_dd  = "%*26c %d %d\n";
  const char *fmt_d   = "%*26c %d\n";
  const char *fmt_s   = "%*26c %s\n";
  char cversion[config_version_length];
  int exit_code = EXIT_SUCCESS;

  if(file != NULL){
    // Look for version marker in a config file
    fscanf(file, fmt_s, cversion);

    // If version agrees, read parameters in current format
    if(!str_validate(cversion, valid_config_version)){

      fscanf(file, fmt_lu,  &cfg->seed);
      fscanf(file, fmt_ddd, &cfg->cells[0], &cfg->cells[1], &cfg->cells[2]);
      fscanf(file, fmt_s,    cfg->symmetry);
      fscanf(file, fmt_dd,  &cfg->first, &cfg->last);
      fscanf(file, fmt_s,    cfg->cfg_channels);
      fscanf(file, fmt_s,    cfg->cfg_slits);
      fscanf(file, fmt_d,   &cfg->rough_inclusions);
      fscanf(file, fmt_d,   &cfg->mk_dimers);
      fscanf(file, fmt_d,   &cfg->dimer_L_sigma);

      fclose(file);
    }else if (!str_validate(cversion, older_config_version)){

      fscanf(file, fmt_lu,  &cfg->seed);
      fscanf(file, fmt_ddd, &cfg->cells[0], &cfg->cells[1], &cfg->cells[2]);
      fscanf(file, fmt_s,    cfg->symmetry);
      fscanf(file, fmt_dd,  &cfg->first, &cfg->last);
      fscanf(file, fmt_dsd, &cfg->mk_channel, cfg->cfg_channels,
             &cfg->num_channels);
      fscanf(file, fmt_dsd, &cfg->mk_slit, cfg->cfg_slits, &cfg->num_slits);
      fscanf(file, fmt_d,   &cfg->mk_dimers);
      fscanf(file, fmt_d,   &cfg->rough_inclusions);

      fclose(file);
    }else{
      // Echo unknown config format
      fprintf(stderr, fmt_unknown_config_format, __func__);
      return EXIT_FAILURE;
    }
  }else{
    fprintf(stderr, fmt_null_ptr, __func__);
    return EXIT_FAILURE;
  }

  // Parameters sanity check
  // * main ini file (TODO: Extend and move this to a function)
  if(cfg->last < cfg->first || cfg->first < 0 || cfg->last < 0){
    fprintf(stderr, fmt_bad_range_values, __func__);
    exit_code = EXIT_FAILURE;
  }

  return exit_code;
}

/*
 * config_inclusions_sanity_check()
 *
 * Prints error messages for wrong inlcusion parameters
 */
int inclusion_count_sanity_check(int num_inclusions, const char *msg)
{
  const char *fmt_switch_sanity_check =
    "  [%s] ERR: %s must be greater than 0\n";
  if(num_inclusions < 1){
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
