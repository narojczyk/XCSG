/*
 * x-mer Crystal Structure Generator (XCSG)
 * Author: Jakub Narojczyk <narojczyk@ifmpan.poznan.pl>
 * (c) 2015-2021
 *
 * The program generates a crystal structure of f.c.c. spheres. The latter can
 * be connected into multimer molecules. The structure can also include
 * additional modifications like nanochannels, nanolayers or point inclusions.
 *
 */

/*
 * TODO:
 *
 */

#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "data.h"
#include "globals.h"
#include "initials.h"
#include "io.h"
#include "utils.h"
#include "structure.h"

int main(int argc, char *argv[])
{
  FILE *f;
  CONFIG cfg;
  MODEL mdl;

  SPH *spheres = NULL;
  DIM3D *dimers = NULL;
  CHA *channels = NULL;
  SLI *slits = NULL;

  char *f_ini = NULL;

  int exit_status = EXIT_SUCCESS; // Initial assumption.
  int Nd2, Ns3, Ns3_odd;
  int zip_Ns3_runs=0, zip_init_sph=-1;
  int Odistrib[6] = {0,0,0,0,0,0};
  int valid_dimer_pair[2];
  int i, s;
  int low_on_dimers = 0;
  int flip_count = 0;
  int fsi, fsn, fsi_chances;

  const char *fmt_open_config_failed =
    " [%s] ERR: cannot open config file: %s\n";
  const char *fmt_open_inclusion_failed =
    " [%s] ERR: cannot open inclusion file: %s\n";
  const char *fmt_missing_inclusion_normal =
    " [%s] WRN: missing normal vector for %s inc. no %d, skipping\n";
  const char *fmt_dimer_distr_header =
    " %-28s %5s %5s %5s %5s %5s %5s\n";
  const char *fmt_wrong_number_of_spheres =
    " [%s] ERR: failed to calculate the number of spheres\n";
  const char *fmt_wrong_container_dimensions =
    " [%s] ERR: failed to calculate continer dimensions\n";

  // Extract program name from the path.
  prog_name = basename(argv[0]);

  // Parse command line options and get the config file name
  parse_options(argc, argv, &f_ini);

  // Open and parse config file
  if((f = fopen(f_ini, "r")) == NULL) {
    fprintf(stderr, fmt_open_config_failed, prog_name, f_ini);
    return EXIT_FAILURE;
  }
  exit_status = parse_config(f, &cfg);
  fclose(f);
  if(exit_status != EXIT_SUCCESS){
      goto cleanup;
  }

  // Allocate and clean memory for channels' data
  channels = malloc(cfg.num_channels * sizeof(CHA));
  memory_clean_channels(channels, cfg.num_channels);

  // Allocate and clean memory for slits' data
  slits = malloc(cfg.num_slits * sizeof(SLI));
  memory_clean_slits(slits, cfg.num_slits);

  // Open and read channel description data
  if( cfg.mk_channel != 0 ){
    if((f = fopen(cfg.cfg_channels, "r")) == NULL) {
      fprintf(stderr, fmt_open_inclusion_failed, prog_name, cfg.cfg_channels);
      return EXIT_FAILURE;
    }
    exit_status = parse_channels(f, channels, cfg.num_channels);
    fclose(f);
    if(exit_status != EXIT_SUCCESS){
      goto cleanup;
    }
  }

  // Open and read slits description data
  if( cfg.mk_slit != 0 ){
    if((f = fopen(cfg.cfg_slits, "r")) == NULL) {
      fprintf(stderr, fmt_open_inclusion_failed, prog_name, cfg.cfg_slits);
      return EXIT_FAILURE;
    }
    exit_status = parse_slits(f, slits, cfg.num_slits);
    fclose(f);
    if(exit_status != EXIT_SUCCESS){
      goto cleanup;
    }
  }

  // Initiate generator with 'seed'
  init_RNG(cfg.seed);

  // Set the number of monomers, dimers, and possibly, further-mers
  mdl.Nsph = number_of_spheres(cfg);
  if(mdl.Nsph == EXIT_FAILURE){
    fprintf(stderr, fmt_wrong_number_of_spheres, prog_name);
    goto cleanup;
  }
  mdl.Ndim = mdl.Nsph / 2;

  // Calculate container dimensions
  if(container_dimensions(&mdl, cfg) == EXIT_FAILURE){
    fprintf(stderr, fmt_wrong_container_dimensions, prog_name);
    goto cleanup;
  }

  // Summary of configuration variables
  display_configuration_summary(cfg, mdl, slits, channels);

  // Allocate memory for spheres and dimers
  spheres = malloc( mdl.Nsph * sizeof(SPH));
  dimers  = malloc( mdl.Ndim * sizeof(DIM3D));

  // Loop over selected set of structures
  for(s=cfg.first; s<=cfg.last; s++){

    // Clean allocated memory
    memory_clean_spheres(spheres, mdl.Nsph);
    memory_clean_dimers(dimers, mdl.Ndim);

    fprintf(stdout, "\n ***\tProcessing structure %d\n",s);
    fprintf(stdout, " Generating pure f.c.c. structure\n");

    // Set fcc structure of spheres
    sph_set_fcc(spheres, mdl.Nsph, cfg.cells);

    // Find neighbors for spheres
    if(find_ngb_spheres(spheres, mdl.Nsph, mdl.box) != 0){
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }

    // Assign lattice indexes to all spheres
    exit_status = sph_assign_lattice_indexes(spheres, mdl.Nsph);
    if(exit_status != EXIT_SUCCESS){
      goto cleanup;
    }

    // Flag all dimers as broken at this point
    for(i=0; i<mdl.Ndim; i++){
      dimers[i].type = 2;
    }

    // Make channel
    if(cfg.mk_channel){
      for(i=0; i<cfg.num_channels; i++){
        // Test channel data
        if(channels[i].nm[0] != 0 || channels[i].nm[1] != 0 || channels[i].nm[2] != 0){
          fprintf(stdout, " Inserting channel %d\n", i);
          make_channel(dimers, spheres, channels[i].nm, channels[i].radius,
                     mdl.box, channels[i].os, channels[i].sph_d, mdl.Nsph);
        }else{
          fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "channel",
                  cfg.num_channels);
        }
      }
    }

    // Make slit
    if(cfg.mk_slit){
      for(i=0; i<cfg.num_slits; i++){
        // Test slit data
        if(slits[i].nm[0] != 0 || slits[i].nm[1] != 0 || slits[i].nm[2] != 0){
          fprintf(stdout, " Inserting slit %d\n", i);
          make_slit(dimers, spheres, mdl.box, slits[i].thickness, slits[i].os,
                    slits[i].sph_d, slits[i].nm, mdl.Nsph);
        }else{
          fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "layer",
                  cfg.num_slits);
        }
      }
    }

    /* NOTE:
     * Brake dimers with spheres of type != 1
     * Set broken dimers as type '2'
     * set free (non-channel) spheres as type '3'
     */
    Nd2 = brake_dimers(dimers, spheres, mdl.Ndim);

    // Check the number of type-3 spheres
    Ns3 = count_typeX_spheres(spheres, 3, mdl.Nsph);

    // Display current statistics
    display_stats(mdl, Nd2, 2 * Nd2 - Ns3, Ns3);

    // Randomly (where possible) connect all neighbouring free spheres
    // into dimers. In critical cases start with spheres with fewest possible
    // connections
    if(Ns3 > 1 && cfg.mk_dimers == 1){

      fprintf(stdout,"\n\n Randomly connecting all possible free spheres\n");
      do{
        // Find a sphere with the lowes count of chanses to form a dimer
        fsi = find_critical_FS(spheres, mdl.Nsph);
        fsi_chances = count_typeX_sp_neighbours(spheres, 3, fsi);

        // If the possibilities are high enough, select sphere randomly
        if(fsi_chances >= 5){
          fsi = draw_sphere_typeX(spheres, 3, mdl.Nsph);
          fsi_chances = count_typeX_sp_neighbours(spheres, 3, fsi);
        }

        // Randomly select neighbour of fsi
        fsn = draw_ngb_sphere_typeX(spheres, 3, fsi);

        // Create a valid dimer from the pair of spheres
        if( fsi != -1 && fsn != -1){
          make_dimer(dimers, spheres, mdl.box, fsi, fsn, mdl.Ndim);
          // Update free spheres count
          Ns3 -= 2;
        }
      }while(fsi_chances > 0);

      test_dimer_distribution(dimers, Odistrib, mdl.Ndim);

      // Check the number of type-2 dimers
      Nd2 = count_typeX_dimers(dimers, 2, mdl.Ndim);

      // Check the number of type-3 spheres
      Ns3 = count_typeX_spheres(spheres, 3, mdl.Nsph);

      // Display current statistics
      display_stats(mdl, Nd2, 2 * Nd2 - Ns3, Ns3);
    }

    // Eliminate remaining free spheres (if any) using zipper
    if(Ns3 > 1 && cfg.mk_dimers == 1){
      fprintf(stdout, "\n\n Connecting remaining free spheres with zipper\n");

      // Check if there are even number of type-3 spheres
      Ns3_odd=0;
      Ns3_odd = Ns3 & 1;

      if(Ns3_odd != 0){
        fprintf(stdout, " NOTE:  Unable to connect all free spheres\n");
      }

      // Calculate the number of zipper runs to perform
      zip_Ns3_runs = (Ns3 - Ns3_odd)/2;

      fprintf(stdout, " Running zipper %d times to eliminate free spheres\n",
              zip_Ns3_runs);

      while(zip_Ns3_runs != 0){
        // Select type-3 sphere as the starting point for zipper
        for(i=0; i<mdl.Nsph; i++){
          if(spheres[i].type == 3){
            zip_init_sph = i;
            break;
          }
        }

        // Start zipper from selected sphere
        fprintf(stdout,"  zipper from type-%1d sph. no. %5d ",
                  spheres[zip_init_sph].type, zip_init_sph);

        fprintf(stdout," completed after %7d steps (%d left)\n",
                  zipper(dimers, spheres, mdl.box, mdl.Ndim, zip_init_sph, mdl.Nsph),
                  --zip_Ns3_runs);
      }

      test_dimer_distribution(dimers, Odistrib, mdl.Ndim);

      // Check the number of type-3 spheres
      Ns3 = count_typeX_spheres(spheres, 3, mdl.Nsph);

      // Check the number of type-2 dimers
      Nd2 = count_typeX_dimers(dimers, 2, mdl.Ndim);

      // Display current statistics
      display_stats(mdl, Nd2,2 * Nd2 - Ns3, Ns3);
    }   // Done eliminating type-3 spheres

    // Check and display structure parameters after inclusions setup (if any)
    if( (mdl.Ndim-Nd2) % 6 == 0){
      fprintf(stdout,"\n Perfect DC orientation possible (%d/dir.)\n",
              (mdl.Ndim-Nd2)/6);
    }else{
      fprintf(stdout,"\n Perfect DC orientation NOT possible (%.2lf/dir.)\n",
              ((double) (mdl.Ndim-Nd2))/6e0);
    }

    fprintf(stdout,fmt_dimer_distr_header, "Initial distribution:",
            "[110]", "[T10]", "[101]", "[T01]", "[011]", "[0T1]");

    // Check if the number of dimers in the system is high enough
    low_on_dimers = (((double) mdl.Ndim - Nd2)/((double) mdl.Ndim)  < 2e-1 ? 1 : 0);
    if(low_on_dimers != 0){
      fprintf(stdout, " Number of dimers to low to get a good structure\n");
    }

    // Reorganize molecules to get best distribution possible
    flip_count = 0;
    if( !validate_distrib(Odistrib, mdl.Ndim-Nd2, flip_count) && !low_on_dimers ){
      do{
        // Find a valid dimer configuration to flip orientations
        find_valid_cluster(dimers, spheres, mdl.box, mdl.Ndim, valid_dimer_pair);

        // Flip dimers
        if(valid_dimer_pair[0] != -1 && valid_dimer_pair[1] != -1){
          flip_dimers(dimers, spheres, mdl.box, Odistrib,
                      valid_dimer_pair[0], valid_dimer_pair[1]);
        }

        // Run zipper every once in a while (zipper length = X*num. of sph.)
        if(flip_count %1000000 == 0){
          do{
            // Select random type-1 sphere to start from
            zip_init_sph = (int) (u_RNG() * mdl.Nsph);
            zip_init_sph = (zip_init_sph < mdl.Nsph ? zip_init_sph : mdl.Nsph - 1);
          }while(spheres[zip_init_sph].type != 1);
          fprintf(stdout," Zipper %7d steps; distr.:",
            zipper(dimers, spheres, mdl.box, mdl.Ndim, zip_init_sph, 100*mdl.Nsph));
          // Display distribution after zipper
          test_dimer_distribution(dimers,  Odistrib, mdl.Ndim);
          display_dimer_distribution(Odistrib);
        }

        flip_count++;
      }while(validate_distrib(Odistrib, mdl.Ndim-Nd2, flip_count) == 0);
    }   // finished with flipping

    // Perform a final test of the structure prior to export
    test_dimer_distribution(dimers,  Odistrib, mdl.Ndim);
    fprintf(stdout," Step %9d distribution :", flip_count);
    display_dimer_distribution(Odistrib);

    // Generate required output files for structure 's'
    if( exp_str_data(cfg, mdl, dimers, spheres, Ns3, Nd2, s) != 0 ){
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }

  } // End structure loop

cleanup:
  // Free resources
  free(spheres);
  free(dimers);
  free(slits);
  free(channels);
  return exit_status;
}


/* vim: set tw=80 ts=2 sw=2 et: */
