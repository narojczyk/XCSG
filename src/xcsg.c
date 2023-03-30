/*
 * x-mer Crystal Structure Generator (XCSG)
 * Author: Jakub Narojczyk <narojczyk@ifmpan.poznan.pl>
 * (c) 2015-2021
 *
 * The program generates a crystal structure of f.c.c. spheres. The latter can
 * be connected into multimer molecules. The structure can also include
 * additional modifications like nanochannels, nanolayers or cluster inclusions.
 *
 */

/*
 * TODO:
 * BUG:  Setting rough_inclusions to 1 will make dimer inclusions regardless of *.dat file setting
 * BUG: reading version string from ini files is unsafe
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
#include "terminators.h"

int main(int argc, char *argv[])
{
  FILE *f = NULL;
  CONFIG cfg = {
    .seed = 123456789,
    .cells = {0, 0, 0},
    .first = 0,
    .last = 0,
    .mk_channel = 0,
    .mk_slit = 0,
    .mk_dimers = 0,
    .mk_inc_dimers = 0,
    .num_channels = 0,
    .num_slits = 0,
    .rough_inclusions = 0,
    .dimer_length = zero,
    .cfg_clusters = "none",
    .cfg_channels = "none",
    .cfg_slits = "none",
    .symmetry = "fcc"
  };

  MODEL mdl = {
    .Nsph = 0,
    .Ndim = 0,
    .mtrx_sph = 0,
    .mtrx_dim = 0,
    .incl_sph = 0,
    .incl_dim = 0,
    .box = {zero, zero, zero}
  };

  SPH *spheres = NULL;
  DIM3D *dimers = NULL;
  // Container struct for particles (valid pointers must be assigned
  // ... when memory has been alocated)
  PARTICLES particles = {NULL, NULL};

  INC *clusters = NULL;
  INC *channels = NULL;
  INC *slits = NULL;

  char *f_ini = NULL;

  int Odistrib[6] = {0,0,0,0,0,0};
  int i, s;
  int low_on_dimers = 0, insert_inclusion_dimers = 0;

  const char *fmt_missing_inclusion_normal =
    " [%s] WRN: missing normal vector for %s inc. no %d, skipping\n";
  const char *fmt_wrong_number_of_spheres =
    " [%s] ERR: failed to calculate the number of spheres\n";
  const char *fmt_wrong_container_dimensions =
    " [%s] ERR: failed to calculate continer dimensions\n";
  const char *fmt_generating_structure = "\n ### Processing structure %d\n";
  const char *fmt_inserting_inclusion = " Inserting %-10s of ID %2d\n";
  const char *fmt_inserting_rough_inclusions =
    "\n Inserting rough inclusions\n";
  const char *fmt_refine_DC_params =
    "\n Improve structure parameters\n"
    " Perfect DC orientation %s possible (%.2lf/dir.)\n";
  const char *fmt_matrix_dimers_to_low =
    " Number of dimers to low to get a good structure\n";

  // Extract program name from the path.
  prog_name = basename(argv[0]);

  // Parse command line options and get the config file name
  parse_options(argc, argv, &f_ini);

  // Open and parse config file
  if( f == NULL ){
    f = open_to_read(f_ini);
    if(parse_config(f, &cfg) != EXIT_SUCCESS){
        return EXIT_FAILURE;
    }
    f = NULL;
  }

  // Try to open and read cluster inclusions' description data
  if((f = try_open_to_read(cfg.cfg_clusters)) != NULL){
    if((clusters = parse_inclusions(f)) != NULL){
      on_exit(releace_memory, clusters);
      cfg.mk_clusters = 1;
      cfg.num_clusters = count_inclusions(clusters);
    }
    fclose(f);
    f = NULL;
  }

  // Try to open and read channel description data
  if((f = try_open_to_read(cfg.cfg_channels)) != NULL){
    if((channels = parse_inclusions(f)) != NULL){
      on_exit(releace_memory, channels);
      cfg.mk_channel = 1;
      cfg.num_channels = count_inclusions(channels);
    }
    fclose(f);
    f = NULL;
  }

  // Try to open and read slits description data
  if((f = try_open_to_read(cfg.cfg_slits)) != NULL){
    if((slits = parse_inclusions(f)) != NULL){
      on_exit(releace_memory, slits);
      cfg.mk_slit = 1;
      cfg.num_slits = count_inclusions(slits);
    }
    fclose(f);
    f = NULL;
  }

  // Basic inclusion parameters sanity check
  if(cfg.mk_clusters && inclusion_count_sanity_check(
    cfg.num_clusters, "number of point inclusions") == EXIT_FAILURE){
    return EXIT_FAILURE;
  }
  if(cfg.mk_channel && inclusion_count_sanity_check(
    cfg.num_channels, "number of channels") == EXIT_FAILURE){
    return EXIT_FAILURE;
  }
  if(cfg.mk_slit && inclusion_count_sanity_check(
    cfg.num_slits, "number of layers") == EXIT_FAILURE){
    return EXIT_FAILURE;
  }

  // Initiate generator with 'seed'
  init_RNG(cfg.seed);

  // Set the number of monomers, dimers, and possibly, further-mers
  mdl.Nsph = number_of_spheres(cfg);
  if(mdl.Nsph == EXIT_FAILURE){
    fprintf(stderr, fmt_wrong_number_of_spheres, prog_name);
    return EXIT_FAILURE;
  }
  mdl.Ndim = mdl.Nsph / 2;

  // Calculate container dimensions
  if(container_dimensions(&mdl, cfg) == EXIT_FAILURE){
    fprintf(stderr, fmt_wrong_container_dimensions, prog_name);
    return EXIT_FAILURE;
  }

  // Summary of configuration variables
  show_cfg_summary(cfg, mdl, slits, channels, clusters);

  // Allocate memory for spheres and dimers
  spheres = malloc( mdl.Nsph * sizeof(SPH));
  on_exit(releace_memory, spheres);
  dimers  = malloc( mdl.Ndim * sizeof(DIM3D));
  on_exit(releace_memory, dimers);
  // Copy pointers into the container structure
  particles.spheres = spheres;
  particles.dimers = dimers;

  // Loop over selected set of structures
  for(s=cfg.first; s<=cfg.last; s++){
    fprintf(stdout, fmt_generating_structure, s);

    // Clean allocated memory
    memory_clean_spheres(spheres, mdl.Nsph);
    memory_clean_dimers(dimers, mdl.Ndim);

// #### CASE 01 > ##############################################################
    // Set requested structure
    if(set_structure(cfg, mdl, spheres) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }
    // Update model structure with position of the sphere with smallest coords
    find_smallest_coordinates(spheres, &mdl);

    // Find neighbors for spheres
    if(find_ngb_spheres(spheres, mdl.Nsph, mdl.box) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }

    // Assign lattice indexes to all spheres
    if(sph_assign_lattice_indexes(spheres, mdl.Nsph) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }
// ############################################################## > CASE 01 ####

// #### CASE 02 > ##############################################################
    // Insert inclusions into sphere system
    // Make point inclusions
    if(cfg.mk_clusters && ! cfg.rough_inclusions){
      for(i=0; i<cfg.num_clusters; i++){
        // Test channel data
        if(clusters[i].radius > 0e0){
          fprintf(stdout, fmt_inserting_inclusion, "cluster", i);
          if(clusters[i].tgt_Nmer == 2) insert_inclusion_dimers = 1;
          make_cluster(mdl, &clusters[i], particles);
        }else{
          fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "cluster",
                  cfg.num_channels);
        }
      }
    }

    // Make channel
    if(cfg.mk_channel && ! cfg.rough_inclusions){
      for(i=0; i<cfg.num_channels; i++){
        // Test channel data
        if(channels[i].nm[0] != 0 || channels[i].nm[1] != 0 ||
          channels[i].nm[2] != 0){
          fprintf(stdout, fmt_inserting_inclusion, "channel", i);
          if(channels[i].tgt_Nmer == 2) insert_inclusion_dimers = 1;
          make_channel(mdl, &channels[i], particles);
        }else{
          fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "channel",
                  cfg.num_channels);
        }
      }
    }

    // Make slit
    if(cfg.mk_slit  && ! cfg.rough_inclusions){
      for(i=0; i<cfg.num_slits; i++){
        // Test slit data
        if(slits[i].nm[0] != 0 || slits[i].nm[1] != 0 || slits[i].nm[2] != 0){
          fprintf(stdout, fmt_inserting_inclusion, "layer", i);
          if(slits[i].tgt_Nmer == 2) insert_inclusion_dimers = 1;
          make_slit(mdl, &slits[i], particles);
        }else{
          fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "layer",
                  cfg.num_slits);
        }
      }
    }

    // Check the number of the respective types of particles
    if(count_particles_by_type(&mdl, particles) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }
// ############################################################## > CASE 02 ####
    // Display current statistics
    show_particle_stats(mdl,cfg);

// #### CASE 06 > ##############################################################
    // Create dimers inside the inclusions
    if(insert_inclusion_dimers){
      // Connect spheres in random pairs ...
      introduce_random_dimers(particles, &mdl, TYPE_INCLUSION_SPHERE,
                              TYPE_INCLUSION_DIMER);
      // ... and use zipper to finish the job if need be.
      introduce_dimers_by_zipper(particles, &mdl, TYPE_INCLUSION_BASE);

      // Display current statistics
      show_particle_stats(mdl,cfg);
    }
// ############################################################## > CASE 06 ####


// #### CASE 03 > ##############################################################
    // Create dimers in the matrix
    if(cfg.mk_dimers == 1 || cfg.rough_inclusions == 1){
      // Inserting dimers - Method #1
      if(mdl.mtrx_sph > 1){
        // Randomly (where possible) connect neighbouring spheres into dimers.
        // In critical cases start with spheres with fewest possible connections.
        introduce_random_dimers(particles, &mdl, TYPE_SPHERE, TYPE_DIMER);

        if(test_dimer_distribution(dimers, Odistrib, mdl.Ndim) == EXIT_FAILURE){
          return EXIT_FAILURE;
        }

        // Display current statistics
        show_particle_stats(mdl, cfg);
      }

      // Inserting dimers - Method #2
      if(mdl.mtrx_sph > 1){
        // Eliminate remaining spheres (if necesarry) using zipper
        introduce_dimers_by_zipper(particles, &mdl, TYPE_MATRIX_BASE);

        if(test_dimer_distribution(dimers, Odistrib, mdl.Ndim) == EXIT_FAILURE){
          return EXIT_FAILURE;
        }

        // Display current statistics
        show_particle_stats(mdl, cfg);
      } // Done eliminating type-spheres

      // ## Make a good DC structure
      fprintf(stdout, fmt_refine_DC_params,
              (mdl.mtrx_dim % 6 == 0) ? "" : "not",
              ((double) mdl.mtrx_dim)/6e0);

      // Check if the number of dimers in the system is high enough
      low_on_dimers = (((double) mdl.mtrx_dim)/((double) mdl.Ndim)
          < minimal_dimer_qty ? 1 : 0);

      if(low_on_dimers != 0){
        fprintf(stdout, fmt_matrix_dimers_to_low);
      }else{
        // Reorganize molecules to get best possible distribution
        if(refine_dimer_distribution(particles, mdl, Odistrib)
            == EXIT_FAILURE){
          return EXIT_FAILURE;
        }
      }

      // Perform a final test of the structure prior to export
      if(test_dimer_distribution(dimers, Odistrib, mdl.Ndim) == EXIT_FAILURE){
        return EXIT_FAILURE;
      }

// #### CASE 04 > ##############################################################
      // Insert inclusions into dimer system
      if(cfg.rough_inclusions){
        fprintf(stdout, fmt_inserting_rough_inclusions);
        // Make channel
        if(cfg.mk_channel){
          for(i=0; i<cfg.num_channels; i++){
            // Test channel data
            if(channels[i].nm[0] != 0 || channels[i].nm[1] != 0 ||
              channels[i].nm[2] != 0){
              fprintf(stdout, fmt_inserting_inclusion, "channel", i);
              make_channel(mdl, &channels[i], particles);
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
            if(slits[i].nm[0] != 0 || slits[i].nm[1] != 0 ||
              slits[i].nm[2] != 0){
              fprintf(stdout, fmt_inserting_inclusion, "layer", i);
              make_slit(mdl, &slits[i], particles);
            }else{
              fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "layer",
                      cfg.num_slits);
            }
          }
        }

        // Brake matrix dimers if not wanted in the model
        if(cfg.mk_dimers == 0){
          // TODO: Move this to a function
          for(i=0; i<mdl.Nsph; i++){
            if(spheres[i].type == TYPE_DIMER){
              spheres[i].type = TYPE_SPHERE;
              dimers[ spheres[i].dim_ind ].type = TYPE_INVALID;
            }
          }
          // Tidy up the dimer array
          int d=0;
          for(i=0; i<mdl.Ndim; i++){
            if(dimers[i].type == TYPE_INVALID) d = i;
            for(int j = i+1; j<mdl.Ndim; j++){
              if(dimers[j].type != TYPE_INVALID){
                dimers[d] = dimers[j];
                dimers[j].type = TYPE_INVALID;
                break;
              }
            }
          }
        }
// ############################################################## > CASE 04 ####
        // Check the number of different spheres
        if(count_particles_by_type(&mdl, particles) == EXIT_FAILURE){
          return EXIT_FAILURE;
        }
        // Display current statistics
        show_particle_stats(mdl, cfg);
      }
    } // end if(cfg.mk_dimers)
// ############################################################## > CASE 03 ####

    // Adjust dimer lengths if required
    update_dimer_lengths(mdl, cfg, particles);


    // Generate required output files for structure 's'
    if(export_structure(cfg, mdl, particles, s) != EXIT_SUCCESS ){
      return EXIT_FAILURE;
    }
  } // End structure loop

  // Succesfully finish program
  return EXIT_SUCCESS;
}

/* vim: set tw=80 ts=2 sw=2 et: */
