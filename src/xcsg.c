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
#include "terminators.h"

int main(int argc, char *argv[])
{
  FILE *f = NULL;
  CONFIG cfg;
  MODEL mdl;

  SPH *spheres = NULL;
  DIM3D *dimers = NULL;
  INC *channels = NULL;
  INC *slits = NULL;

  char *f_ini = NULL;

  int zipper_runs=0, zip_init_sph=-1;
  int Odistrib[6] = {0,0,0,0,0,0};
  int flipable_dimers[2] = {-1, -1};
  int i, s;
  int low_on_dimers = 0;
  int flip_count = 0;
  int s_id, s_ngb_id, s_ngb_qty;

  const char *fmt_missing_inclusion_normal =
    " [%s] WRN: missing normal vector for %s inc. no %d, skipping\n";
  const char *fmt_dimer_distr_header =
    " %-28s %5s %5s %5s %5s %5s %5s\n";
  const char *fmt_wrong_number_of_spheres =
    " [%s] ERR: failed to calculate the number of spheres\n";
  const char *fmt_wrong_container_dimensions =
    " [%s] ERR: failed to calculate continer dimensions\n";
  const char *fmt_generating_structure = "\n ### Processing structure %d\n";
  const char *fmt_inserting_inclusion = " Inserting %-10s of ID %2d\n";
  const char *fmt_mk_dimers_at_random =
    " Randomly join spheres into dimers\n";
  const char *fmt_mk_dimers_with_zipper =
    " Connecting remaining spheres with zipper\n";
  const char *fmt_required_zipper_runs =
    " Running zipper %d times to create dimers\n";
  const char *fmt_ziper_run =
    " Zipper from sphere ID %5d completed after %7d steps\n";
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

  // Allocate and clean memory for channels' data
  channels = malloc(cfg.num_channels * sizeof(INC));
  on_exit(releace_memory, channels);
  memory_clean_inclusion(channels, cfg.num_channels);

  // Allocate and clean memory for slits' data
  slits = malloc(cfg.num_slits * sizeof(INC));
  on_exit(releace_memory, slits);
  memory_clean_inclusion(slits, cfg.num_slits);

  // Open and read channel description data
  if( cfg.mk_channel != 0 && f == NULL ){
    f = open_to_read(cfg.cfg_channels);
    if(parse_inclusions(f, channels, cfg.num_channels) != EXIT_SUCCESS){
      return EXIT_FAILURE;
    }
    f = NULL;
  }

  // Open and read slits description data
  if( cfg.mk_slit != 0 && f == NULL ){
    f = open_to_read(cfg.cfg_slits);
    if(parse_inclusions(f, slits, cfg.num_slits) != EXIT_SUCCESS){
      return EXIT_FAILURE;
    }
    f = NULL;
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
  display_configuration_summary(cfg, mdl, slits, channels);

  // Allocate memory for spheres and dimers
  spheres = malloc( mdl.Nsph * sizeof(SPH));
  on_exit(releace_memory, spheres);
  dimers  = malloc( mdl.Ndim * sizeof(DIM3D));
  on_exit(releace_memory, dimers);

  // Loop over selected set of structures
  for(s=cfg.first; s<=cfg.last; s++){
    fprintf(stdout, fmt_generating_structure, s);

    // Clean allocated memory
    memory_clean_spheres(spheres, mdl.Nsph);
    memory_clean_dimers(dimers, mdl.Ndim);

    // Set requested structure
    if(set_structure(cfg, mdl, spheres) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }

    // Find neighbors for spheres
    if(find_ngb_spheres(spheres, mdl.Nsph, mdl.box) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }

    // Assign lattice indexes to all spheres
    if(sph_assign_lattice_indexes(spheres, mdl.Nsph) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }

    // Make channel
    if(cfg.mk_channel){
      for(i=0; i<cfg.num_channels; i++){
        // Test channel data
        if(channels[i].nm[0] != 0 || channels[i].nm[1] != 0 || channels[i].nm[2] != 0){
          fprintf(stdout, fmt_inserting_inclusion, "channel", i);
          make_channel(mdl, dimers, spheres, &channels[i]);
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
          fprintf(stdout, fmt_inserting_inclusion, "layer", i);
          make_slit(mdl, dimers, spheres, &slits[i]);
        }else{
          fprintf(stderr, fmt_missing_inclusion_normal, prog_name, "layer",
                  cfg.num_slits);
        }
      }
    }

    // Check the number of different spheres
    if(count_particles_by_type(&mdl, spheres, dimers) == EXIT_FAILURE){
      return EXIT_FAILURE;
    }

    // Display current statistics
    display_stats(mdl,cfg);

    // Create dimers
    if(mdl.mtrx_sph > 1 && cfg.mk_dimers == 1){
      // Randomly (where possible) connect all neighbouring spheres into dimers.
      // In critical cases start with spheres with fewest possible connections.
      fprintf(stdout, fmt_mk_dimers_at_random);
      do{
        // Find a sphere with the lowes count of chanses to form a dimer
        s_id = find_critical_sphere(spheres, mdl.Nsph);
        s_ngb_qty =
          count_typeX_sp_neighbours(spheres, TYPE_SPHERE, s_id, mdl.Nsph);
        // If the possibilities are high enough, select sphere randomly
        if(s_ngb_qty >= 5){
          s_id = draw_sphere_typeX(spheres, TYPE_SPHERE, mdl.Nsph);
          s_ngb_qty =
            count_typeX_sp_neighbours(spheres, TYPE_SPHERE, s_id, mdl.Nsph);
        }

        // Randomly select neighbour of s_id
        s_ngb_id = draw_ngb_sphere_typeX(spheres, TYPE_SPHERE, s_id);

        // Create a valid dimer from the pair of spheres
        if( s_id != -1 && s_ngb_id != -1){
          make_dimer(dimers, spheres, mdl, s_id, s_ngb_id);
          // Update free spheres count
          mdl.mtrx_sph -= 2;
        }
      }while(s_ngb_qty > 0);

      if(test_dimer_distribution(dimers, Odistrib, mdl.Ndim) == EXIT_FAILURE){
        return EXIT_FAILURE;
      }

      if(count_particles_by_type(&mdl, spheres, dimers) == EXIT_FAILURE){
        return EXIT_FAILURE;
      }

      // Display current statistics
      display_stats(mdl, cfg);

      // Eliminate remaining spheres (if necesarry) using zipper
      if(mdl.mtrx_sph > 1){
        fprintf(stdout, fmt_mk_dimers_with_zipper);

        // Calculate the number of zipper runs to perform
        zipper_runs = (mdl.mtrx_sph - (mdl.mtrx_sph & 1))/2;
        fprintf(stdout, fmt_required_zipper_runs, zipper_runs);

        while(zipper_runs != 0){
          // Select type-sphere object as the starting point for zipper
          for(i=0; i<mdl.Nsph; i++){
            if(spheres[i].type == TYPE_SPHERE){
              zip_init_sph = i;
              break;
            }
          }

          // Start zipper from selected sphere
          fprintf(stdout,"  (%d)", zipper_runs);
          fprintf(stdout, fmt_ziper_run, zip_init_sph,
                    zipper(mdl, dimers, spheres, zip_init_sph, lazy));
          zipper_runs--;
        }

        if(test_dimer_distribution(dimers, Odistrib, mdl.Ndim) == EXIT_FAILURE){
          return EXIT_FAILURE;
        }

        if(count_particles_by_type(&mdl, spheres, dimers) == EXIT_FAILURE){
          return EXIT_FAILURE;
        }

        // Display current statistics
        display_stats(mdl, cfg);
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
      }

      // Reorganize molecules to get best possible distribution
      flip_count = 0;
      fprintf(stdout,fmt_dimer_distr_header, "Initial distribution:",
              "[110]", "[T10]", "[101]", "[T01]", "[011]", "[0T1]");
      if(!validate_distrib(Odistrib, mdl.mtrx_dim, flip_count) &&
         !low_on_dimers ){
        do{
          // Find a valid dimer configuration to flip orientations
          find_valid_cluster(mdl, dimers, spheres, flipable_dimers);
          // Flip dimers
          if(flipable_dimers[0] != -1 && flipable_dimers[1] != -1){
            flip_dimers(mdl, dimers, spheres, Odistrib, flipable_dimers);
          }
          // Run zipper every once in a while (zipper length = X*num. of sph.)
          if(flip_count % run_zipper == 0){
            do{
              // Select random sphere (forming any dimer) to start from
              zip_init_sph = (int) (u_RNG() * mdl.Nsph);
              zip_init_sph = (zip_init_sph < mdl.Nsph ? zip_init_sph : mdl.Nsph - 1);
            }while(spheres[zip_init_sph].type != TYPE_SPHERE_DIMER);
            fprintf(stdout, fmt_ziper_run, zip_init_sph,
                    zipper(mdl, dimers, spheres, zip_init_sph, diligent));
            // Display distribution after zipper
            if(test_dimer_distribution(dimers,  Odistrib, mdl.Ndim)
              == EXIT_FAILURE){
              return EXIT_FAILURE;
            }
            fprintf(stdout," %-27s :", "Zipper distribution");
            display_dimer_distribution(Odistrib);
          }

          flip_count++;
        }while(validate_distrib(Odistrib, mdl.mtrx_dim, flip_count) == 0);
      }   // finished with flipping

      // Perform a final test of the structure prior to export
      if(test_dimer_distribution(dimers,  Odistrib, mdl.Ndim) == EXIT_FAILURE){
        return EXIT_FAILURE;
      }
    } // end if(cfg.mk_dimers)

    // Generate required output files for structure 's'
    if( exp_str_data(cfg, mdl, dimers, spheres, s) != EXIT_SUCCESS ){
      return EXIT_FAILURE;
    }
  } // End structure loop

  // Succesfully finish program
  return EXIT_SUCCESS;
}

/* vim: set tw=80 ts=2 sw=2 et: */
