/*
 * mms
 * f.c.c. structure of DC dimers generator/translator
 * Author: Jakub Narojczyk <narojczyk@ifmpan.poznan.pl>
 *
 * The program generates a DC structure of dimers on f.c.c. spheres with the
 * presence of (optional) nano-channel of spheres.
 *
 *
 * (c) 2015
 */

/*
 *
 * TODO:
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

  SPH *spheres = NULL;
  DIM3D *dimers = NULL;
  CHA *channels = NULL;
  SLI *slits = NULL;

  char *f_ini = NULL;

  int exit_status;
  int Ns, Nd, Nd2, Ns3, Ns3_odd;
  int zip_Ns3_runs=0, zip_init_sph=-1;
  int Odistrib[6] = {0,0,0,0,0,0};
  int valid_dimer_pair[2];
  int i, s;
  int low_on_dimers = 0;
  int flip_count = 0;

  double box_edge[3];

  int fsi, fsn, fsi_chances;


  // Extract program name from the path.
  prog_name = basename(argv[0]);

  // For now, assume everything will be fine.
  exit_status = EXIT_SUCCESS;

  // Parse command line options and get the config file name
  parse_options(argc, argv, &f_ini);

  // Open and parse config file
  if((f = fopen(f_ini, "r")) == NULL) {
    fprintf(stderr, "  [%s]: error: cannot open config file: %s\n",
            prog_name, f_ini);
    return EXIT_FAILURE;
  }
  exit_status = parse_config(f);
  fclose(f);
  if(exit_status != EXIT_SUCCESS){
      goto cleanup;
  }

  // Allocate memory for channels' data
  channels = malloc( i_n_channels * sizeof(CHA));

  // Allocate memory for slits' data
  slits = malloc( i_n_slits * sizeof(SLI));

  // Open and read channel description data
  if( i_make_channel != 0 ){
    // Clean alocated memory for channels
    memory_clean_channels(channels, i_n_channels);

    if((f = fopen(i_Fchannels, "r")) == NULL) {
      fprintf(stderr, "  [%s]: error: cannot open channels file: %s\n",
              prog_name, i_Fchannels);
      return EXIT_FAILURE;
    }
    exit_status = parse_channels(f, channels);
    fclose(f);
    if(exit_status != EXIT_SUCCESS){
      goto cleanup;
    }
  }

  // Open and read slits description data
  if( i_make_slit != 0 ){
    // Clean alocated memory for channels
    memory_clean_slits(slits, i_n_slits);

    if((f = fopen(i_Fslits, "r")) == NULL) {
      fprintf(stderr, "  [%s]: error: cannot open slits file: %s\n",
              prog_name, i_Fslits);
      return EXIT_FAILURE;
    }
    exit_status = parse_slits(f, slits);
    fclose(f);
    if(exit_status != EXIT_SUCCESS){
      goto cleanup;
    }
  }

  // Initiate generator with 'seed'
  init_RNG(i_seed);

  // Set the number of spheres and dimres
  Ns = 4 * i_edge_fcc_N[0] * i_edge_fcc_N[1] * i_edge_fcc_N[2];
  Nd = Ns / 2;

  // Calculate cube edge
  box_edge[0] = i_edge_fcc_N[0] * sqrt(two);
  box_edge[1] = i_edge_fcc_N[1] * sqrt(two);
  box_edge[2] = i_edge_fcc_N[2] * sqrt(two);

  // Summary of configuration variables
  fprintf(stdout,"\tConfig summary\n");
  fprintf(stdout," System size (cells)  : %d by %d by %d\n",
          i_edge_fcc_N[0], i_edge_fcc_N[1], i_edge_fcc_N[2]);
  fprintf(stdout," MT19937 seed         : %8lu\n", i_seed);
  fprintf(stdout," Str. index range     : %d to %d (total %d files)\n",
          i_iDCfrom, i_iDCto, i_iDCto-i_iDCfrom+1);

  fprintf(stdout,"\n\tOther parameters\n");
  fprintf(stdout," Coordinates range    : %.16le %.16le (x)\n",
          -box_edge[0]/two,box_edge[0]/two);
  fprintf(stdout,"                        %.16le %.16le (y)\n",
          -box_edge[1]/two,box_edge[1]/two);
  fprintf(stdout,"                        %.16le %.16le (z)\n",
          -box_edge[2]/two,box_edge[2]/two);
  fprintf(stdout," Box dimensions       : %.16le (x)\n", box_edge[0]);
  fprintf(stdout,"                        %.16le (y)\n", box_edge[1]);
  fprintf(stdout,"                        %.16le (z)\n", box_edge[2]);
  fprintf(stdout," Number of dimers     : %d\n",Nd);
  fprintf(stdout," Number of spheres    : %d\n",Ns);
  fprintf(stdout," Number of channels   : %d\n",i_n_channels);
  if (i_make_channel){
    fprintf(stdout," Channels data from   : %s\n",i_Fchannels);
    fprintf(stdout," No.\t channell offset\t\t chanel normal\t\tradius"
      "\tsph. diameter\n");
    for(i=0; i<i_n_channels; i++){
      fprintf(stdout," %d | %lf %lf %lf | %lf %lf %lf | %lf | %lf\n", i,
              channels[i].os[0], channels[i].os[1], channels[i].os[2],
              channels[i].nm[0], channels[i].nm[1], channels[i].nm[2],
              channels[i].radius, channels[i].sph_d);
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout," Number of slits      : %d\n",i_n_slits);
  if (i_make_slit){
    fprintf(stdout," Slit data from       : %s\n",i_Fslits);
    fprintf(stdout," No.\tslit offset\t\t\tslit normal\t\tthick."
      "\t sph. diam.\n");
    for(i=0; i<i_n_slits; i++){
      fprintf(stdout," %d | %lf %lf %lf | %lf %lf %lf | %lf | %lf\n", i,
              slits[i].os[0], slits[i].os[1], slits[i].os[2],
              slits[i].nm[0], slits[i].nm[1], slits[i].nm[2],
              slits[i].thickness, slits[i].sph_d);
    }
    fprintf(stdout,"\n");
  }

  // Allocate memory for spheres and dimers
  spheres = malloc( Ns * sizeof(SPH));
  dimers  = malloc( Nd * sizeof(DIM3D));

  // Loop over selected set of structures
  for(s=i_iDCfrom; s<=i_iDCto; s++){

    // Clean allocated memory
    memory_clean_spheres(spheres, Ns);
    memory_clean_dimers(dimers, Nd);

    fprintf(stdout,"\n ***\tProcessing structure %d\n",s);
    fprintf(stdout," Generating pure f.c.c. structure\n");

    // Set fcc structure of spheres
    sph_set_fcc( spheres, Ns, i_edge_fcc_N);

    // Find neighbors for spheres
    if(find_ngb_spheres(spheres, Ns, box_edge) != 0){
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }

    // Assign lattice indexes to all spheres
    sph_assign_lattice_indexes(spheres, 0, Ns);
    sph_assign_lattice_indexes(spheres, 1, Ns);
    sph_assign_lattice_indexes(spheres, 2, Ns);

    // Flag all dimers as broken at this point
    for(i=0; i<Nd; i++){
      dimers[i].type = 2;
    }

    // Make channel
    if(i_make_channel){
      for(i=0; i<i_n_channels; i++){
        // Test channel data
        if(channels[i].nm[0] != 0 || channels[i].nm[1] != 0 || channels[i].nm[2] != 0){
          fprintf(stdout, " Inserting channel %d\n", i);
          make_channel(dimers, spheres, channels[i].nm, channels[i].radius,
                     box_edge, channels[i].os, channels[i].sph_d, Ns);
        }else{
          fprintf(stderr,
                  " [ERR] Missing normal vector for channel %d, skipping\n",
                  i_n_channels);
        }
      }
    }

    // Make slit
    if(i_make_slit){
      for(i=0; i<i_n_slits; i++){
        // Test slit data
        if(slits[i].nm[0] != 0 || slits[i].nm[1] != 0 || slits[i].nm[2] != 0){
          fprintf(stdout, " Inserting slit %d\n", i);
          make_slit(dimers, spheres, box_edge, slits[i].thickness, slits[i].os,
                    slits[i].sph_d, slits[i].nm, Ns);
        }else{
          fprintf(stderr,
                  " [ERR] Missing normal vector for channel %d, skipping\n",
                  i_n_slits);
        }
      }
    }

    /* NOTE:
     * Brake dimers with spheres of type != 1
     * Set broken dimers as type '2'
     * set free (non-channel) spheres as type '3'
     */
    Nd2 = brake_dimers(dimers, spheres, Nd);

    // Check the number of type-3 spheres
    Ns3 = count_typeX_spheres(spheres, 3, Ns);

    // Display current statistics
    display_stats(Nd-Nd2, Nd2, 2 * Nd2 - Ns3, Ns3, Ns);

    // Randomly (where possible) connect all neighbouring free spheres
    // into dimers. In critical cases start with spheres with fewest possible
    // connections
    if(Ns3 > 1 && i_fs_connect == 1){

      fprintf(stdout,"\n\n Randomly connecting all possible free spheres\n");
      do{
        // Find a sphere with the lowes count of chanses to form a dimer
        fsi = find_critical_FS(spheres, Ns);
        fsi_chances = count_typeX_sp_neighbours(spheres, 3, fsi);

        // If the possibilities are high enough, select sphere randomly
        if(fsi_chances >= 5){
          fsi = draw_sphere_typeX(spheres, 3, Ns);
          fsi_chances = count_typeX_sp_neighbours(spheres, 3, fsi);
        }

        // Randomly select neighbour of fsi
        fsn = draw_ngb_sphere_typeX(spheres, 3, fsi);

        // Create a valid dimer from the pair of spheres
        if( fsi != -1 && fsn != -1){
          make_dimer(dimers, spheres, box_edge, fsi, fsn, Nd);
          // Update free spheres count
          Ns3 -= 2;
        }
      }while(fsi_chances > 0);

      test_dimer_distribution(dimers, Odistrib, Nd);

      // Check the number of type-2 dimers
      Nd2 = count_typeX_dimers(dimers, 2, Nd);

      // Check the number of type-3 spheres
      Ns3 = count_typeX_spheres(spheres, 3, Ns);

      // Display current statistics
      display_stats(Nd-Nd2, Nd2, 2 * Nd2 - Ns3, Ns3, Ns);
    }

    // Eliminate remaining free spheres (if any) using zipper
    if(Ns3 > 1 && i_fs_connect == 1){
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
        for(i=0; i<Ns; i++){
          if(spheres[i].type == 3){
            zip_init_sph = i;
            break;
          }
        }

        // Start zipper from selected sphere
        fprintf(stdout,"  zipper from type-%1d sph. no. %5d ",
                  spheres[zip_init_sph].type, zip_init_sph);

        fprintf(stdout," completed after %7d steps (%d left)\n",
                  zipper(dimers, spheres, box_edge, Nd, zip_init_sph, Ns),
                  --zip_Ns3_runs);
      }

      test_dimer_distribution(dimers, Odistrib, Nd);

      // Check the number of type-3 spheres
      Ns3 = count_typeX_spheres(spheres, 3, Ns);

      // Check the number of type-2 dimers
      Nd2 = count_typeX_dimers(dimers, 2, Nd);

      // Display current statistics
      display_stats(Nd-Nd2, Nd2,2 * Nd2 - Ns3, Ns3, Ns);
    }   // Done eliminating type-3 spheres

    // Check and display structure parameters after inclusions setup (if any)
    if( (Nd-Nd2) % 6 == 0){
      fprintf(stdout,"\n Perfect DC orientation possible (%d/dir.)\n",
              (Nd-Nd2)/6);
    }else{
      fprintf(stdout,"\n Perfect DC orientation NOT possible (%.2lf/dir.)\n",
              ((double) (Nd-Nd2))/6e0);
    }



    fprintf(stdout," %-28s %5s %5s %5s %5s %5s %5s\n",
          "Initial distribution:", "[110]", "[i10]", "[101]", "[i01]", "[011]", "[0i1]");
    // Displayed at the first flip or prior to structure export (if perfect
    // distribution is alredy present)

    // Check if the number of dimers in the system is high enough
    low_on_dimers = (((double) Nd - Nd2)/((double) Nd) < 2e-1 ? 1 : 0);
    if(low_on_dimers != 0){
      fprintf(stdout, " Number of dimers to low to get a good structure\n");
    }

    // Reorganize molecules to get best distribution possible
    flip_count = 0;
    if( !validate_distrib(Odistrib, Nd-Nd2, flip_count) && !low_on_dimers ){
      do{
        // Find a valid dimer configuration to flip orientations
        find_valid_cluster(dimers, spheres, box_edge, Nd, valid_dimer_pair);

        // Flip dimers
        if(valid_dimer_pair[0] != -1 && valid_dimer_pair[1] != -1){
          flip_dimers(dimers, spheres, box_edge, Odistrib,
                      valid_dimer_pair[0], valid_dimer_pair[1]);
        }

        // Run zipper every once in a while (zipper length = X*num. of sph.)
        if(flip_count %1000000 == 0){
          do{
            // Select random type-1 sphere to start from
            zip_init_sph = (int) (u_RNG() * Ns);
            zip_init_sph = (zip_init_sph < Ns ? zip_init_sph : Ns - 1);
          }while(spheres[zip_init_sph].type != 1);
          fprintf(stdout," Zipper %7d steps; distr.:",
            zipper(dimers, spheres, box_edge, Nd, zip_init_sph, 100*Ns));
          // Display distribution after zipper
          test_dimer_distribution(dimers,  Odistrib, Nd);
          display_dimer_distribution(Odistrib);
        }

        flip_count++;
      }while(validate_distrib(Odistrib, Nd-Nd2, flip_count) == 0);
    }   // finished with flipping

    // Perform a final test of the structure prior to export
    test_dimer_distribution(dimers,  Odistrib, Nd);
    fprintf(stdout," Step %9d distribution :", flip_count);
    display_dimer_distribution(Odistrib);

    // Generate required output files for structure 's'
    if( exp_str_data(dimers, spheres, box_edge, Ns, Nd, Ns3, Nd2, s) != 0 ){
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
