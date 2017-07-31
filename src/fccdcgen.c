/*
 * fccdcgen
 * f.c.c. structure of DC dimers generator/translator
 * Author: Jakub Narojczyk <narojczyk@ifmpan.poznan.pl>
 * Author: Mikolaj Kowalik <kowalik@ifmpan.poznan.pl>
 *
 * The program generates a DC structure of dimers on f.c.c. spheres with the
 * presence of (optional) nano-channel of spheres.
 *
 * The program translates an existing DC structure generated by 'dscgen' by
 * Mikolaj Kowalik to conform the new data structure of the simulation program
 *
 * (c) 2015
 */

/*
 *
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

  SPH *spheres = NULL;
  DIM3D *dimers = NULL;
  CHA *channels = NULL;

  char *f_ini = NULL;
  char f_inp[15];

  const char *fo_iDC = "initdc3d.%04d";

  int exit_status;
  int Ns, Nd, Nd2, Ns3, Ns3_odd;
  int zip_Ns3_runs=0, zip_init_sph=-1;
  int Odistrib[6] = {0,0,0,0,0,0};
  int valid_dimer_pair[2];
  int i, s;
  int low_on_dimers = 0;
  int flip_count = 0;

  double cube_edge[3];
  

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
  parse_config(f);
  fclose(f);
  
  // Allocate memory for channels' data
  channels = malloc( i_n_channels * sizeof(CHA));
  
  // Open and read channel description data
  if( i_make_channel != 0 ){
    if((f = fopen(i_chdesc_file, "r")) == NULL) {
      fprintf(stderr, "  [%s]: error: cannot open channels file: %s\n",
              prog_name, i_chdesc_file);
      return EXIT_FAILURE;
    }
    exit_status = parse_channels(f, channels);
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
  cube_edge[0] = i_edge_fcc_N[0] * sqrt(two);
  cube_edge[1] = i_edge_fcc_N[1] * sqrt(two);
  cube_edge[2] = i_edge_fcc_N[2] * sqrt(two);

  // Summary of configuration variables
  fprintf(stdout,"\tConfig summary\n");
  fprintf(stdout," System size (cells)  : %d by %d by %d\n",
          i_edge_fcc_N[0], i_edge_fcc_N[1], i_edge_fcc_N[2]);
  fprintf(stdout," MT19937 seed         : %8lu\n", i_seed);
  fprintf(stdout," Str. index range     : %d to %d (total %d files)\n",
          i_iDCfrom, i_iDCto, i_iDCto-i_iDCfrom+1);

  fprintf(stdout,"\n\tOther parameters\n");
  fprintf(stdout," Coordinates range    : %.16le %.16le (x)\n",
          -cube_edge[0]/two,cube_edge[0]/two);
  fprintf(stdout,"                        %.16le %.16le (y)\n",
          -cube_edge[1]/two,cube_edge[1]/two);
  fprintf(stdout,"                        %.16le %.16le (z)\n",
          -cube_edge[2]/two,cube_edge[2]/two);
  fprintf(stdout," Box dimensions       : %.16le (x)\n", cube_edge[0]);
  fprintf(stdout,"                        %.16le (y)\n", cube_edge[1]);
  fprintf(stdout,"                        %.16le (z)\n", cube_edge[2]);
  fprintf(stdout," Number of dimers     : %d\n",Nd);
  fprintf(stdout," Number of spheres    : %d\n",Ns);
  fprintf(stdout," Number of channels   : %d\n",i_n_channels);
  if (i_make_channel){
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

  // Allocate memory for spheres and dimers
  spheres = malloc( Ns * sizeof(SPH));
  dimers  = malloc( Nd * sizeof(DIM3D));

  // Loop over selected set of structures
  for(s=i_iDCfrom; s<=i_iDCto; s++){

    // Clean allocated memory
    memory_clean_spheres(spheres, Ns);
    memory_clean_dimers(dimers, Nd);

    fprintf(stdout,"\n ***\tProcessing structure %d\n",s);
    // Regarding the ini settings, load input structure or grnerate a new one
    if(i_iDCfrom >= 0 && i_iDCto >= i_iDCfrom){
      sprintf(f_inp, fo_iDC, s);
      fprintf(stdout," Reading structure from file %s\n",f_inp);
      // Load existing DC structure
      if((f = fopen(f_inp, "r")) == NULL) {
        fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
                prog_name, f_inp);
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }
      if(load_dcsgen(f, dimers, cube_edge, Nd) != Nd){
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }
    }else{
      fprintf(stdout," Generating new structure\n");
      // Set fcc structure of spheres
      sph_set_fcc( spheres, Ns, i_edge_fcc_N);
      // TODO: Set initial dimer structure here

    }

    // Bind dimers to spheres and vice versa
    bind_spheres_to_dimers(dimers, spheres, Nd);

    // Generate sphere positions for all dimers
    for(i=0; i<Nd; i++){
      update_sphere_positions(dimers, spheres, cube_edge, i);
    }

    // Find neighbors for spheres
    if(find_ngb_spheres(spheres, Ns, cube_edge) != 0){
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }

    // Make channel    
    if(i_make_channel){
      for(i=0; i<i_n_channels; i++){
        // Test channel data
        if(channels[i].nm[0]+channels[i].nm[1]+channels[i].nm[2] != 0){
          fprintf(stdout, " Inserting %d channel(s)\n", i_n_channels);
          make_channel(dimers, spheres, channels[i].nm, channels[i].radius, 
                     cube_edge, channels[i].os, channels[i].sph_d, Nd);
        }else{
          fprintf(stderr, 
                  " [ERR] Missing normal vector for channel %d, skipping\n", 
                  i_n_channels);
        }
      }       
    }

    // Make slit
    if(i_make_slit && i_normal[0]+i_normal[1]+i_normal[2] > 0){
      make_slit(dimers, spheres, i_slit_Th, i_normal, Nd);
    }

    /* NOTE:
     * Brake dimers with spheres of type != 1
     * Set broken dimers as type '2'
     * set free (non-channel) spheres as type '3'
     */
    Nd2 = brake_dimers(dimers, spheres, Nd);

    // Check the number of type-3 spheres
    Ns3 = 0;
    for(i=0; i<Ns; i++){
      if(spheres[i].type == 3){
        Ns3++;
      }
    }

    fprintf(stdout, " Dimers broken   (if any): %4d %6.2lf %%\n",
            Nd2, (1e2*Nd2)/(1e0*Nd) );
    fprintf(stdout, " Channel spheres (if any): %4d %6.2lf %%\n",
            2 * Nd2 - Ns3, (1e2*(2 * Nd2 - Ns3))/(1e0*Ns));
    fprintf(stdout, " Free spheres    (if any): %4d %6.2lf %%\n",
            Ns3, (1e2*Ns3)/(1e0*Ns) );

    if(Ns3 > 0 && i_fs_connect == 1){
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
                  zipper(dimers, spheres, cube_edge, Nd, zip_init_sph, Ns),
                  --zip_Ns3_runs);
      }

      // Check the number of type-3 spheres
      Ns3 = 0;
      for(i=0; i<Ns; i++){
        if(spheres[i].type == 3){
          Ns3++;
        }
      }
      // Check the number of type-2 dimers
      Nd2=0;
      for(i=0; i<Nd; i++){
        if(dimers[i].type == 2){
          Nd2++;
        }
      }

      fprintf(stdout, "\n Dimers valid            : %4d %6.2lf %%\n",
              Nd-Nd2, (1e2*(Nd-Nd2))/(1e0*Nd) );
      fprintf(stdout, " Dimers broken           : %4d %6.2lf %%\n",
              Nd2, (1e2*Nd2)/(1e0*Nd) );
      fprintf(stdout, " Free spheres            : %4d %6.2lf %%\n",
              Ns3, (1e2*Ns3)/(1e0*Ns) );
    }   // Done eliminating type-3 spheres

    // Check and display structure parameters after channel setup (if any)
    if( (Nd-Nd2) % 6 == 0){
      fprintf(stdout,"\n Perfect DC orientation possible (%d/dir.)\n",
              (Nd-Nd2)/6);
    }else{
      fprintf(stdout,"\n Perfect DC orientation NOT possible (%.2lf/dir.)\n",
              ((double) (Nd-Nd2))/6e0);
    }

    test_dimer_distribution(dimers, Odistrib, Nd);

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
        find_valid_cluster(dimers, spheres, cube_edge, Nd, valid_dimer_pair);

        // Flip dimers
        if(valid_dimer_pair[0] != -1 && valid_dimer_pair[1] != -1){
          flip_dimers(dimers, spheres, cube_edge, Odistrib, 
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
            zipper(dimers, spheres, cube_edge, Nd, zip_init_sph, 100*Ns));
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
    if( exp_str_data(dimers, spheres, cube_edge, Ns, Nd, s) != 0 ){
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }
        
  } // End structure loop

cleanup:
  // Free resources
  free(spheres);
  free(dimers);
  free(channels);
  return exit_status;
}


/* vim: set tw=80 ts=2 sw=2 et: */
