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

int main(int argc, char *argv[])
{

  FILE *f;
  
  SPH *spheres;
  DIM3D *dimers;

  char *f_ini = NULL;
  char f_out[15];
  char f_inp[15];
  
  const char *fo_iDC = "initdc3d.%03d";
  const char *fo_exp_dim = "d3d%05d.csv";
  const char *fo_exp_sph = "s3d%05d.csv";

  int exit_status;
  int Ns, Nd;
  int bd;
  
  int i, m, n;
  int s;
  /*
  int h, v, a;
  int next, ngb, master;
  int c=0, exp_c=0;
  int pc_of = -1, pc_with = -1;
  int target_index, initial_N;
  int terminate_simulation = 0;
  int cwd = 0;  //colliding walker deleted [0|1]
  int survivors = 0, surv_n = 0;

//   int shits;

  double pc_time;
  double sum_R = zero;
  double v1_x, v1_y, v2_x, v2_y;
  double x1, y1, x2, y2;
  double R1, R2;
  double del_vx, del_vy, del_x, del_y, relativeV_SQ, toSQRT;
  double tau = zero, d_tau = zero, ts1, ts2;

  double rr, vr;
  double rho = zero;
  double tvx, tvy;
  double pmarker = zero;
*/
  double cube_edge;
  
  // Extract program name from the path.
  prog_name = basename(argv[0]);

  // For now, assume everything will be fine.
  exit_status = EXIT_SUCCESS;

  // Parse command line options and get the config file name
  parse_options(argc, argv, &f_ini);

  // Open and parse config file
  if((f = fopen(f_ini, "r")) == NULL) {
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            prog_name, f_ini);
    return EXIT_FAILURE;
  }
  parse_config(f);
  fclose(f);
  
  // Initiate generator with 'seed'
  init_MT19937(i_seed);
  
  // Set the number of spheres and dimres
  Ns = 4 * i_edge_fcc_N * i_edge_fcc_N * i_edge_fcc_N;
  Nd = Ns / 2;
  
  // Calculate cube edge
  cube_edge = i_edge_fcc_N * sqrt(two);
  
    // < DEBUG
  printf("Config variables:\n");
  printf("edge cells: %d\n",i_edge_fcc_N);
  printf("box edge:  %.16le\n",cube_edge);
  printf("half edge: %.16le\n",cube_edge/two);
  printf("coordinates range: %.16le %.16le\n",-cube_edge/two,cube_edge/two);
  printf("channel %d %d %d\n", i_channel[0], i_channel[1], i_channel[2]);
  printf("channel radius %.16le\n",i_channel_R);
  printf("Load initial DC from index: %d to index: %d (total %d files)\n",
          i_iDCfrom, i_iDCto, i_iDCto-i_iDCfrom+1);
  
  printf("Ns = %d, Nd = %d\n",Ns, Nd);
  // DEBUG >
  
  // Allocate memory for spheres and dimers
  spheres = malloc( Ns * sizeof(SPH));
  dimers  = malloc( Nd * sizeof(DIM3D));
  
  // Clean allocated memory
  memory_clean_spheres(spheres, Ns);
  memory_clean_dimers(dimers, Nd);
  
  
  
  // Loop over selected set of structures
  for(s=i_iDCfrom; s<=i_iDCto; s++){
    
    fprintf(stdout," Processing structure %d\n",s);
    // Regarding the ini settings, load input structure or grnerate a new one
    if(i_iDCfrom >= 0 && i_iDCto >= i_iDCfrom){
      sprintf(f_inp, fo_iDC, s);
      fprintf(stdout," Reading file %s\n",f_inp);
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
      // TODO: Set initial arrangemente into dimers   
      
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
    if(i_channel[0]+i_channel[1]+i_channel[2] > 0){
      fprintf(stdout, " Creating channel in the direction: %d %d %d\n",
        i_channel[0],i_channel[1],i_channel[2]);    
      make_channel(dimers, spheres, i_channel, i_channel_R, cube_edge, Nd);
    } 
    
    /* TODO: 
     * Brake dimers with spheres of type != 1
     * Set broken dimers as type '2'
     * set free (non-channel) spheres as type '3'
     */
    bd = brake_dimers(dimers, spheres, Nd);
    fprintf(stdout, " Dimers broken due to channel (if any): %d\n", bd);
    
    // Set the file name for dimer data and open the file for write
    sprintf(f_out, fo_exp_dim, s);
    fprintf(stdout," Writting dimer  data to file %s\n",f_out);
    
    if((f = fopen(f_out, "w")) == NULL) {
      fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
              prog_name, f_out);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }
    // Export dimer datat to file
    export_dimers(f, dimers, Nd);
    fclose(f);
    
    // Set the file name for dimer data and open the file for write
    sprintf(f_out, fo_exp_sph, s);
    fprintf(stdout," Writting sphere data to file %s\n",f_out);
        
    if((f = fopen(f_out, "w")) == NULL) {
      fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
              prog_name, f_out);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }
    // Export sphere datat to file
    export_spheres(f, spheres, Ns);
    fclose(f);

  #ifdef DATA_VISGL_OUTPUT
    // Export data in data_visGL format 
    if( export_to_GLviewer(spheres, cube_edge, s, Ns) != 0 ){
      exit_status=EXIT_FAILURE;
      goto cleanup;    
    }
  #endif
    
  } // End structure loop
 



  
cleanup:
  // Free resources
/*  
#ifdef DATA_VISGL_OUTPUT

  for(i=0; i<initial_N; i++){
    fclose(t[i]);
  }
  free(t);
#endif

#ifdef PUSH_MSD
  fclose(g);
#endif
*/
  free(spheres);
  free(dimers);
  return exit_status;
}


/* vim: set tw=80 ts=2 sw=2 et: */
