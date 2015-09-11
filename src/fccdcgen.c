/*
 * fccdcgen
 * f.c.c. structure of DC dimers generator/translator
 * Author: Jakub Narojczyk <narojczyk@ifmpan.poznan.pl>
 * Author: Mikolaj Kowalik <kowalik@ifmpan.poznan.pl>
 * 
 * The program generates a DC structure of dimers on f.c.c. spheres with the
 * presence of (optional) nano-chanel of spheres.
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

int main(int argc, char *argv[])
{

  FILE *f;
  
  SPH *spheres;
  DIM3D *dimers;

  char *f_ini = NULL;
  char f_out[15];
  char f_inp[15];
  
  const char *fo_iDC = "initdc3d.%03d";

  int exit_status;
  int Ns, Nd;
  
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
  printf("chanel %d %d %d\n", i_chanel[0], i_chanel[1], i_chanel[2]);
  printf("chanel radius %.16le\n",i_chanel_R);
  printf("Load initial DC from index: %d to index: %d (total %d files)\n",
          i_iDCfrom, i_iDCto, i_iDCto-i_iDCfrom+1);
  
  printf("Ns = %d, Nd = %d\n",Ns, Nd);
  // DEBUG >
  
  // Allocate memory for spheres and dimers
  spheres = malloc( Ns * sizeof(SPH));
  dimers  = malloc( Nd * sizeof(DIM3D));
  
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

    
    
//     for(i=0;i<Nd;i++){
//    printf("%3d (%4d %4d)\n",i,dimers[i].sph_ind[0], dimers[i].sph_ind[1]);
//   }
    
  }
 /*
  * obsolete debug code
  for(i=0;i<Ns;i++){
   printf("% lf % lf % lf\t %4d\n",spheres[i].r[0],spheres[i].r[1],spheres[i].r[2],i);
  }
*/
#ifdef DATA_VISGL_OUTPUT
  // Export data in data_visGL format 
  if( export_to_GLviewer(spheres, cube_edge, Ns) != 0 ){
    exit_status=EXIT_FAILURE;
    goto cleanup;    
  }
#endif
  
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
