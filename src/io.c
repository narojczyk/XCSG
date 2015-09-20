/*
 * io.c
 *
 * File defines I/O functions used in the program
 */

#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
// #include <string.h>
// #include "config.h"
#include "data.h"
#include "io.h"

/*
 * expot_structure()
 * Writes complete information about the structure to a file.
 */
int export_structure(FILE *file, DIM3D *dim, SPH *sph, int nd)
{
  const char *exp_f_dim = 
  "%5d  %2d % .16le % .16le % .16le % .16le % .16le % .16le %.16le, %5d %5d\n";
//   const char *exp_f_dim_ngb = 
  const char *exp_f_sph = "%5d % .16le % .16le % .16le %.16le %d\n";
  int i, atom0, atom1;

  for(i=0; i<nd; i++){
    // Export information about molecule
    atom0 = dim[i].sph_ind[0];
    atom1 = dim[i].sph_ind[1];
    if(fprintf(file, exp_f_dim, 
        i, dim[i].type,
        dim[i].R[0], dim[i].R[1], dim[i].R[2], 
        dim[i].O[0], dim[i].O[1], dim[i].O[2],
        dim[i].L,
        atom0, atom1) == EOF){
      fprintf(stderr,"  [%s]: error: exporting positions failed\n", __func__);
      return EXIT_FAILURE;
    }
    /*
    if(fprintf(file, exp_f_dim, i, 
        sph[i].r[0], sph[i].r[1], sph[i].r[2], 
        sph[i].d, sph[i].type ) == EOF){
      fprintf(stderr,"  [%s]: error: exporting positions failed\n", __func__);
      return EXIT_FAILURE;
    }*/
  }
  return 0;
}

/*
 * load_dcsgen()
 * Load initial DC configuration from dcsgen, program by Mikolaj Kowalik
 */
int load_dcsgen(FILE *file, DIM3D *dim, double box_x, int nd)
{
  // Template structure instance for data input
  DIM3D tm; 
  int i, j;
  
  // Initiate fields not read from file with defaults
  tm.type = 1;
  tm.sph_ind[0] = tm.sph_ind[1] = -1;
  for(i=0; i<22; i++){
    tm.ngb[i][0] = -1;
    tm.ngb[i][1] = -1;
  }
  tm.L = 1e0;

  while(fscanf(file, "%d %*d %*d %lf %lf %lf %lf %lf %lf",
               &i, tm.R, tm.R + 1, tm.R + 2, tm.O, tm.O + 1, tm.O + 2) != EOF) {
    // Decrement 'i' to count from 0
    if(--i < nd){
      for(j=0; j<3; j++){
        // Scale molecule positions to cube dimensions
        tm.R[j] *= box_x;
      }
      // Copy data from template to data structure
      dim[i] = tm;
    }
  }

  // Print errors when input structure size is different than ini parameters
  if(++i != nd) {
    fprintf(stderr, "  [%s]: error: input structure ", __func__);
    if(i > nd) {
      fprintf(stderr, "too large;\n");
    } else if(i < nd) {
    fprintf(stderr, "too small;\n");
    }

    fprintf(stderr,
      "  [%s]: error: input file contains %d molecules,\n", __func__, i);
    fprintf(stderr,
      "  [%s]: error: allocated memory for %d molecules\n", __func__, nd);
  }
  
  // Return the number of data lines imported
  return i;
}

/*
 * export_spheres(f,sp,ns)
 * Export positions and other data for spheres
 */
int export_spheres(FILE *file, SPH *sph, int ns)
{
  const char *exp_form = "%5d % .16le % .16le % .16le %.16le %d\n";
  int i;

  for(i=0; i<ns; i++){
    if(fprintf(file, exp_form, i, 
        sph[i].r[0], sph[i].r[1], sph[i].r[2], 
        sph[i].d, sph[i].type ) == EOF){
      fprintf(stderr,"  [%s]: error: exporting positions failed\n", __func__);
      return EXIT_FAILURE;
    }
  }
  return 0;
}

/*
 * export_to_GLviewer()
 * Export required data for viewing the structure by data_visGL
 */
int export_to_GLviewer(SPH *sph, double box_x, int ns)
{
  FILE *f;
  const char* data_visGL="vcontrol.ini";
  const char* data_spheresGL="gl_spheres.csv";
  
  // Prepare control file for data_visGL program
  if((f = fopen(data_visGL, "w")) == NULL) {
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__,data_visGL);
    return EXIT_FAILURE;
  }

  // Write to data_visGL
  fprintf(f, "Type of data            : dim3dDCS\n");
  fprintf(f, "Number of input files   : %d\n", 1);
  fprintf(f, "Input lines per file    : %d\n", ns);
  fprintf(f, "Number of searchers     : %d\n", 0);
  fprintf(f, "Number of targets       : %d\n", 0);
  fprintf(f, "Box dimensions          : %.16le %.16le\n", box_x, box_x);
  fprintf(f, "Searcher radius         : %.16le\n", 0e0);
  fprintf(f, "Target radius           : %.16le\n", 0e0);
  
  // Close data_visGL
  fclose(f);
  
  // Open file for spheres data
  if((f = fopen(data_spheresGL, "w")) == NULL){
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__, data_spheresGL);
    return EXIT_FAILURE;
  }
  
  // Export structure data to data_spheresGL
  if( export_spheres(f, sph, ns) != 0 ){
    return EXIT_FAILURE;
  }

  // Close data_spheresGL
  fclose(f);
  
  return 0;
}