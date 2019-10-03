/*
 * io.c
 *
 * File defines I/O functions used in the program
 */

#include <stdio.h>
#include <stdlib.h>
// #include <math.h>
// #include <string.h>
#include "config.h"
#include "data.h"
#include "io.h"

/*
 * export_structure_data()
 * Use export_spheres(), export_dimers(), and export_to_GLviewer() to write the
 * complete set of informations for the final structure.
 */
int exp_str_data(DIM3D *dim, SPH *sph, double box[3], int ns, int nd, int strn)
{
  FILE *file;
  char f_out[15];
  const char *fo_exp_dim = "d3d%05d.csv";
  const char *fo_exp_sph = "s3d%05d.csv";

  // Set the file name for dimer data and open the file for write
  sprintf(f_out, fo_exp_dim, strn);
  fprintf(stdout,"\n Writting dimer  data to file %s\n",f_out);

  if((file = fopen(f_out, "w")) == NULL) {
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__, f_out);
    return EXIT_FAILURE;
  }
  // Export dimer datat to file
  export_dimers(file, dim, nd);
  fclose(file);

  // Set the file name for dimer data and open the file for write
  sprintf(f_out, fo_exp_sph, strn);
  fprintf(stdout," Writting sphere data to file %s\n",f_out);

  if((file = fopen(f_out, "w")) == NULL) {
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__, f_out);
    return EXIT_FAILURE;
  }
  // Export sphere datat to file
  export_spheres(file, sph, ns);
  fclose(file);


#ifdef DATA_VISGL_OUTPUT
  // Export data in data_visGL format
  if( export_to_GLviewer(dim, sph, box, strn, ns, nd) != 0 ){
    return EXIT_FAILURE;
  }
#endif

  return 0;
}

/*
 * export_dimers()
 * Writes information about dimers to a file. Perform the write operation twice
 * (add more if required) in order to separate dimers of different type.
 */
int export_dimers(FILE *file, DIM3D *dim, int nd)
{
  const char *exp_f_dim =
  "%5d  %2d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %5d %5d\n";
  int i,k=0;

  for(i=0; i<nd; i++){
    // Export information about molecule
    if(dim[i].type == 1){
      if(fprintf(file, exp_f_dim,
          k++, dim[i].type,
          dim[i].R[0], dim[i].R[1], dim[i].R[2],
          dim[i].O[0], dim[i].O[1], dim[i].O[2],
          dim[i].L,
          dim[i].sph_ind[0], dim[i].sph_ind[1]) == EOF){
        fprintf(stderr,"  [%s]: error: exporting dimer data failed\n",__func__);
        return EXIT_FAILURE;
      }
    }
  }

  for(i=0; i<nd; i++){
    // Export information about molecule
    if(dim[i].type == 2){
      if(fprintf(file, exp_f_dim,
          k++, dim[i].type,
          dim[i].R[0], dim[i].R[1], dim[i].R[2],
          dim[i].O[0], dim[i].O[1], dim[i].O[2],
          dim[i].L,
          dim[i].sph_ind[0], dim[i].sph_ind[1]) == EOF){
        fprintf(stderr,"  [%s]: error: exporting dimer data failed\n",__func__);
        return EXIT_FAILURE;
      }
    }
  }
  return 0;
}

/*
 * export_spheres(f,sp,ns)
 * Export positions and other data for spheres
 */
int export_spheres(FILE *file, SPH *sph, int ns)
{
  const char *exp_f_sph_0 = "%5d %2d % .16le % .16le % .16le %.16le ";
  const char *exp_f_sph_1 = "%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d  ";
  const char *exp_f_sph_2 = "%2d %2d %2d\n";
  int i;

  for(i=0; i<ns; i++){
    // Write sphere positions and properties
    if(fprintf(file, exp_f_sph_0,
        i, sph[i].type,
        sph[i].r[0], sph[i].r[1], sph[i].r[2], sph[i].d) == EOF){
      fprintf(stderr,"  [%s]: error: exporting sphere data-1 failed\n", __func__);
      return EXIT_FAILURE;
    }

    // Write sphere neighbors table
    if(fprintf(file, exp_f_sph_1,
        sph[i].ngb[0], sph[i].ngb[1], sph[i].ngb[2], sph[i].ngb[3],
        sph[i].ngb[4], sph[i].ngb[5], sph[i].ngb[6], sph[i].ngb[7],
        sph[i].ngb[8], sph[i].ngb[9], sph[i].ngb[10], sph[i].ngb[11]) == EOF){
      fprintf(stderr,"  [%s]: error: exporting sphere data-2 failed\n", __func__);
      return EXIT_FAILURE;
    }

    // Write sphere lattice indexes
    if(fprintf(file, exp_f_sph_2,
      sph[i].lattice_ind[0], sph[i].lattice_ind[1], sph[i].lattice_ind[2]) == EOF){
      fprintf(stderr,"  [%s]: error: exporting sphere data-3 failed\n", __func__);
    return EXIT_FAILURE;
      }
  }
  return 0;
}

/*
 * export_to_GLviewer()
 * Export required data for viewing the structure by data_visGL
 */
int export_to_GLviewer(DIM3D *dim, SPH *sph, double box[3], int strn, int ns,
                       int nd)
{
  FILE *f;
  int i;
  char gl_f_out[16];
  const char *exp_form = "%5d % .16le % .16le % .16le %.16le %d\n";
  const char *exp_form_d =
    "%5d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %d\n";
  const char *data_visGL="vcontrol.ini";
  const char *data_spheresGL="gl_s3d%05d.csv";
  const char *data_dimersGL="gl_d3d%05d.csv";

  // Prepare control file for data_visGL program
  if((f = fopen(data_visGL, "w")) == NULL) {
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__,data_visGL);
    return EXIT_FAILURE;
  }

  // Write to data_visGL
  fprintf(f, "Type of data            : dim3dDCS\n");
  fprintf(f, "Number of input files   : %d\n", 2);
  fprintf(f, "Input lines per file    : %d\n", ns);
  fprintf(f, "Structure index         : %d\n", strn);
  fprintf(f, "(# not used          #) : %d\n", 0);
  fprintf(f, "Box dimensions          : %.16le %.16le %.16le\n",
          box[0], box[1], box[2]);
  fprintf(f, "(# not used          #) : %.16le\n", 0e0);
  fprintf(f, "(# not used          #) : %.16le\n", 0e0);

  // Close data_visGL
  fclose(f);

  // Open file for spheres data
  sprintf(gl_f_out, data_spheresGL, strn);
  fprintf(stdout," Writting sphere data to file %s (visualization)\n",gl_f_out);
  if((f = fopen(gl_f_out, "w")) == NULL){
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__, gl_f_out);
    return EXIT_FAILURE;
  }

  // Export structure data to data_spheresGL
  for(i=0; i<ns; i++){
    if(fprintf(f, exp_form, i,
        sph[i].r[0], sph[i].r[1], sph[i].r[2],
        sph[i].d, sph[i].type ) == EOF){
      fprintf(stderr,"  [%s]: error: exporting positions failed\n", __func__);
      return EXIT_FAILURE;
    }
  }

  // Close data_spheresGL
  fclose(f);

  // Open file for dimer data
  sprintf(gl_f_out, data_dimersGL, strn);
  fprintf(stdout," Writting dimer  data to file %s (visualization)\n",gl_f_out);
  if((f = fopen(gl_f_out, "w")) == NULL){
    fprintf(stderr, "  [%s]: error: cannot open config file %s\n",
            __func__, gl_f_out);
    return EXIT_FAILURE;
  }

  // Export structure data to data_dimersGL
  for(i=0; i<nd; i++){
    if(fprintf(f, exp_form_d, i,
        dim[i].R[0], dim[i].R[1], dim[i].R[2],
        dim[i].O[0], dim[i].O[1], dim[i].O[2],
        dim[i].L, dim[i].type ) == EOF){
      fprintf(stderr,"  [%s]: error: exporting positions failed\n", __func__);
      return EXIT_FAILURE;
    }
  }

  // Close data_dimersGL
  fclose(f);

  return 0;
}
