/*
 * io.c
 *
 * File defines I/O functions used in the program
 */

#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "data.h"
#include "io.h"

extern char *prog_name;
extern unsigned long int i_seed;
extern int i_edge_fcc_N[3];
extern int i_iDCfrom;
extern int i_iDCto;
extern int i_make_channel;
extern int i_make_slit;
extern int i_n_channels;
extern int i_n_slits;
extern char i_Fchannels[41];
extern char i_Fslits[41];
extern const double two;

void display_configuration_summary(SLI *slits, CHA *channels,
                                   double box_edge[3], int Ns, int Nd)
{
  int i;
  const char *fmt_sees = " %-21s: %.16le %.16le (%1s)\n";
  const char *fmt_ses  = " %-21s: %.16le (%1s)\n";
  const char *fmt_sd   = " Number of %-11s: %d\n";
  const char *fmt_ss   = " %-8s data from   : %s\n";
  const char *fmt_inlusion_data_header =
    " No.| %s offset\t\t  | %s normal\t       | %s\t  | %s\n";
  const char *fmt_inclusion_data =
    " %2d | %lf %lf % lf | %lf %lf %lf | %lf | %lf\n";

  fprintf(stdout," ### %s - configuration summary\n", prog_name);
  fprintf(stdout,"\n ## Ini file parameters\n");
  fprintf(stdout," System size (cells)  : %d by %d by %d\n",
          i_edge_fcc_N[0], i_edge_fcc_N[1], i_edge_fcc_N[2]);
  fprintf(stdout," PRNG seed            : %8lu\n", i_seed);
  fprintf(stdout," Structure indices    : from %d to %d (%d files total)\n",
          i_iDCfrom, i_iDCto, i_iDCto-i_iDCfrom+1);

  fprintf(stdout,"\n ## Derived parameters\n");
  fprintf(stdout, fmt_sees, "Coordinates ranges",
          -box_edge[0]/two,box_edge[0]/two, "x");
  fprintf(stdout, fmt_sees, "", -box_edge[1]/two,box_edge[1]/two, "y");
  fprintf(stdout, fmt_sees, "", -box_edge[2]/two,box_edge[2]/two, "z");
  fprintf(stdout, fmt_ses, "Box dimensions", box_edge[0], "x");
  fprintf(stdout, fmt_ses, "", box_edge[1], "y");
  fprintf(stdout, fmt_ses, "", box_edge[2], "z");
  fprintf(stdout, fmt_sd, "dimers", Nd);
  fprintf(stdout, fmt_sd, "spheres", Ns);
  fprintf(stdout, fmt_sd, "channels", i_n_channels);

  fprintf(stdout,"\n ## Inclusions settings\n");
  if (i_make_channel){
    fprintf(stdout, fmt_ss, "Channels",i_Fchannels);
    fprintf(stdout, fmt_inlusion_data_header,
            "channel", "channel", "radius", "sph. diam.");
    for(i=0; i<i_n_channels; i++){
      fprintf(stdout, fmt_inclusion_data, i,
              channels[i].os[0], channels[i].os[1], channels[i].os[2],
              channels[i].nm[0], channels[i].nm[1], channels[i].nm[2],
              channels[i].radius, channels[i].sph_d);
    }
    fprintf(stdout,"\n");
  }
  fprintf(stdout, fmt_sd, "layers",i_n_slits);
  if (i_make_slit){
    fprintf(stdout, fmt_ss, "Layers", i_Fslits);
    fprintf(stdout, fmt_inlusion_data_header,
            "layer", "layer", "thick.", "sph. diam.");
    for(i=0; i<i_n_slits; i++){
      fprintf(stdout, fmt_inclusion_data, i,
              slits[i].os[0], slits[i].os[1], slits[i].os[2],
              slits[i].nm[0], slits[i].nm[1], slits[i].nm[2],
              slits[i].thickness, slits[i].sph_d);
    }
    fprintf(stdout,"\n");
  }
}

/*
void print_inclusion_data(const char *fmt_ss, const char *fmt_inclusion_data,
                          const char *inc_header,
                          const char *inc_type, const char *inc_size,
                          char fileName[], int n_inclusions
                          double os[3], double nm[3], double thickness,
                          double sph_d)
{
  int i;
  fprintf(stdout, fmt_ss, inc_header, fileName);
    fprintf(stdout, fmt_inlusion_data_header,
            inc_type, inc_type, inc_size, "sph. diam.");
    for(i=0; i<n_inclusions; i++){
      fprintf(stdout, fmt_inclusion_data, i,
              os[0], os[1], os[2],
              nm[0], nm[1], nm[2],
              thickness, sph_d);
    }
    fprintf(stdout,"\n");
}*/

/*
 * export_structure_data()
 * Use export_spheres(), export_dimers(), and export_to_GLviewer() to write the
 * complete set of informations for the final structure.
 */
int exp_str_data(DIM3D *dim, SPH *sph, double box[3], 
                 int ns, int nd, int ns3, int nd2, int strn)
{
  FILE *file;
  char f_out[128];
  const char *fmt_exp_dsc = "s3d0_summary_%05d.ini";
  const char *fmt_exp_sph = "s3d1_monomers_%05d.csv";
  const char *fmt_exp_dim = "s3d2_dimers_%05d.csv";
  const char *fmt_exporting_failed = "  [%s] ERR: writting to %s failed\n";
  const char *fmt_message = " Writting %-6s data to %s\n";

  // Set the file name for info data and open the file for write
  sprintf(f_out, fmt_exp_dsc, strn);
  fprintf(stdout, fmt_message, "ini", f_out);

  if((file = fopen(f_out, "w")) == NULL) {
    fprintf(stderr, fmt_exporting_failed, __func__, f_out);
    return EXIT_FAILURE;
  }
  // Export str. descrip. data to file
  fprintf(file,"Number of spheres       : %d\n", ns);
  // TODO: this is valid only for dimers
  fprintf(file,"Number of particles     : %d\n", nd + nd2); 
  fprintf(file,"Number of x-mers        : %d %d\n",2*nd2, nd - nd2);
  fprintf(file,"Box matrix h00 h01 h02  : %.16le %.16le %.16le\n",
          box[0], 0e0, 0e0);
  fprintf(file,"Box matrix     h11 h12  : %.16le %.16le\n",
          box[1], 0e0);
  fprintf(file,"Box matrix         h22  : %.16le\n", box[2]);
  fclose(file);
  
  // Set the file name for dimer data and open the file for write
  sprintf(f_out, fmt_exp_sph, strn);
  fprintf(stdout, fmt_message, "sphere", f_out);
  
  if((file = fopen(f_out, "w")) == NULL) {
    fprintf(stderr, fmt_exporting_failed, __func__, f_out);
    return EXIT_FAILURE;
  }
  // Export sphere datat to file
  export_spheres(file, sph, ns);
  fclose(file);

  // Set the file name for dimer data and open the file for write
  sprintf(f_out, fmt_exp_dim, strn);
  fprintf(stdout, fmt_message, "dimer", f_out);

  if((file = fopen(f_out, "w")) == NULL) {
    fprintf(stderr, fmt_exporting_failed, __func__, f_out);
    return EXIT_FAILURE;
  }
  // Export dimer datat to file
  export_dimers(file, dim, nd);
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
  const char *fmt_exp_dim =
    "%5d  %3d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %5d %5d\n";
  const char *fmt_exporting_failed =
    "  [%s] ERR: exporting dimer data failed\n";
  int i, k=0;

  for(i=0; i<nd; i++){
    // Export information about molecule
    if(dim[i].type == 1){
      if(fprintf(file, fmt_exp_dim,
          k++, 2,
          dim[i].R[0], dim[i].R[1], dim[i].R[2],
          dim[i].O[0], dim[i].O[1], dim[i].O[2],
          dim[i].L,
          dim[i].sph_ind[0], dim[i].sph_ind[1]) == EOF){
        fprintf(stderr, fmt_exporting_failed, __func__);
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
  const char *fmt_exp_sph_0 = "%5d %3d % .16le % .16le % .16le %.16le ";
  const char *fmt_exp_sph_1 =
    "%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d  ";
  const char *fmt_exp_sph_2 = "  %2d %2d %2d\n";
  const char *fmt_exporting_failed = "  [%s] ERR: exporting sphere %s failed\n";
  int i;
  int newtype;
  for(i=0; i<ns; i++){
    // Modify type for new program
    if(sph[i].type == 3){
      newtype = 1;
    }else if(sph[i].type == 2){
      newtype = 101;
    }else if(sph[i].type == 1){
      newtype = 2;
    }else{
      newtype = -1;
    }
    // Write sphere positions and properties
    if(fprintf(file, fmt_exp_sph_0,
        i, newtype,
        sph[i].r[0], sph[i].r[1], sph[i].r[2], sph[i].d) == EOF){
      fprintf(stderr, fmt_exporting_failed, "positions", __func__);
      return EXIT_FAILURE;
    }

    // Write sphere neighbors table
    if(fprintf(file, fmt_exp_sph_1,
        sph[i].ngb[0], sph[i].ngb[1], sph[i].ngb[2], sph[i].ngb[3],
        sph[i].ngb[4], sph[i].ngb[5], sph[i].ngb[6], sph[i].ngb[7],
        sph[i].ngb[8], sph[i].ngb[9], sph[i].ngb[10], sph[i].ngb[11]) == EOF){
      fprintf(stderr, fmt_exporting_failed, "neighbors", __func__);
      return EXIT_FAILURE;
    }

    // Write sphere lattice indexes
    if(fprintf(file, fmt_exp_sph_2,
      sph[i].lattice_ind[0], sph[i].lattice_ind[1],
      sph[i].lattice_ind[2]) == EOF){
        fprintf(stderr, fmt_exporting_failed, "lattice idnices", __func__);
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
  FILE *file;
  int i;
  char f_GLout[64];
  const char *exp_form = "%5d % .16le % .16le % .16le %.16le %d\n";
  const char *exp_form_d =
    "%5d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %d\n";
  const char *fmt_opening_failed = "  [%s] ERR: cannot open file %s\n";
  const char *fmt_message = " Writting %-6s data to file %s (visualization)\n";
  const char *fmt_exporting_failed = "  [%s] ERR: exporting %s data failed\n";
  const char *data_visGL = "vcontrol.ini";
  const char *data_spheresGL = "gl_s3d%05d.csv";
  const char *data_dimersGL = "gl_d3d%05d.csv";

  // Prepare control file for data_visGL program
  if((file = fopen(data_visGL, "w")) == NULL) {
    fprintf(stderr, fmt_opening_failed, __func__,data_visGL);
    return EXIT_FAILURE;
  }

  // Write to data_visGL
  fprintf(file, "Type of data            : dim3dDCS\n");
  fprintf(file, "Number of input files   : %d\n", 2);
  fprintf(file, "Input lines per file    : %d\n", ns);
  fprintf(file, "Structure index         : %d\n", strn);
  fprintf(file, "(# not used          #) : %d\n", 0);
  fprintf(file, "Box dimensions          : %.16le %.16le %.16le\n",
          box[0], box[1], box[2]);
  fprintf(file, "(# not used          #) : %.16le\n", 0e0);
  fprintf(file, "(# not used          #) : %.16le\n", 0e0);
  fclose(file);

  // Open file for spheres data
  sprintf(f_GLout, data_spheresGL, strn);
  fprintf(stdout, fmt_message, "sphere", f_GLout);
  if((file = fopen(f_GLout, "w")) == NULL){
    fprintf(stderr, fmt_opening_failed, __func__, f_GLout);
    return EXIT_FAILURE;
  }

  // Export structure data to data_spheresGL
  for(i=0; i<ns; i++){
    if(fprintf(file, exp_form, i,
        sph[i].r[0], sph[i].r[1], sph[i].r[2],
        sph[i].d, sph[i].type ) == EOF){
      fprintf(stderr, fmt_exporting_failed, "sphere", __func__);
      return EXIT_FAILURE;
    }
  }
  fclose(file);

  // Open file for dimer data
  sprintf(f_GLout, data_dimersGL, strn);
  fprintf(stdout, fmt_message, "dimer",f_GLout);
  if((file = fopen(f_GLout, "w")) == NULL){
    fprintf(stderr, fmt_opening_failed, __func__, f_GLout);
    return EXIT_FAILURE;
  }

  // Export structure data to data_dimersGL
  for(i=0; i<nd; i++){
    if(fprintf(file, exp_form_d, i,
        dim[i].R[0], dim[i].R[1], dim[i].R[2],
        dim[i].O[0], dim[i].O[1], dim[i].O[2],
        dim[i].L, dim[i].type ) == EOF){
      fprintf(stderr, fmt_exporting_failed, "dimer", __func__);
      return EXIT_FAILURE;
    }
  }
  fclose(file);

  return 0;
}
