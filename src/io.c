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

extern const double two;
extern const double one;

void display_configuration_summary(CONFIG cfg, MODEL md, SLI *slits,
                                   CHA *channels){
  int i;
  const char *fmt_sees = " %-23s: %.16le %.16le (%1s)\n";
  const char *fmt_ses  = " %-23s: %.16le (%1s)\n";
  const char *fmt_sd   = " Max number of %-9s: %d\n";
  const char *fmt_ssDat   = " %-10s data from   : %s\n";
  const char *fmt_ss   = " %-23s: %s\n";
  const char *fmt_slu  = " %-23s: %8lu\n";
  const char *fmt_inlusion_data_header =
    " No.| %s offset\t\t  | %s normal\t       | %s\t  | %s\n";
  const char *fmt_inclusion_data =
    " %2d | %lf %lf % lf | %lf %lf %lf | %lf | %lf\n";

  fprintf(stdout," # %s - configuration summary\n", prog_name);
  fprintf(stdout,"\n ## Ini file parameters\n");
  fprintf(stdout," System size (cells)    : %d by %d by %d\n",
          cfg.cells[0], cfg.cells[1], cfg.cells[2]);
  fprintf(stdout, fmt_ss, "Symmetry", cfg.symmetry);
  fprintf(stdout," Structure indices      : from %d to %d (%d files total)\n",
          cfg.first, cfg.last, cfg.last-cfg.first+1);
  fprintf(stdout, fmt_sd, "channels", cfg.num_channels);
  fprintf(stdout, fmt_sd, "layers", cfg.num_slits);
  fprintf(stdout, fmt_slu, "PRNG seed", cfg.seed);

  fprintf(stdout,"\n ## Derived parameters\n");
  fprintf(stdout, fmt_sees, "Coordinates ranges",
          -md.box[0]/two, md.box[0]/two, "x");
  fprintf(stdout, fmt_sees, "", -md.box[1]/two,md.box[1]/two, "y");
  fprintf(stdout, fmt_sees, "", -md.box[2]/two,md.box[2]/two, "z");
  fprintf(stdout, fmt_ses, "Box dimensions", md.box[0], "x");
  fprintf(stdout, fmt_ses, "", md.box[1], "y");
  fprintf(stdout, fmt_ses, "", md.box[2], "z");
  fprintf(stdout, fmt_sd, "dimers", md.Nsph);
  fprintf(stdout, fmt_sd, "spheres", md.Ndim);

  fprintf(stdout,"\n ## Inclusions settings\n");
  if (cfg.mk_channel){
    fprintf(stdout, fmt_ssDat, "Channels",cfg.cfg_channels);
    fprintf(stdout, fmt_inlusion_data_header,
            "channel", "channel", "radius", "sph. diam.");
    for(i=0; i<cfg.num_channels; i++){
      fprintf(stdout, fmt_inclusion_data, i,
              channels[i].os[0], channels[i].os[1], channels[i].os[2],
              channels[i].nm[0], channels[i].nm[1], channels[i].nm[2],
              channels[i].radius, channels[i].sph_d);
    }
    fprintf(stdout,"\n");
  }

  if (cfg.mk_slit){
    fprintf(stdout, fmt_ssDat, "Layers", cfg.cfg_slits);
    fprintf(stdout, fmt_inlusion_data_header,
            "layer", "layer", "thick.", "sph. diam.");
    for(i=0; i<cfg.num_slits; i++){
      fprintf(stdout, fmt_inclusion_data, i,
              slits[i].os[0], slits[i].os[1], slits[i].os[2],
              slits[i].nm[0], slits[i].nm[1], slits[i].nm[2],
              slits[i].thickness, slits[i].sph_d);
    }
    fprintf(stdout,"\n");
  }
}

/*
 * display_stats()
 */
void display_stats(MODEL md)
{
  const char *fmt_sde = " %-17s (if any): %4d %6.2lf %%\n";

  double fNsph = one * md.Nsph;
  double hdim = (double) (200 * md.mtrx_dim); // (hekto dimers x2)
  double hsph = (double) (100 * md.mtrx_sph); // (hekto spheres)
  double hisp = (double) (100 * md.incl_sph); // (hekto inc. spheres)

  fprintf(stdout, fmt_sde, "Free spheres",      md.mtrx_sph, hsph/fNsph);
  fprintf(stdout, fmt_sde, "Inclusion spheres", md.incl_sph, hisp/fNsph);
  fprintf(stdout, fmt_sde, "Dimers valid",      md.mtrx_dim, hdim/fNsph);
}

/*
 * export_structure_data()
 * Use export_spheres(), export_dimers(), and export_to_GLviewer() to write the
 * complete set of informations for the final structure.
 */
int exp_str_data(CONFIG cf, MODEL md, DIM3D *dim, SPH *sph,
                 int ns3, int nd2, int strn)
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
  fprintf(file,"Number of spheres       : %d\n", md.Nsph);
  // TODO: this value is now erronusly calculated
  fprintf(file,"Number of particles     : %d\n", md.Ndim + nd2);
  // TODO: the 2*nd gives wrong value for the number of spheres
  fprintf(file,"Number of x-mers        : %d %d\n",2*nd2, md.Ndim - nd2);
  fprintf(file,"Box matrix h00 h01 h02  : %.16le %.16le %.16le\n",
          md.box[0], 0e0, 0e0);
  fprintf(file,"Box matrix     h11 h12  : %.16le %.16le\n",
          md.box[1], 0e0);
  fprintf(file,"Box matrix         h22  : %.16le\n", md.box[2]);
  fclose(file);
  
  // Set the file name for sphere data and open the file for write
  sprintf(f_out, fmt_exp_sph, strn);
  fprintf(stdout, fmt_message, "sphere", f_out);
  
  if((file = fopen(f_out, "w")) == NULL) {
    fprintf(stderr, fmt_exporting_failed, __func__, f_out);
    return EXIT_FAILURE;
  }
  // Exports complete Nsph sphere data regardles of particles they form
  // (additional constraints and bonds are exported to files for the given
  // molecules)
  export_spheres(file, sph, md.Nsph);
  fclose(file);

  if(cf.mk_dimers){
    // Set the file name for dimer data and open the file for write
    sprintf(f_out, fmt_exp_dim, strn);
    fprintf(stdout, fmt_message, "dimer", f_out);

    if((file = fopen(f_out, "w")) == NULL) {
      fprintf(stderr, fmt_exporting_failed, __func__, f_out);
      return EXIT_FAILURE;
    }
    // Export dimer datat to file
    export_dimers(file, dim, md.Ndim);
    fclose(file);
  }

#ifdef DATA_VISGL_OUTPUT
  // Export data in data_visGL format
  if( export_to_GLviewer(md, dim, sph, strn) != 0 ){
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
  extern const int TYPE_DIMER;
  const char *fmt_exp_dim =
    "%5d  %3d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %5d %5d\n";
  const char *fmt_exporting_failed =
    "  [%s] ERR: exporting dimer data failed\n";
  int i, k=0;

  for(i=0; i<nd; i++){
    // Export information about molecule
    if(dim[i].type == TYPE_DIMER){
      if(fprintf(file, fmt_exp_dim,
          k++, dim[i].type,
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
  for(i=0; i<ns; i++){
    // Write sphere positions and properties
    if(fprintf(file, fmt_exp_sph_0,
        i, sph[i].type,
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
 * double box[3],  int ns,  int nd
 */
int export_to_GLviewer(MODEL md, DIM3D *dim, SPH *sph, int strn){
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
  fprintf(file, "Input lines per file    : %d\n", md.Nsph);
  fprintf(file, "Structure index         : %d\n", strn);
  fprintf(file, "(# not used          #) : %d\n", 0);
  fprintf(file, "Box dimensions          : %.16le %.16le %.16le\n",
          md.box[0], md.box[1], md.box[2]);
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
  for(i=0; i<md.Nsph; i++){
    if(fprintf(file, exp_form, i,
        sph[i].r[0], sph[i].r[1], sph[i].r[2],
        sph[i].d,
        legacy_GLexport_sphere_type_converter(sph[i].type) ) == EOF){
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
  for(i=0; i<md.Ndim; i++){
    if(fprintf(file, exp_form_d, i,
        dim[i].R[0], dim[i].R[1], dim[i].R[2],
        dim[i].O[0], dim[i].O[1], dim[i].O[2],
        dim[i].L,
        legacy_GLexport_dimer_type_converter(dim[i].type) ) == EOF){
      fprintf(stderr, fmt_exporting_failed, "dimer", __func__);
      return EXIT_FAILURE;
    }
  }
  fclose(file);

  return 0;
}

int legacy_GLexport_dimer_type_converter(int type){
  extern const int TYPE_DIMER;
  if(type == TYPE_DIMER) return 1;
  return type;
}

int legacy_GLexport_sphere_type_converter(int type){
  extern const int TYPE_SPHERE;
  extern const int TYPE_SPHERE_DIMER;
  extern const int TYPE_INCLUSION_SPHERE;

  if(type == TYPE_SPHERE) return 3;
  if(type == TYPE_INCLUSION_SPHERE) return 2;
  if(type == TYPE_SPHERE_DIMER) return 1;
  return -1;
/*  // conversion legend - remove when verified that it works
    if(sph[i].type == 3){
      newtype = 1;
    }else if(sph[i].type == 2){
      newtype = 101;
    }else if(sph[i].type == 1){
      newtype = 2;
    }else{
      newtype = -1;
    }*/
}
