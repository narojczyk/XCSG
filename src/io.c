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

static int export_dimers(FILE *file, DIM3D *dim, int nd);
static int export_spheres(FILE *file, SPH *sph, int ns);
static int povray_export_spheres(FILE *file, SPH *sph, int ns);
static int povray_export_dimers(FILE *file, DIM3D *dim, int nd);
static int export_to_GLviewer(MODEL md, DIM3D *dim, SPH *sph, int strn);
static int legacy_GLexport_dimer_type_converter(int type);
static int legacy_GLexport_sphere_type_converter(int type);
// These function can close program on error
static FILE* open_file(const char *file, const char *mode, int strict);

/* # SEC ############## CONSOLE OUTPUT ###################################### */

void show_cfg_summary(CONFIG cfg, MODEL md, INC *slits, INC *channels,
                      INC *clusters){
  INC *cl = NULL, *ch = NULL, *la = NULL;
  int i;
  const char *fmt_sees = " %-23s: %.16le %.16le (%1s)\n";
  const char *fmt_ses  = " %-23s: %.16le (%1s)\n";
  const char *fmt_ssDat   = " %-10s data from   : %s\n";
  const char *fmt_ss   = " %-23s: %s\n";
  const char *fmt_slu  = " %-23s: %8lu\n";
  const char *fmt_inlusion_data_header =
    " No.| %s offset\t\t  | %s normal\t       | %s\t  | %s | %s\n";
  const char *fmt_inclusion_data =
    " %2d | %lf %lf % lf | %lf %lf %lf | %lf | %lf | %d \n";

  fprintf(stdout," # %s - configuration summary\n", prog_name);
  fprintf(stdout,"\n ## Ini file parameters\n");
  fprintf(stdout," System size (cells)    : %d by %d by %d (%d spheres)\n",
          cfg.cells[0], cfg.cells[1], cfg.cells[2], md.Nsph);
  fprintf(stdout, fmt_ss, "Symmetry", cfg.symmetry);
  fprintf(stdout," Structure indices      : from %d to %d (%d files total)\n",
          cfg.first, cfg.last, cfg.last-cfg.first+1);
  fprintf(stdout, fmt_ss, "Clusters file", cfg.cfg_clusters);
  fprintf(stdout, fmt_ss, "Channels file", cfg.cfg_channels);
  fprintf(stdout, fmt_ss, "Layers file", cfg.cfg_slits);
  fprintf(stdout, fmt_slu, "PRNG seed", cfg.seed);

  fprintf(stdout,"\n ## Derived parameters\n");
  fprintf(stdout, fmt_sees, "Coordinates ranges",
          -md.box[0]/two, md.box[0]/two, "x");
  fprintf(stdout, fmt_sees, "", -md.box[1]/two,md.box[1]/two, "y");
  fprintf(stdout, fmt_sees, "", -md.box[2]/two,md.box[2]/two, "z");
  fprintf(stdout, fmt_ses, "Box dimensions", md.box[0], "x");
  fprintf(stdout, fmt_ses, "", md.box[1], "y");
  fprintf(stdout, fmt_ses, "", md.box[2], "z");

  fprintf(stdout,"\n ## Inclusions settings\n");
  if (cfg.mk_clusters){
    fprintf(stdout, fmt_ssDat, "Cluster inclusions",cfg.cfg_clusters);
    fprintf(stdout, fmt_inlusion_data_header,
            "cluster", "(ignored)", "radius", "sp. dia.", "n-mer");
    for(i=0; i<cfg.num_clusters; i++){
      cl = &clusters[i];
      fprintf(stdout, fmt_inclusion_data, i,
              (*cl).os[0],  (*cl).os[1], (*cl).os[2],
              (*cl).nm[0],  (*cl).nm[1], (*cl).nm[2],
              (*cl).radius, (*cl).sph_d, (*cl).tgt_Nmer);
    }
    fprintf(stdout,"\n");
  }

  if (cfg.mk_channel){
    fprintf(stdout, fmt_ssDat, "Channels",cfg.cfg_channels);
    fprintf(stdout, fmt_inlusion_data_header,
            "channel", "channel", "radius", "sp. dia.", "n-mer");
    for(i=0; i<cfg.num_channels; i++){
      ch = &channels[i];
      fprintf(stdout, fmt_inclusion_data, i,
              (*ch).os[0],  (*ch).os[1], (*ch).os[2],
              (*ch).nm[0],  (*ch).nm[1], (*ch).nm[2],
              (*ch).radius, (*ch).sph_d, (*ch).tgt_Nmer);
    }
    fprintf(stdout,"\n");
  }

  if (cfg.mk_slit){
    fprintf(stdout, fmt_ssDat, "Layers", cfg.cfg_slits);
    fprintf(stdout, fmt_inlusion_data_header,
            "layer", "layer", "thick.", "sp. dia.", "n-mer");
    for(i=0; i<cfg.num_slits; i++){
      la = &slits[i];
      fprintf(stdout, fmt_inclusion_data, i,
             (*la).os[0],     (*la).os[1],  (*la).os[2],
             (*la).nm[0],     (*la).nm[1],  (*la).nm[2],
             (*la).thickness, (*la).sph_d,  (*la).tgt_Nmer);
    }
    fprintf(stdout,"\n");
  }
}

/*
 * show_particle_stats()
 */
void show_particle_stats(MODEL md, CONFIG cfg)
{
  const char *fmt_sde = " %-17s: %4d %6.2lf%%\n";

  double fNsph = one * md.Nsph;
  double hdim = (double) (200 * md.mtrx_dim); // (hekto dimers x2)
  double hsph = (double) (100 * md.mtrx_sph); // (hekto spheres)
  double hidm = (double) (200 * md.incl_dim); // (hekto inc. dimers x2)
  double hisp = (double) (100 * md.incl_sph); // (hekto inc. spheres)

  fprintf(stdout, "\n");
  fprintf(stdout, fmt_sde, "Matrix spheres",    md.mtrx_sph, hsph/fNsph);
  if(cfg.mk_dimers){
    fprintf(stdout, fmt_sde, "Matrix dimers",      md.mtrx_dim, hdim/fNsph);
  }
  if(cfg.mk_clusters | cfg.mk_channel | cfg.mk_slit){
    fprintf(stdout, fmt_sde, "Inclusion spheres", md.incl_sph, hisp/fNsph);
  }
  if(md.incl_dim){
    fprintf(stdout, fmt_sde, "Inclusion dimers", md.incl_dim, hidm/fNsph);
  }
}

/* # SEC ############## FILE OUTPUT ######################################### */

/*
 * export_structure()
 * Use export_spheres(), export_dimers(), ... (possibly other -mers())
 * to write the complete set of informations for the final structure.
 * Use export_to_GLviewer() and povray_export_spheres() to export sphere data
 * in formats readable to graphic viewers. TODO: add povray_export_dimers()
 */
int export_structure(CONFIG cf, MODEL md, PARTICLES pts, int strn)
{
#ifdef DEBUG_MODE
  extern const char *fmt_dbg_opening_file;
#endif
  extern const char *fmt_write_notify;
  SPH *sph = pts.spheres;
  DIM3D *dim = pts.dimers;
  FILE *file = NULL;
  char f_out[128];
  const char *fmt_exp_dsc = "s3d0_summary_%05d.ini";
  const char *fmt_exp_sph = "s3d1_monomers_%05d.csv";
  const char *fmt_exp_sph_pov = "v_s3d1_monomers_%05d.pov";
  const char *fmt_exp_dim = "s3d2_dimers_%05d.csv";
  const char *fmt_exp_dim_pov = "v_s3d2_dimers_%05d.pov";
  const char *fmt_particle_qty_d  = "%-24s: %5d\n";
  const char *fmt_particle_qty_dd = "%-24s: %5d %5d\n";

  // Set the file name for info data and open the file for write
  sprintf(f_out, fmt_exp_dsc, strn);
  fprintf(stdout, "\n");

#ifdef DEBUG_MODE
  fprintf(stdout, fmt_dbg_opening_file, __func__, f_out);
#endif
  file = open_to_write(f_out);

  fprintf(stdout, fmt_write_notify, "ini", f_out);
  // Export str. descrip. data to file
  fprintf(file, fmt_particle_qty_d, "Number of spheres", md.Nsph);
  fprintf(file, fmt_particle_qty_dd, "Matrix particles",
            md.mtrx_sph, md.mtrx_dim);
  fprintf(file, fmt_particle_qty_dd, "Inclusion particles",
            md.incl_sph, md.incl_dim);
  fprintf(file,"Box matrix h00 h01 h02  : %.16le %.16le %.16le\n",
          md.box[0], 0e0, 0e0);
  fprintf(file,"Box matrix     h11 h12  : %.16le %.16le\n",
          md.box[1], 0e0);
  fprintf(file,"Box matrix         h22  : %.16le\n", md.box[2]);
  fclose(file);
  
  // Set the file name for sphere data and open the file for write
  sprintf(f_out, fmt_exp_sph, strn);

#ifdef DEBUG_MODE
  fprintf(stdout, fmt_dbg_opening_file, __func__, f_out);
#endif

  // Exports complete Nsph sphere data regardles of particles they form
  // (additional constraints and bonds are exported to files for the given
  // molecules)
  file = open_to_write(f_out);
  fprintf(stdout, fmt_write_notify, "sphere", f_out);
  if(export_spheres(file, sph, md.Nsph) != EXIT_SUCCESS){
    return EXIT_FAILURE;
  }
  file = NULL;

  if(cf.mk_dimers || md.incl_dim){
    // Set the file name for dimer data and open the file for write
    sprintf(f_out, fmt_exp_dim, strn);
  #ifdef DEBUG_MODE
    fprintf(stdout, fmt_dbg_opening_file, __func__, f_out);
  #endif
    file = open_to_write(f_out);
    // Export dimer datat to file
    fprintf(stdout, fmt_write_notify, "dimer", f_out);
    if(export_dimers(file, dim, md.Ndim) != EXIT_SUCCESS){
      return EXIT_FAILURE;
    }
    file = NULL;
  }

#ifdef DATA_VISGL_OUTPUT
  // Export data in data_visGL format
  if( export_to_GLviewer(md, dim, sph, strn) != EXIT_SUCCESS ){
    return EXIT_FAILURE;
  }
#endif

  // Export data in POV-Ray format
  // export sphere data
  sprintf(f_out, fmt_exp_sph_pov, strn);
#ifdef DEBUG_MODE
  fprintf(stdout, fmt_dbg_opening_file, __func__, f_out);
#endif
  file = open_to_write(f_out);
  fprintf(stdout, fmt_write_notify, "sphere", f_out);
  if(povray_export_spheres(file, sph, md.Nsph) != EXIT_SUCCESS){
    return EXIT_FAILURE;
  }

  // export dimer data
  if(cf.mk_dimers || md.incl_dim){
    sprintf(f_out, fmt_exp_dim_pov, strn);
  #ifdef DEBUG_MODE
    fprintf(stdout, fmt_dbg_opening_file, __func__, f_out);
  #endif
    file = open_to_write(f_out);
    fprintf(stdout, fmt_write_notify, "dimer", f_out);
    if(povray_export_dimers(file, dim, md.incl_dim + md.mtrx_dim) != EXIT_SUCCESS){
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

/*
 * povray_export_spheres(file, sph, ns)
 * Export sphere data in format readable by POV-Ray
 */
static int povray_export_spheres(FILE *file, SPH *sph, int ns){
  extern const int TYPE_MATRIX_LIMIT;
  extern const int TYPE_INCLUSION_BASE;
  extern const char *fmt_writting_failed;
  const char *fmt_data_corruption =
    " [%s] ERR: incomplete sphere data (%d/%d) exported to file\n";
  const char *fmt_exp_sph_0 =
    " sph(<% .6le, % .6le, % .6le> + TL, %s, %s%04d) // sph %5d type %d\n";
  int i, exp_count = 0;

  fprintf(file, "#declare matrixSpheres = union{\n");

  for(i=0; i<ns; i++){
    // Write sphere positions and properties
    if(sph[i].type <= TYPE_MATRIX_LIMIT){
      if(fprintf(file, fmt_exp_sph_0, sph[i].r[0], sph[i].r[1], sph[i].r[2],
        "Dmtr", "Tmtr", sph[i].type, i, sph[i].type) == EOF){
        fprintf(stderr, fmt_writting_failed, __func__, "sphere (povray)");
        fclose(file);
        return EXIT_FAILURE;
      }
      exp_count++;
    }
  }
  fprintf(file, "}\n\n#declare inclusionSpheres = ");

  // If all spheres are exported, there is no inclusion particles to export
  if(exp_count == ns){
    fprintf(file, "0;\n");
    fclose(file);
    return EXIT_SUCCESS;
  }

  // In not, proceed with writing inclusion spheres
  fprintf(file, "union{\n");
  for(i=0; i<ns; i++){
    // Write sphere positions and properties
    if(sph[i].type > TYPE_INCLUSION_BASE){
      if(fprintf(file, fmt_exp_sph_0, sph[i].r[0], sph[i].r[1], sph[i].r[2],
        "Dinc", "Tinc", sph[i].type, i, sph[i].type) == EOF){
        fprintf(stderr, fmt_writting_failed, __func__, "sphere (povray)");
        fclose(file);
        return EXIT_FAILURE;
      }
      exp_count++;
    }
  }
  fprintf(file, "}\n");
  fclose(file);

  if(exp_count != ns){
    fprintf(stderr, fmt_data_corruption, __func__, exp_count, ns);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/*
 * povray_export_dimers(file, sph, nd)
 * Export dimer data in format readable by POV-Ray
 */
static int povray_export_dimers(FILE *file, DIM3D *dim, int nd){
  extern const int TYPE_MATRIX_LIMIT;
  extern const int TYPE_INCLUSION_BASE;
  extern const char *fmt_writting_failed;
  const char *fmt_data_corruption =
    " [%s] ERR: incomplete dimer data (%d/%d) exported to file\n";
  const char *fmt_exp_dim_0 =
    " dimSkeleton(<% .6le, % .6le, % .6le>, <% .6le, % .6le, % .6le>, <% .6le, % .6le, % .6le>, SR, TL, %s%04d) // dim %5d type %d\n";
  int i,j, exp_count = 0;
  double c0[3], c1[3];

  fprintf(file, "#declare matrixDimers = union{\n");
  for(i=0; i<nd; i++){
    // Write dimer positions and properties
    if(dim[i].type <= TYPE_MATRIX_LIMIT){
      for(j=0; j<3; j++){
        c0[j] = dim[i].R[j] + dim[i].O[j] * dim[i].L / 2e0;
        c1[j] = dim[i].R[j] - dim[i].O[j] * dim[i].L / 2e0;
      }
      if(fprintf(file, fmt_exp_dim_0,
          c0[0], c0[1], c0[2],
          c1[0], c1[1], c1[2],
          dim[i].R[0], dim[i].R[1], dim[i].R[2],
         "TmtrDm", dim[i].type, i, dim[i].type) == EOF){
        fprintf(stderr, fmt_writting_failed, __func__, "dimer (povray)");
        fclose(file);
        return EXIT_FAILURE;
      }
      exp_count++;
    }
  }
  fprintf(file, "}\n\n#declare inclusionDimers = ");

  // If all dimers are exported, there is no inclusion particles to export
  if(exp_count == nd){
    fprintf(file, "0;\n");
    fclose(file);
    return EXIT_SUCCESS;
  }

  // In not, proceed with writing inclusion dimers
  fprintf(file, "union{\n");
  for(i=0; i<nd; i++){
    // Write dimer positions and properties
    if(dim[i].type > TYPE_INCLUSION_BASE){
      for(j=0; j<3; j++){
        c0[j] = dim[i].R[j] + dim[i].O[j] * dim[i].L / 2e0;
        c1[j] = dim[i].R[j] - dim[i].O[j] * dim[i].L / 2e0;
      }
      if(fprintf(file, fmt_exp_dim_0,
          c0[0], c0[1], c0[2],
          c1[0], c1[1], c1[2],
          dim[i].R[0], dim[i].R[1], dim[i].R[2],
         "TincDm", dim[i].type, i, dim[i].type) == EOF){
        fprintf(stderr, fmt_writting_failed, __func__, "dimer (povray)");
        fclose(file);
        return EXIT_FAILURE;
      }
      exp_count++;
    }
  }
  fprintf(file, "}\n");
  fclose(file);

  if(exp_count != nd){
    fprintf(stderr, fmt_data_corruption, __func__, exp_count, nd);
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/*
 * export_dimers()
 * Writes information about dimers to a file. Perform the write operation twice
 * (add more if required) in order to separate dimers of different type.
 */
static int export_dimers(FILE *file, DIM3D *dim, int nd)
{
  extern const int TYPE_DIMER;
  extern const int TYPE_INCLUSION_DIMER;
  extern const char *fmt_writting_failed;
  const char *fmt_exp_dim =
    "%5d  %3d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %5d %5d\n";
  int i, k=0;

  for(i=0; i<nd; i++){
    // Export information about molecule
    if(dim[i].type == TYPE_DIMER || dim[i].type == TYPE_INCLUSION_DIMER){
      if(fprintf(file, fmt_exp_dim,
          k++, dim[i].type,
          dim[i].R[0], dim[i].R[1], dim[i].R[2],
          dim[i].O[0], dim[i].O[1], dim[i].O[2],
          dim[i].L,
          dim[i].sph_ind[0], dim[i].sph_ind[1]) == EOF){
        fprintf(stderr, fmt_writting_failed, __func__, "dimer");
        fclose(file);
        return EXIT_FAILURE;
      }
    }
  }
  fclose(file);
  return EXIT_SUCCESS;
}

/*
 * export_spheres(f,sp,ns)
 * Export positions and other data for spheres
 */
static int export_spheres(FILE *file, SPH *sph, int ns)
{
  extern const char *fmt_writting_failed;
  const char *fmt_exp_sph_0 = "%5d %3d % .16le % .16le % .16le %.16le ";
  const char *fmt_exp_sph_1 =
    "%5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d %5d  ";
  const char *fmt_exp_sph_2 = "  %2d %2d %2d\n";
  int i;
  for(i=0; i<ns; i++){
    // Write sphere positions and properties
    if(fprintf(file, fmt_exp_sph_0,
        i, sph[i].type,
        sph[i].r[0], sph[i].r[1], sph[i].r[2], sph[i].d) == EOF){
      fprintf(stderr, fmt_writting_failed, __func__, "sphere positions");
      fclose(file);
      return EXIT_FAILURE;
    }

    // Write sphere neighbors table
    if(fprintf(file, fmt_exp_sph_1,
        sph[i].ngb[0], sph[i].ngb[1], sph[i].ngb[2], sph[i].ngb[3],
        sph[i].ngb[4], sph[i].ngb[5], sph[i].ngb[6], sph[i].ngb[7],
        sph[i].ngb[8], sph[i].ngb[9], sph[i].ngb[10], sph[i].ngb[11]) == EOF){
      fprintf(stderr, fmt_writting_failed, __func__, "sphere neighbors");
      fclose(file);
      return EXIT_FAILURE;
    }

    // Write sphere lattice indexes
    if(fprintf(file, fmt_exp_sph_2,
      sph[i].lattice_ind[0], sph[i].lattice_ind[1],
      sph[i].lattice_ind[2]) == EOF){
        fprintf(stderr, fmt_writting_failed, __func__,
                "sphere lattice idnices");
        fclose(file);
        return EXIT_FAILURE;
    }
  }
  fclose(file);
  return EXIT_SUCCESS;
}

/*
 * export_to_GLviewer()
 * Export required data for viewing the structure by data_visGL
 * double box[3],  int ns,  int nd
 */
static int export_to_GLviewer(MODEL md, DIM3D *dim, SPH *sph, int strn){
#ifdef DEBUG_MODE
  extern const char *fmt_dbg_opening_file;
#endif
  extern const char *fmt_writting_failed;
  extern const char *fmt_write_notify;
  FILE *file = NULL;
  int i;
  char f_GLout[64];
  const char *exp_form = "%5d % .16le % .16le % .16le %.16le %d\n";
  const char *exp_form_d =
    "%5d % .16le % .16le % .16le % .16le % .16le % .16le %.16le %d\n";
  const char *data_visGL = "vcontrol.ini";
  const char *data_spheresGL = "gl_s3d%05d.csv";
  const char *data_dimersGL = "gl_d3d%05d.csv";

  // Prepare control file for data_visGL program
#ifdef DEBUG_MODE
  fprintf(stdout, fmt_dbg_opening_file, __func__, data_visGL);
#endif
  file = open_to_write(data_visGL);

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
#ifdef DEBUG_MODE
  fprintf(stdout, fmt_dbg_opening_file, __func__, f_GLout);
#endif
  file = open_to_write(f_GLout);

  // Export structure data to data_spheresGL
  fprintf(stdout, fmt_write_notify, "sphere", f_GLout);
  for(i=0; i<md.Nsph; i++){
    if(fprintf(file, exp_form, i,
        sph[i].r[0], sph[i].r[1], sph[i].r[2],
        sph[i].d,
        legacy_GLexport_sphere_type_converter(sph[i].type) ) == EOF){
      fprintf(stderr, fmt_writting_failed, __func__, "sphere");
      return EXIT_FAILURE;
    }
  }
  fclose(file);

  // Open file for dimer data
  sprintf(f_GLout, data_dimersGL, strn);
#ifdef DEBUG_MODE
  fprintf(stdout, fmt_dbg_opening_file, __func__, f_GLout);
#endif
  file = open_to_write(f_GLout);

  // Export structure data to data_dimersGL
  fprintf(stdout, fmt_write_notify, "dimer",f_GLout);
  for(i=0; i<md.Ndim; i++){
    if(fprintf(file, exp_form_d, i,
        dim[i].R[0], dim[i].R[1], dim[i].R[2],
        dim[i].O[0], dim[i].O[1], dim[i].O[2],
        dim[i].L,
        legacy_GLexport_dimer_type_converter(dim[i].type) ) == EOF){
      fprintf(stderr, fmt_writting_failed, __func__, "dimer");
      return EXIT_FAILURE;
    }
  }
  fclose(file);

  return EXIT_SUCCESS;
}

/* # SUB-SEC ########## FILE OUTPUT - DATA CONVERTERS ####################### */

/*
 * Old (but still used) export functions
 */
static int legacy_GLexport_dimer_type_converter(int type){
  extern const int TYPE_DIMER;
  extern const int TYPE_INCLUSION_DIMER;
  if(type == TYPE_DIMER) return 1;
  if(type == TYPE_INCLUSION_DIMER) return 1;
  return type;
}

static int legacy_GLexport_sphere_type_converter(int type){
  extern const int TYPE_SPHERE;
  extern const int TYPE_SPHERE_DIMER;
  extern const int TYPE_INCLUSION_SPHERE;
  extern const int TYPE_INCLUSION_SPHERE_DIMER;

  if(type == TYPE_SPHERE) return 3;
  if(type == TYPE_INCLUSION_SPHERE) return 2;
  if(type == TYPE_INCLUSION_SPHERE_DIMER) return 2;
  if(type == TYPE_SPHERE_DIMER) return 1;
  return -1;
}

/* # SEC ############## LOW LEVEL FILE OPERATIONS ########################### */

/*
 * open_to_read(file)
 * Open file to read and return the pointer to file. Terminate program on error.
 */
FILE* try_open_to_read(const char *file){
  return open_file(file, "r", 0);
}

FILE* open_to_read(const char *file){
  return open_file(file, "r", 1);
}

FILE* open_to_write(const char *file){
  return open_file(file, "w", 1);
}

static FILE* open_file(const char *file, const char *mode, int strict){
  const char *fmt_open_failed = " [%s] ERR: cannot open file: %s (%c)\n";
  FILE *p = NULL;

  if((p = fopen(file, mode)) == NULL && strict == 1){
    fprintf(stderr, fmt_open_failed, __func__, file, mode);
    exit(EXIT_FAILURE);
  }

  return p;
}


