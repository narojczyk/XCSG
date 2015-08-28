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
#include "globals.h"
#include "data.h"

#include "initials.h"

char *prog_name;

int main(int argc, char *argv[])
{
#ifdef PUSH_TRAJECTORIES
  FILE **t;
  const char* data_visGL="vcontrol.ini";
#endif

  FILE *f;
  
  SPH *spheres;

  char *f_ini = NULL;
  char t_out[15];

  int exit_status;
  int Ns, Nd;
  
  int i, m, n;
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
  
  // Set the number of spheres and dimres
  Ns = 4 * i_edge_fcc_N * i_edge_fcc_N * i_edge_fcc_N;
  Nd = Ns / 2;
  
  // Allocate memory for wolkers
  spheres = malloc( Ns * sizeof(SPH));
  
  // < DEBUG
  printf("Config variables:\n");
  printf("edge cells: %d\n",i_edge_fcc_N);
  printf("chanel %d %d %d\n", i_chanel[0], i_chanel[1], i_chanel[2]);
  printf("chanel radius %.16le\n",i_chanel_R);
  
  printf("Ns = %d, Nd = %d\n",Ns, Nd);
  // DEBUG >
  
  // Set fcc structure of spheres
  sph_set_fcc( spheres, Ns, i_edge_fcc_N);
  
  for(i=0;i<Ns;i++){
  printf("%lf %lf %lf\t %4d\n",spheres[i].r[0],spheres[i].r[1],spheres[i].r[2],i);
  }


cleanup:
  // Free resources
#ifdef PUSH_TRAJECTORIES
  for(i=0; i<initial_N; i++){
    fclose(t[i]);
  }
  free(t);
#endif

#ifdef PUSH_MSD
  fclose(g);
#endif

  free(spheres);
  return exit_status;
}


/* vim: set tw=80 ts=2 sw=2 et: */
