/*
 *
 * TODO:
 * 3. fix elastic colisions for friendly colisions
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

#ifdef PUSH_MSD
  const char* msdcsv="msd.csv";
  FILE *g;
#endif

  FILE *f;

  Walker *walkers, *johnn, *martin, *searcher, *target;

  char *f_ini = NULL;
  char t_out[15];

  int exit_status;
  int N;
  int i, m, n;
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

//  free(walkers);
  return exit_status;
}


/* vim: set tw=80 ts=2 sw=2 et: */
