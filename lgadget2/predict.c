#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef HPM
#include <libhpm.h>
#endif

/*! \file predict.c
 *  \brief drifts particles
 */

/* This function drifts all particles from the current time
 * to the future:  time0 - > time1
 *
 * Also, the tree nodes are drifted and updated accordingly.
 *
 * NOTE: For periodic boundary conditions we do the mapping of
 * coordinates onto the interval [0, All.BoxSize] only before the
 * domain decomposition, or for output to snapshot files.  This
 * simplifies dynamic tree updates, and allows the domain
 * decomposition to be carried out only every once in a while.  In the
 * force computation and neighbour search the correct nearest image
 * will nevertheless always be selected.
 */
void move_particles(int time0, int time1)
{
  int i, j;
  double dt_drift;
  double t0, t1;
  double iso = All.TimeBegin * exp(0.5 * (time0 + time1) * All.Timebase_interval);
  double anifac[3];
  eval_aniss(iso, anifac);
  for(i=0; i<3; ++i)
      anifac[i] = iso*iso / (anifac[i]*anifac[i]);

#ifdef LIGHTCONE
  double ThisBoxTime;
  double ThisBox_r;
  double ThisBox_r_prev;
  int imageCoverage;
  float LCPrevPos[3];
#endif

  t0 = second();

#ifdef HPM
  hpmStart(9, "Predict");
#endif

  dt_drift = get_drift_factor(time0, time1);

#ifdef LIGHTCONE
  ThisBoxTime = All.TimeBegin * exp(time0*All.Timebase_interval);
  ThisBox_r_prev = RofZ(1.0/ThisBoxTime - 1.0);

  ThisBoxTime = All.TimeBegin * exp(time1*All.Timebase_interval);
  ThisBox_r = RofZ(1.0/ThisBoxTime - 1.0);

  if(ThisTask == 0)
    printf("checking lightcone for particle with time %g, %g %g...\n", ThisBoxTime, ThisBox_r, All.BoxSize);

  if(ThisBox_r <= All.BoxSize*sqrt(3.0) && time0 != time1)
    {
      //init light cone interp functions
      init_lc_interp(time0,time1);

      if(ThisBox_r_prev < All.BoxSize)
	imageCoverage = 0;
      else
	imageCoverage = 1;

      for(i = 0; i < NumPart; i++)
	{
	  LCPrevPos[0] = P[i].Pos[0];
	  LCPrevPos[1] = P[i].Pos[1];
	  LCPrevPos[2] = P[i].Pos[2];

	  for(j = 0; j < 3; j++)
	    P[i].Pos[j] += P[i].Vel[j] * dt_drift * anifac[j];  /* correct for lightcones? */

	  check_particle(i, LCPrevPos, imageCoverage);
	}
    }
  else
    {
#endif
      for(i = 0; i < NumPart; i++)
	{
	  for(j = 0; j < 3; j++)
	    P[i].Pos[j] += P[i].Vel[j] * dt_drift * anifac[j];
	}
#ifdef LIGHTCONE
    }
#endif

#ifdef HPM
  hpmStop(9);
#endif

  t1 = second();

  All.CPU_Predict += timediff(t0, t1);

#ifdef LIGHTCONE
#ifndef LIGHTCONE_INTERNAL_TIMER
  LCProf.CPU_CheckLC = timediff(t0, t1);
#endif
#endif
}



/*  This function makes sure that all particles coordinates (Pos) are
 *  mapped onto the interval [0, BoxSize].
 *  After this function has been called, a new domain decomposition
 *  should be done, or at least a new force-tree needs to be constructed.
 */
void do_box_wrapping(void)
{
  int i, j;
  double boxsize;

  boxsize = All.BoxSize;

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize;

	while(P[i].Pos[j] >= boxsize)
	  P[i].Pos[j] -= boxsize;
      }
}
