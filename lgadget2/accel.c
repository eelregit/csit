#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef HPM
#include <libhpm.h>
#endif


/*! \file accel.c
 *  \brief driver routines to carry out force computation
 */


/* This routine computes the accelerations for all active particles.
 * First, the gravitational forces are computed (this also
 * reconstructs the tree, if needed). Also note that the gas-particle
 * tree will in any case be updated in its geometrical properties.
 *
 * If gas particles are presented, the `interior' of the local domain
 * is determined. This region is guaranteed to contain only particles
 * local to the processor. This information will be used to reduce
 * communication in the hydro part.  The density for active SPH
 * particles is computed next. If the number of neighbours should be
 * outside the allowed bounds, it will be readjusted by the function
 * ensure_neighbours(). Finally, the hydrodynamical forces are added.
 */
void compute_accelerations(int mode)
{
  double tstart, tend;

  if(ThisTask == 0)
    {
      printf("Start force computation...\n");
      fflush(stdout);
    }

  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      tstart = second();
#ifdef HPM
      hpmStart(2, "LongRange"); 
#endif
      pmforce_periodic();
#ifdef HPM
      hpmStop(2); 
#endif
      tend = second();
      All.CPU_PM += timediff(tstart, tend);
    }

  tstart = second();		/* measure the time for the full force computation */
#ifdef HPM
  hpmStart(3, "ShortRange"); 
#endif
  gravity_tree();		/* computes gravity accel. */
#ifdef HPM
  hpmStop(3);
#endif
  tend = second();
  All.CPU_Gravity += timediff(tstart, tend);

 
#ifdef MAKEGLASS
  glass_step();
#endif

  if(ThisTask == 0)
    {
      printf("force computation done.\n");
      fflush(stdout);
    }
}


#ifdef MAKEGLASS
void glass_step(void)
{
  int i, j;
  double disp, dispmax, globmax, dmean, fac, disp2sum, globdisp2sum;

  for(i = 0, dispmax = 0, disp2sum = 0; i < NumPart; i++)
    {
      disp = sqrt(P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2]);

      disp *= 2.0 / (3 * All.Hubble * All.Hubble);

      disp2sum += disp * disp;

      if(disp > dispmax)
	dispmax = disp;
    }

  MPI_Allreduce(&dispmax, &globmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&disp2sum, &globdisp2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(All.PartMass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  if(globmax > dmean)
    fac = dmean / globmax;
  else
    fac = 1.0;

  if(ThisTask == 0)
    {
      printf("\nglass-making:  dmean= %g  global disp-maximum= %g  rms= %g\n\n",
	     dmean, globmax, sqrt(globdisp2sum / All.TotNumPart));
      fflush(stdout);
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  P[i].Pos[j] -= fac * P[i].Vel[j] * 2.0 / (3 * All.Hubble * All.Hubble);
	  P[i].Vel[j] = 0;
	}
    }

}
#endif
