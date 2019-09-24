#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

/*! \file global.c
 *  \brief computes some global quantities
 */

/* This routine computes various global properties of the particle distribution and stores the result in the
 * struct `SysState'.  Currently, not all the information that's computed here is actually used (e.g. momentum
 * is not really used anywhere), just the energies are written to a log-file every once in a while.
 */
double compute_mean_rms_velocity(void)
{
  int i, j;
  double rms_sum, rms_tot;

  for(i = 0, rms_sum = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  rms_sum += P[i].Vel[j] * P[i].Vel[j];
	}
    }

  /* some the stuff over all processors */
  MPI_Allreduce(&rms_sum, &rms_tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return sqrt(rms_tot / All.TotNumPart);
}
