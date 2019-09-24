#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file init.c
 *  \brief initialization tasks for a simulation run
 */


/*
 *  init() reads in the initial conditions, and allocates storage for the tree(s).  An intial domain
 *  decomposition and force computation is done, followed by a second domain decomposition based on the
 *  initial work estimates. Then the first particle timesteps are determined. The simulation is set up for the
 *  timestep iteration in run().
 */
void init(void)
{
  int i, j;

  All.Time = All.TimeBegin;

  switch (All.ICFormat)
    {
    case 1:
#if (MAKEGLASS > 1) 
      seed_glass();
#else
      read_ic(All.InitCondFile);
#endif
      break;
    default:
      if(ThisTask == 0)
	printf("ICFormat=%d not supported.\n", All.ICFormat);
      endrun(0);
    }

  All.Time = All.TimeBegin;
  eval_aniss(All.Time, All.TimeAni);


  All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
  All.Ti_Current = 0;

  All.PM_Ti_endstep = All.PM_Ti_begstep = 0;

  All.NumCurrentTiStep = 0;	/* setup some counters */
  All.SnapshotFileCount = 0;
  if(RestartFlag == 2)
    All.SnapshotFileCount = atoi(All.InitCondFile + strlen(All.InitCondFile) - 3) + 1;

  All.TotNumOfForces = 0;

  All.TwoPointFlag = 0;
  All.PowerSpecFlag = 0;

  check_omega();

  All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;

  /*  change to new velocity variable */

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      P[i].Vel[j] *= sqrt(All.TimeAni[j]) * All.TimeAni[j];

  for(i = 0; i < NumPart; i++)	/*  start-up initialization */
    {
      P[i].Ti_endstep = 0;
      P[i].Ti_begstep = 0;

      P[i].OldAcc = 0;
      P[i].GravCost = 1;
    }

}


/* This routine computes the mass content of the box and
 * compares it to the specified value of Omega.
 * If discrepant, the run is terminated.
 */
void check_omega(void)
{
  double mass = 0, masstot, omega;
  int i;

  for(i = 0; i < NumPart; i++)
    mass = NumPart * All.PartMass;

  MPI_Allreduce(&mass, &masstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  omega =
    masstot / (All.BoxSize * All.BoxSize * All.BoxSize) / (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G));

  if(fabs(omega - All.Omega0) > 1.0e-3)
    {
      if(ThisTask == 0)
	{
	  printf("\n\nI've found something odd!\n");
	  printf
	    ("The mass content accounts only for Omega=%g,\nbut you specified Omega=%g in the parameterfile.\n",
	     omega, All.Omega0);
	  printf("\nI better stop.\n");

	  fflush(stdout);
	}
      endrun(1);
    }
}




#if (MAKEGLASS > 1)
void seed_glass(void)
{
  int i, k, n_for_this_task, iseed;
  double Range[3], LowerBound[3];
  double drandom;
  long long IDstart;

  All.TotNumPart = MAKEGLASS;
  All.PartMass = All.Omega0 * (3 * All.Hubble * All.Hubble / (8 * M_PI * All.G))
    * (All.BoxSize * All.BoxSize * All.BoxSize) / All.TotNumPart;

  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */
  allocate_memory();

  header.npartTotal[1] = All.TotNumPart;
  header.mass[1] = All.PartMass;

  if(ThisTask == 0)
    {
      printf("\nPartMass= %g\n", All.PartMass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
    }


  /* set the number of particles assigned locally to this task */
  n_for_this_task = All.TotNumPart / NTask;

  fprintf(stdout, "ThisTask %d n_for=%d\n", ThisTask, n_for_this_task);
  fflush(stdout);


  if(ThisTask == NTask - 1)
    n_for_this_task = All.TotNumPart - (NTask - 1) * n_for_this_task;

  NumPart = 0;
  IDstart = 1 + (All.TotNumPart / NTask) * ThisTask;

  /* split the temporal domain into Ntask slabs in z-direction */

  Range[0] = Range[1] = All.BoxSize;
  Range[2] = All.BoxSize / NTask;
  LowerBound[0] = LowerBound[1] = 0;
  LowerBound[2] = ThisTask * Range[2];


  fprintf(stdout, "ThisTask %d zmin=%g range=%g\n", ThisTask, LowerBound[2], Range[2]);
  fflush(stdout);

  srand48(ThisTask);

  for(i = 0; i < n_for_this_task; i++)
    {

      for(k = 0; k < 3; k++)
	{
	  iseed = IDstart + i;

	  drandom = drand48();

	  P[i].Pos[k] = LowerBound[k] + Range[k] * drandom;
	  P[i].Vel[k] = 0;
	}

      P[i].ID = IDstart + i;

      if(i < 10)
	fprintf(stdout, "ThisTask %d P[i].ID= %d%09d\n", ThisTask, 
		(int) (P[i].ID / 1000000000),
		(int) (P[i].ID % 1000000000));
      fflush(stdout);

      NumPart++;
    }
}
#endif

