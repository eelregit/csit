#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>

#include "allvars.h"
#include "proto.h"


/*! \file run.c
 *  \brief  iterates over timesteps, main loop
 */


/* This routine contains the main simulation loop that iterates over
 * the single timesteps. The loop terminates when the cpu-time
 * limit is reached, when a `stop' file is found in the output 
 * directory, or when the simulation ends because we arrived
 * at TimeMax
 */
void run(void)
{
  FILE *fd;
  int stopflag = 0;
  char buf[200], stopfname[200];
  double t0, t1;

  sprintf(stopfname, "%sstop", All.OutputDir);

  do				/* main loop */
    {
      t0 = second();

      find_next_sync_point_and_drift();	/* find next synchronization point and drift particles to this time.
					 * If needed, this function will also write an output file
					 * at the desired time.
					 */
      every_timestep_stuff();	/* write some info to log-files */


      DomainDecomposition(0);	/* do domain decomposition */


      compute_accelerations(0);	/* compute accelerations for 
				 * the particles that are to be advanced,
				 * kick them on the fly, and assign new timesteps
				 */

      if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics)	/* check whether we want a full statistics */
	{
	  energy_statistics();	/* compute and output energy statistics */
	  All.TimeLastStatistics += All.TimeBetStatistics;
	}

      /* advance particles and compute new timesteps for them */

      All.NumCurrentTiStep++;

      if(ThisTask == 0)		/* Check whether we need to interrupt the run */
	{
	  if((fd = fopen(stopfname, "r")))	/* Is the stop-file present? If yes, interrupt the run. */
	    {
	      fclose(fd);
	      stopflag = 1;
	      unlink(stopfname);
	    }

	  if(CPUThisRun > 0.85 * All.TimeLimitCPU)	/* are we running out of CPU-time ? If yes, interrupt run. */
	    {
	      printf("reaching time-limit. stopping.\n");
	      printf("restart=1\n");
	      stopflag = 2;
	    }
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag)
	{
	  restart(0);		/* write restart file */
	  MPI_Barrier(MPI_COMM_WORLD);

	  if(stopflag == 2 && All.ResubmitOn && ThisTask == 0)
	    {
	      close_outputfiles();
	      sprintf(buf, "%s", All.ResubmitCommand);
	      system(buf);
	    }
	  return;
	}

      if(ThisTask == 0)		/* is it time to write a restart-file? (for security) */
	{
	  if((CPUThisRun - All.TimeLastRestartFile) >= All.CpuTimeBetRestartFile)
	    {
	      All.TimeLastRestartFile = CPUThisRun;
	      stopflag = 3;
	    }
	  else
	    stopflag = 0;
	}

      MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);

      if(stopflag == 3)
	{
	  restart(0);		/* write an occasional restart file */
	  stopflag = 0;
	}

      t1 = second();

      All.CPU_Total += timediff(t0, t1);
      CPUThisRun += timediff(t0, t1);
    }
  while(All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax);

  restart(0);

  savepositions(All.SnapshotFileCount++);	/* write a last snapshot file at
						 * final time (will be overwritten if
						 * All.TimeMax is increased and the
						 * run is continued) 
						 */

#ifdef LIGHTCONE
  Finalize_LightCone();
#endif

#if defined(TWOPOINT) || defined (POWERSPEC)
  compute_accelerations(0);                     /* compute accelerations once more since this
						 * will trigger power spectrum computation and
						 * two-point correlation function at the final time 
						 */
#endif

}


void find_next_sync_point_and_drift(void)
{
  int n, min, min_glob, flag;
  double timeold;
  double t0, t1;

  t0 = second();

  timeold = All.Time;

  for(n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
    if(min > P[n].Ti_endstep)
      min = P[n].Ti_endstep;

  MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  /* We check whether this is a full step where all particles are synchronized */
  flag = 1;
  for(n = 0; n < NumPart; n++)
    if(P[n].Ti_endstep > min_glob)
      flag = 0;

  MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  if(min_glob > All.PM_Ti_endstep)
    {
      min_glob = All.PM_Ti_endstep;
      Flag_FullStep = 1;
    }

  t1 = second();

  All.CPU_Predict += timediff(t0, t1);

  while(min_glob >= All.Ti_nextoutput && All.Ti_nextoutput >= 0)
    {
      move_particles(All.Ti_Current, All.Ti_nextoutput);

      All.Ti_Current = All.Ti_nextoutput;

      All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
      eval_aniss(All.Time, All.TimeAni);

      savepositions(All.SnapshotFileCount++);	/* write snapshot file */

      All.Ti_nextoutput = find_next_outputtime(All.Ti_nextoutput + 1);
    }

  move_particles(All.Ti_Current, min_glob);

  All.Ti_Current = min_glob;

  All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
  eval_aniss(All.Time, All.TimeAni);

  All.TimeStep = All.Time - timeold;
  
#ifdef LIGHTCONE
  //this code helps to profile the light cone generator
  //actual profile info is output to lightcone.txt in every_timestep_stuff below
  //could do it here to be cleaner, but that would disrupt the structure of the code too much
  //only doing reduce calls here so that light cone profile info falls in correct order when looking at STDOUT
  long long NumLCTests,NumLCFound;
  double CPU_CheckLC;
  double CPU_IO;

  //in some MPI libraries these reduce calls may have the side effect of producing effectively a barrier,
  //  though one should not count on this!!
  //if you remove this reduce statement and comment out code marked below - then it will just profile task 0
  MPI_Reduce(&(LCProf.NumLCTests),&NumLCTests,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&(LCProf.NumLCFound),&NumLCFound,1,MPI_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&(LCProf.CPU_CheckLC),&CPU_CheckLC,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&(LCProf.CPU_IO),&CPU_IO,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  
  if(ThisTask == 0)
    {
      if(NumLCTests > 0)
        printf("all tasks tested parts %lld times for LC crossing in %g sec (%g tests/sec, found %lld parts in LC).\n",
	       NumLCTests,CPU_CheckLC/NTask,((double) (NumLCTests))/CPU_CheckLC/NTask,NumLCFound);
      
      //save info on task zero to be output below in every_timestep_stuff()
      //comment this out if you remove the MPI_Reduce calls below
      LCProf.NumLCTests = NumLCTests;
      LCProf.NumLCFound = NumLCFound;
      LCProf.CPU_CheckLC = CPU_CheckLC;
      LCProf.CPU_IO = CPU_IO;
    }
  else
    {
      //reset profile info for other tasks here
      //task 0 profile info is reset in every_timestep_stuff()
      LCProf.NumLCTests = 0;
      LCProf.NumLCFound = 0;
      LCProf.CPU_CheckLC = 0.0;
      LCProf.CPU_IO = 0.0;
    }
#endif
  
  if(Flag_FullStep)
    find_dt_displacement_constraint();
}




/*! This function computes an upper limit ('dt_displacement') to the global timestep of the system based on
 *  the rms velocities of particles. For cosmological simulations, the criterion used is that the rms
 *  displacement should be at most a fraction MaxRMSDisplacementFac of the mean particle separation. Note that
 *  the latter is estimated using the assigned particle masses, separately for each particle type. If comoving
 *  integration is not used, the function imposes no constraint on the timestep.
 */
void find_dt_displacement_constraint(void)
{
  int i;
  double v, v_sum;
  double dt, dmean, asmth = 0;
  double hfac, hubble_a;

  hubble_a = hubble_function(All.Time);

  hfac = All.Time * All.Time * hubble_a;

  All.Dt_displacement = All.MaxSizeTimestep;



  v = 0;

  for(i = 0; i < NumPart; i++)
    {
      v += P[i].Vel[0] * P[i].Vel[0] + P[i].Vel[1] * P[i].Vel[1] + P[i].Vel[2] * P[i].Vel[2];
    }

  MPI_Allreduce(&v, &v_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  dmean = pow(All.PartMass / (All.Omega0 * 3 * All.Hubble * All.Hubble / (8 * M_PI * All.G)), 1.0 / 3);

  dt = All.MaxRMSDisplacementFac * hfac * dmean / sqrt(v_sum / All.TotNumPart);

  asmth = All.Asmth[0];
  if(asmth < dmean)
    dt = All.MaxRMSDisplacementFac * hfac * asmth / sqrt(v_sum / All.TotNumPart);

  if(ThisTask == 0)
    printf("dmean=%g asmth=%g mass=%g a=%g  sqrt(<p^2>)=%g  dlogmax=%g\n",
	   dmean, asmth, All.PartMass, All.Time, sqrt(v_sum / All.TotNumPart), dt);

  if(dt < All.Dt_displacement)
    All.Dt_displacement = dt;

  if(ThisTask == 0)
    printf("displacement time constraint: %g  (%g)\n", All.Dt_displacement, All.MaxSizeTimestep);
}





/* this function returns the next output time
 * that is equal or larger to ti_curr
 */
int find_next_outputtime(int ti_curr)
{
  int i, ti, ti_next;
  double next, time;

  ti_next = -1;

  if(All.OutputListOn)
    {
      for(i = 0; i < All.OutputListLength; i++)
	{
	  time = All.OutputListTimes[i];

	  ti = log(time / All.TimeBegin) / All.Timebase_interval;

	  if(ti >= ti_curr)
	    {
	      if(ti_next == -1)
		ti_next = ti;

	      if(ti_next > ti)
		ti_next = ti;
	    }
	}
    }
  else
    {
      time = All.TimeOfFirstSnapshot;

      do
	{
	  ti = log(time / All.TimeBegin) / All.Timebase_interval;

	  time *= All.TimeBetSnapshot;
	}
      while(ti < ti_curr);

      ti_next = ti;
    }

  if(ti_next == -1)
    ti_next = 2 * TIMEBASE;	/* this will prevent any further output */

  next = All.TimeBegin * exp(ti_next * All.Timebase_interval);

  if(ThisTask == 0)
    printf("Setting next time for snapshot file to Time_next= %g\n\n", next);

  return ti_next;
}




/* This routine writes one line for every timestep to two log-files.
 * In FdInfo, we just list the timesteps that have been done,
 * while in FdCPU the cumulative cpu-time consumption in various parts
 * of the code is stored.
 */
void every_timestep_stuff(void)
{
  double z;

#ifdef DARKENERGY
  double hubble_a;
#endif
  
  if(ThisTask == 0)
    {
      z = 1.0 / (All.Time) - 1;
      
#ifdef LIGHTCONE
      //profile the Light cone generator per time step
      fprintf(FdLC,"%d %g %g %lld %10.2f %10.2f %lld\n",All.NumCurrentTiStep,All.Time,RofZ(z),
	      LCProf.NumLCTests,LCProf.CPU_CheckLC/NTask,LCProf.CPU_IO/NTask,LCProf.NumLCFound);
      fflush(FdLC);
      
      //reset the profile info for task 0 here, rest of tasks are done above...not clean, but it works
      LCProf.NumLCTests = 0;
      LCProf.NumLCFound = 0;
      LCProf.CPU_CheckLC = 0.0;
      LCProf.CPU_IO = 0.0;
#endif
      
      fprintf(FdInfo, "\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n",
	      All.NumCurrentTiStep, All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
      printf("\n\nBegin Step %d, Time: %g, Redshift: %g, Systemstep: %g, Dloga: %g\n", All.NumCurrentTiStep,
	     All.Time, z, All.TimeStep, log(All.Time) - log(All.Time - All.TimeStep));
      fflush(FdInfo);


      fprintf(FdCPU, "Step %d, Time: %g\n", All.NumCurrentTiStep, All.Time);

      fprintf(FdCPU,
	      "%10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n",
	      All.CPU_Total, All.CPU_Gravity, All.CPU_Domain,
	      All.CPU_Predict, All.CPU_FoF, All.CPU_Snapshot, All.CPU_TreeWalk, All.CPU_TreeConstruction,
	      All.CPU_CommSum, All.CPU_Imbalance, All.CPU_PM, All.CPU_Peano,
	      All.CPU_TwoPoint, All.CPU_PowerSpec);


      fflush(FdCPU);

#ifdef DARKENERGY
      hubble_a = hubble_function(All.Time);
      fprintf(FdDE, "%d %g %e ", All.NumCurrentTiStep, All.Time, hubble_a);
#ifndef TIMEDEPDE
      fprintf(FdDE, "%e ", All.DarkEnergyParam);
#else
      fprintf(FdDE, "%e %e ", get_wa(All.Time), DarkEnergy_a(All.Time));
#endif
      fprintf(FdDE, "\n");
      fflush(FdDE);
#endif
      
    }

  set_random_numbers();
}


/* This routine first calls a computation of various global quantities
 * of the particle distribution, and then writes some statistics
 * about the energies in the various particle components to the 
 * file FdEnergy.
 */
void energy_statistics(void)
{
  double vel_rms;

  vel_rms = compute_mean_rms_velocity();

  if(ThisTask == 0)
    {
      fprintf(FdEnergy, "%g %g \n", All.Time, vel_rms);
      fflush(FdEnergy);
    }
}
