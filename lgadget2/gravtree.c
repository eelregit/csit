#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#ifdef HPM
#include <libhpm.h>
#endif

/*! \file gravtree.c
 *  \brief gets short-range tree force for all active particles
 */


static double factor, hubble_a;

static double LocalTreeAllocFactor = 0;

/*  This function computes the gravitational forces for all active particles.  A new tree is constructed, if
 *  the number of force computations since it's last construction exceeds some fraction of the total particle
 *  number, otherwise tree nodes are dynamically updated if needed.
 */
void gravity_tree(void)
{
  long long ntot, ntotleft;
  int numForceUpdate;
  int numnodes, nexportsum = 0;
  int i, j, k, iter = 0, nac, acindex;
  int *numlist, *numnodeslist, *nrecv, maxnumnodes, nexport;
  int *noffset, *nbuffer, *nsend, *nsend_local, *ndonelist;
  int ngrp, place;
  int ndone, maxfill, maxnodes;
  int level, sendTask, recvTask;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double costtotal, *costtreelist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm, megs, megsum;
  size_t allbytes;
  FLOAT acc[3];
  MPI_Status status;


  /* set new softening lengths */
  set_softenings();

#ifdef ONLY_PM
  pm_only();
  return;
#endif

  factor = 1 / (All.Time * All.Time);
  hubble_a = hubble_function(All.Time);


  /* Determine 'numForceUpdate', i.e. the number of particles on this processor that want a force update */
  for(i = 0, numForceUpdate = 0; i < NumPart; i++)
    {
      if(P[i].Ti_endstep == All.Ti_Current)
	numForceUpdate++;
    }

  /* Note: We will assume that numForceUpdate still fits into a 32 bit
     integer, but the sum over all processors can be larger. We therefore need
     to use 64-bit for "ntot", meaning that we manually construct the sum and
     cannot simply use MPI_Allreduce */

  numlist = mymalloc(sizeof(int) * NTask);
  MPI_Allgather(&numForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  myfree(numlist);


  if(ThisTask == 0)
    {
      printf("Tree construction...\n");
      fflush(stdout);
    }

  tstart = second();

  if(LocalTreeAllocFactor == 0)
    LocalTreeAllocFactor = All.TreeAllocFactor;

  do
    {
      maxnodes = LocalTreeAllocFactor * NumPart;

      allbytes = force_treeallocate(maxnodes, NumPart);

      numnodes = force_treebuild();

      if(numnodes >= maxnodes)
	{
	  LocalTreeAllocFactor *= 1.05;
	  printf("on task=%d: need to increase treeallocfactor to %g\n", ThisTask, LocalTreeAllocFactor);
	  fflush(stdout);

	  if(LocalTreeAllocFactor > 2.0)
	    {
	      printf("Task=%d: Failure to construct tree. Something's wrong here.\n", ThisTask);
	      endrun(12);
	    }

	  force_treefree();
	}

      if(numnodes < 0.9 * maxnodes)
	LocalTreeAllocFactor /= 1.05;
    }
  while(numnodes >= maxnodes);

  megs = allbytes / (1024.0 * 1024.0);
  MPI_Reduce(&megs, &megsum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  if(ThisTask == 0)
    {
      printf("Use %g MByte for BH-tree.\n", megsum / NTask);
      printf("Tree construction done. (took %g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }


  if(ThisTask == 0)
    {
      printf("Begin tree force.\n");
      fflush(stdout);
    }


  allocate_commbuffers();
  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);
  ndonelist = mymalloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */
  costtotal = 0;


  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
#ifdef HPM
      hpmStart(5, "TreeWalk");
#endif
      for(nexport = 0, ndone = 0, nac = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    costtotal += force_treeevaluate_shortrange(i, 0, acc);

	    for(j = 0, acindex = -1; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    if(acindex < 0)
		      {
			acindex = nac;
			GravDataAccTable[acindex].Pindex = i;
			GravDataAccTable[acindex].Acc[0] = acc[0];
			GravDataAccTable[acindex].Acc[1] = acc[1];
			GravDataAccTable[acindex].Acc[2] = acc[2];
			nac++;
		      }

		    for(k = 0; k < 3; k++)
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
		    GravDataGet[nexport].OldAcc = P[i].OldAcc;

		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = acindex;
		    GravDataIndexTable[nexport].SortIndex = nexport;
		    nexport++;
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }

	    if(acindex < 0)
	      kick_particle(i, acc);
	  }
#ifdef HPM
      hpmStop(5);
#endif
      tend = second();

      timetree += timediff(tstart, tend);

      if(ThisTask == 0)
	{
	  printf("on task 0, did %d local particles in %g sec.\n", ndone, timediff(tstart, tend));
	  fflush(stdout);
	}

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);


      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);


	  tstart = second();
#ifdef HPM
	  hpmStart(5, "TreeWalk");
#endif
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
	      costtotal += force_treeevaluate_shortrange(j, 1, acc);
	    }
#ifdef HPM
	  hpmStop(5);
#endif

	  tend = second();
	  timetree += timediff(tstart, tend);

	  if(ThisTask == 0)
	    {
	      printf("on task 0, did %d imported particles in %g sec.\n", nbuffer[ThisTask],
		     timediff(tstart, tend));
	      fflush(stdout);
	    }

	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      /* process the results */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    GravDataAccTable[place].Acc[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      /* now kick the exported particles */
      for(j = 0; j < nac; j++)
	{
	  kick_particle(GravDataAccTable[j].Pindex, GravDataAccTable[j].Acc);
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  myfree(ndonelist);
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);
  free_commbuffers();

  if(ThisTask == 0)
    {
      printf("tree force is done.\n");
      fflush(stdout);
    }

  tstart = second();
  force_costevaluate();		/* assign accumulated cost to particles */
  tend = second();

  if(ThisTask == 0)
    {
      printf("done cost evaluation (%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }

#ifdef TWOPOINT
  if(All.TwoPointFlag)
    {
      twopoint();
      All.TwoPointFlag = 0;
    }
#endif

  force_treefree();		/* free memory for tree */


  /* This will switch to the relative opening criterion */
  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;


  /*  gather some diagnostic performance information */

  timetreelist = mymalloc(sizeof(double) * NTask);
  timecommlist = mymalloc(sizeof(double) * NTask);
  costtreelist = mymalloc(sizeof(double) * NTask);
  numnodeslist = mymalloc(sizeof(int) * NTask);
  nrecv = mymalloc(sizeof(int) * NTask);

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled(PE0): %g\n", maxnumnodes, numnodes / ((double) maxnodes));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }

  myfree(nrecv);
  myfree(numnodeslist);
  myfree(costtreelist);
  myfree(timecommlist);
  myfree(timetreelist);
}




void kick_particle(int i, FLOAT * acc)
{
  double ax, ay, az, aax, aay, aaz, ac, dt, dt_gravkick;
  int ti_step, ti_min, tstart, tend;

  ax = acc[0];
  ay = acc[1];
  az = acc[2];

  if(All.PM_Ti_begstep == All.Ti_Current)
    {
      P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az + P[i].GravCost);
      P[i].GravCost = 0.01;
    }

  ax *= All.G * All.PartMass;
  ay *= All.G * All.PartMass;
  az *= All.G * All.PartMass;

  aax = factor * ax;
  aay = factor * ay;
  aaz = factor * az;
  ac = sqrt(aax * aax + aay * aay + aaz * aaz);

  dt = sqrt(2 * All.ErrTolIntAccuracy * All.Time * All.ComovSoftening / ac);

  /* convert the physical timestep to dloga if needed */

  dt *= hubble_a;

  if(dt >= All.MaxSizeTimestep)
    dt = All.MaxSizeTimestep;

  if(dt >= All.Dt_displacement)
    dt = All.Dt_displacement;

  if(dt < All.MinSizeTimestep)
    {
#ifndef NOSTOP_WHEN_BELOW_MINTIMESTEP
      printf("warning: Timestep wants to be below the limit `MinSizeTimestep'\n");
      printf("Part-ID= %d%09d dt=%g ac=%g xyz=(%g|%g|%g)\n",
	     (int) (P[i].ID / 1000000000),
	     (int) (P[i].ID % 1000000000), dt, ac, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
      fflush(stdout);
      endrun(888);
#endif
      dt = All.MinSizeTimestep;
    }

  ti_step = dt / All.Timebase_interval;

  /* make it a power 2 subdivision */
  ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  ti_step = ti_min;

  /* do synchronization */
  if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
    {
      if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
	ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
    }

  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
    ti_step = 0;

  tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;
  tend = P[i].Ti_endstep + ti_step / 2;

  P[i].Ti_begstep = P[i].Ti_endstep;
  P[i].Ti_endstep = P[i].Ti_begstep + ti_step;

  /* do the kick */
#ifdef MAKEGLASS
  dt_gravkick = 1;
#else
  dt_gravkick = get_gravkick_factor(tstart, tend);
#endif

  P[i].Vel[0] += ax * dt_gravkick;
  P[i].Vel[1] += ay * dt_gravkick;
  P[i].Vel[2] += az * dt_gravkick;
}



void pm_only(void)
{
  double dt;
  int i, ti_step, ti_min, tstart, tend;

  dt = All.MaxSizeTimestep;

  if(dt >= All.Dt_displacement)
    dt = All.Dt_displacement;

  ti_step = dt / All.Timebase_interval;

  /* make it a power 2 subdivision */
  ti_min = TIMEBASE;
  while(ti_min > ti_step)
    ti_min >>= 1;
  ti_step = ti_min;
  
  //FIXME - variable i is used without a value set - setting to zero as a fix
  i = 0;
  
  /* do synchronization */
  if(ti_step > (P[i].Ti_endstep - P[i].Ti_begstep))	/* timestep wants to increase */
    {
      if(((TIMEBASE - P[i].Ti_endstep) % ti_step) > 0)
	ti_step = P[i].Ti_endstep - P[i].Ti_begstep;	/* leave at old step */
    }

  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
    ti_step = 0;

  for(i=0; i < NumPart; i++)
    {
      tstart = (P[i].Ti_begstep + P[i].Ti_endstep) / 2;
      tend = P[i].Ti_endstep + ti_step / 2;

      P[i].Ti_begstep = P[i].Ti_endstep;
      P[i].Ti_endstep = P[i].Ti_begstep + ti_step;
    }
}



/*  This function sets the (comoving) softening length of all particle species in the table
 *  All.SofteningTable[...]  We check here that the proper softening length is bounded by the ..MaxPhys
 *  values.
 */
void set_softenings(void)
{
  if(All.Softening * All.Time > All.SofteningMaxPhys)
    All.ComovSoftening = All.SofteningMaxPhys / All.Time;
  else
    All.ComovSoftening = All.Softening;
}


int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}
