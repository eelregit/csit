#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"

#ifdef HPM
#include <libhpm.h>
#endif


#define REDUC_FAC      0.98


/*! \file domain.c
 *  \brief establishes domain decomposition 
 */


static int *toGo;
static int *local_toGo;
static int *list_NumPart;
static int *list_load;
static double *list_work;
static int *startlist, *endlist;

/* toGo[task*NTask + partner] gives the number of particles in
   task 'task' that have to go to 'partner' */


static int maxload;




/* This function converts a float in the range [0,1[ to an integer in the
 * range [0,DOMAINGRID-1]. The reason for doing it in this elaborate way
 * is basically because on the Intel platform, a correct IEEE float to int
 * conversion is very expensive because it would require to (re)set the 
 * default rounding mode, which can only be done by flushing the FPU state.
 * As a result, most optimizing compilers do a somewhat fuzzy float-to-int
 * conversion, which can sometimes lead to a round-up, at other times to 
 * a round-down... We here need however reproducible rounding, hence we 
 * do things ourselves by exploiting float-number storage directly.
 *
 * The representation of a IEEE double-number is:
 * sign[1 bit]  exp[11 bits]  mantissa[52 bits]
 *
 * The value of the double number is computed as:
 *     (-1)^sign * 1.mantissa * 2^(exp - 1023)
 *
 * The representation of a IEEE float-number is:
 * sign[1 bit]  exp[8 bits]  mantissa[23 bits]
 *
 * The value of the double number is computed as:
 *     (-1)^sign * 1.mantissa * 2^(exp - 127)
 *
 * Note: For values in the range [1,2[, the exponent is constant!
 */
#define DOUBLE_to_DOMAINGRID(y) ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - DOMAINLEVELS)))

/* This is the main routine for the domain decomposition. The code will 
 * decompose each particle type separately. It will try to balance the 
 * work-load as defined by the sum of the P[i]-GravCost in each domain.
 * The decomposition will respect the maximum memory-imbalance given
 * by the value of PartAllocFactor.
 *
 */
void DomainDecomposition(int geom_flag)
{
  int i;
  double t0, t1;

  t0 = second();

#ifdef HPM
  hpmStart(8, "Domain");
#endif

#ifndef PHDOMAINDECOMP  
  geom_flag = 1;
#endif
  
  if(ThisTask == 0)
    {
      if(geom_flag != 0)
	printf("begin equal-volume domain decomposition... \n");
      else
	printf("begin domain decomposition... \n");
      fflush(stdout);
    }

  reallocate_particle_memory_MaxPart();

  do_box_wrapping();		/* map the particles back onto the box */

  toGo = mymalloc(sizeof(int) * NTask * NTask);
  local_toGo = mymalloc(sizeof(int) * NTask);
  list_NumPart = mymalloc(sizeof(int) * NTask);
  list_load = mymalloc(sizeof(int) * NTask);
  list_work = mymalloc(sizeof(double) * NTask);
  domain_allocatebuffer();

  MPI_Allgather(&NumPart, 1, MPI_INT, list_NumPart, 1, MPI_INT, MPI_COMM_WORLD);

  maxload = All.MaxPart * REDUC_FAC;

  domain_decomposeType(geom_flag);

  if(All.PM_Ti_endstep == All.Ti_Current)
    {
      for(i = 0; i < NumPart; i++)
	P[i].GravCost = 0.01;
    }

  domain_freebuffer();
  myfree(list_work);
  myfree(list_load);
  myfree(list_NumPart);
  myfree(local_toGo);
  myfree(toGo);

  reallocate_particle_memory_NumPart();

#ifdef HPM
  hpmStop(8);
#endif

  t1 = second();
  All.CPU_Domain += timediff(t0, t1);

  if(ThisTask == 0)
    {
      printf("domain decomposition done. (total %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }

  if(Flag_FullStep || geom_flag)	/* redo it only on full steps or before output */
    {
      t0 = second();
      peano_hilbert_order();
      t1 = second();
      All.CPU_Peano += timediff(t0, t1);
    }
}



/*  This function does the domain decomposition for a single (or several) particle type(s)
 */
void domain_decomposeType(int geom_flag)
{
  int i, status;
  long long sumtogo, sumload;
  int ngrp, task, partner, sendcount, recvcount;
  int maxload;
  double sumwork, maxwork;
  double t0, t1;

  /* determine global dimensions of domain grid */
  domain_findExtent();

  /* determine cost distribution in domain grid */
  domain_sumCost();

  startlist = &DomainStartList[0];
  endlist = &DomainEndList[0];

  if(geom_flag == 0)
    {
      /* find the split of the domain grid recursively */
      status = domain_findSplit(0, NTask, 0, DOMAINGRID * DOMAINGRID * DOMAINGRID - 1);
      if(status != 0)
	{
	  if(ThisTask == 0)
	    printf("\nNo domain decomposition that stays within memory bounds is possible.\n");
	  endrun(0);
	}

      /* now try to improve the work-load balance of the split */
      domain_shiftSplit();
    }
  else
    {
      status = domain_findSplit_equalvolume(0, NTask, 0, DOMAINGRID * DOMAINGRID * DOMAINGRID - 1);
      if(status != 0)
	{
	  /* now fall back to equal particle load */
	  status = domain_findSplit(0, NTask, 0, DOMAINGRID * DOMAINGRID * DOMAINGRID - 1);
	  if(status != 0)
	    {
	      if(ThisTask == 0)
		printf("\nNo domain decomposition that stays within memory bounds is possible.\n");
	      endrun(0);
	    }
	}
    }

  DomainMyStart = startlist[ThisTask];
  DomainMyLast = endlist[ThisTask];

  if(ThisTask == 0)
    {
      sumload = maxload = 0;
      sumwork = maxwork = 0;
      for(i = 0; i < NTask; i++)
	{
	  sumload += list_load[i];
	  sumwork += list_work[i];

	  if(list_load[i] > maxload)
	    maxload = list_load[i];

	  if(list_work[i] > maxwork)
	    maxwork = list_work[i];
	}

      printf("work-load balance=%g   memory-balance=%g\n",
	     maxwork / (sumwork / NTask), maxload / (((double) sumload) / NTask));
      fflush(stdout);
    }

  /* determine for each cpu how many particles have to be shifted to other cpus */
  domain_countToGo();

  for(i = 0, sumtogo = 0; i < NTask * NTask; i++)
    sumtogo += toGo[i];

  while(sumtogo > 0)
    {
      if(ThisTask == 0)
	{
	  printf("exchange of %d%09d particles ", (int) (sumtogo / 1000000000), (int) (sumtogo % 1000000000));
	  fflush(stdout);
	}

      t0 = second();

      for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	{
	  for(task = 0; task < NTask; task++)
	    {
	      partner = task ^ ngrp;

	      if(partner < NTask && task < partner)
		{
		  domain_findExchangeNumbers(task, partner, &sendcount, &recvcount);

		  list_NumPart[task] += recvcount - sendcount;
		  list_NumPart[partner] -= recvcount - sendcount;

		  toGo[task * NTask + partner] -= sendcount;
		  toGo[partner * NTask + task] -= recvcount;

		  if(task == ThisTask)	/* actually carry out the exchange */
		    domain_exchangeParticles(partner, sendcount, recvcount);
		  if(partner == ThisTask)
		    domain_exchangeParticles(task, recvcount, sendcount);
		}
	    }
	}

      t1 = second();
      if(ThisTask == 0)
	{
	  printf("  (took %g sec)\n", timediff(t0, t1));
	  fflush(stdout);
	}

      for(i = 0, sumtogo = 0; i < NTask * NTask; i++)
	sumtogo += toGo[i];
    }

}

/*! This function tries to find a split point in a range of cells in the domaingrid.  The range of cells
 *  starts at 'first', and ends at 'last' (inclusively). The number of cpus that holds the range is 'ncpu',
 *  with the first cpu given by 'cpustart'. If more than 2 cpus are to be split, the function calls itself
 *  recursively. The division tries to achieve a best particle-load balance under the constraint that
 *  'maxload' and 'maxloadsph' may not be exceeded, and that each cpu holds at least one cell from the
 *  domaingrid. If such a decomposition cannot be achieved, a non-zero error code is returned.
 *
 *  After successful completion, DomainMyStart[] and DomainMyLast[] contain the first and last cell of the
 *  domaingrid assigned to the local task for the given type. Also, DomainTask[] contains for each cell the
 *  task it was assigned to.
 */
int domain_findSplit(int cpustart, int ncpu, int first, int last)
{
  int i, split, ok_left, ok_right;
  long long load, load_leftOfSplit;
  int ncpu_leftOfSplit;
  int *loadtab, *tasktab;
  double maxAvgLoad_CurrentSplit, maxAvgLoad_NewSplit;

  loadtab = &DomainCount[0];
  tasktab = &DomainTask[0];

  ncpu_leftOfSplit = ncpu / 2;

  for(i = first, load = 0; i <= last; i++)
    load += loadtab[i];

  split = first + ncpu_leftOfSplit;

  for(i = first, load_leftOfSplit = 0; i < split; i++)
    load_leftOfSplit += loadtab[i];

  /* find the best split point in terms of work-load balance */

  while(split < last - (ncpu - ncpu_leftOfSplit - 1))
    {
      maxAvgLoad_CurrentSplit =
	dmax(load_leftOfSplit / ncpu_leftOfSplit, (load - load_leftOfSplit) / (ncpu - ncpu_leftOfSplit));

      maxAvgLoad_NewSplit =
	dmax((load_leftOfSplit + loadtab[split]) / ncpu_leftOfSplit,
	     (load - load_leftOfSplit - loadtab[split]) / (ncpu - ncpu_leftOfSplit));

      if(maxAvgLoad_NewSplit <= maxAvgLoad_CurrentSplit)
	{
	  load_leftOfSplit += loadtab[split];
	  split++;
	}
      else
	break;
    }


  /* we will now have to check whether this solution is possible given the restrictions on the maximum load */

  for(i = first, load_leftOfSplit = 0; i < split; i++)
    load_leftOfSplit += loadtab[i];

  if(((double) load_leftOfSplit) / ncpu_leftOfSplit > maxload ||
     ((double) (load - load_leftOfSplit)) / (ncpu - ncpu_leftOfSplit) > maxload)
    {
      /* we did not find a viable split */
      return -1;
    }

  if(ncpu_leftOfSplit >= 2)
    ok_left = domain_findSplit(cpustart, ncpu_leftOfSplit, first, split - 1);
  else
    ok_left = 0;

  if((ncpu - ncpu_leftOfSplit) >= 2)
    ok_right = domain_findSplit(cpustart + ncpu_leftOfSplit, ncpu - ncpu_leftOfSplit, split, last);
  else
    ok_right = 0;

  if(ok_left == 0 && ok_right == 0)
    {
      /* found a viable split */

      if(ncpu_leftOfSplit == 1)
	{
	  for(i = first; i < split; i++)
	    tasktab[i] = cpustart;

	  list_load[cpustart] = load_leftOfSplit;
	  startlist[cpustart] = first;
	  endlist[cpustart] = split - 1;
	}

      if((ncpu - ncpu_leftOfSplit) == 1)
	{
	  for(i = split; i <= last; i++)
	    tasktab[i] = cpustart + ncpu_leftOfSplit;

	  list_load[cpustart + ncpu_leftOfSplit] = load - load_leftOfSplit;
	  startlist[cpustart + ncpu_leftOfSplit] = split;
	  endlist[cpustart + ncpu_leftOfSplit] = last;
	}
      
      return 0;
    }

  /* we did not find a viable split */
  return -1;
}


int domain_findSplit_equalvolume(int cpustart, int ncpu, int first, int last)
{
  int i, split, ok_left, ok_right, task;
  long long load, load_leftOfSplit;
  int ncpu_leftOfSplit;
  int *loadtab, *tasktab;

  loadtab = &DomainCount[0];
  tasktab = &DomainTask[0];

  ncpu_leftOfSplit = ncpu / 2;

  for(i = first, load = 0; i <= last; i++)
    load += loadtab[i];


  split = first + ncpu_leftOfSplit * ((last - first + 1) / ncpu);

  

  /* we now check whether this solution is possible given the restrictions on the maximum load */

  for(i = first, load_leftOfSplit = 0; i < split; i++)
    load_leftOfSplit += loadtab[i];

  if(((double) load_leftOfSplit) / ncpu_leftOfSplit > maxload ||
     ((double) (load - load_leftOfSplit)) / (ncpu - ncpu_leftOfSplit) > maxload)
    {
      /* we did not find a viable split */
      return -1;
    }

  if(ncpu_leftOfSplit >= 2)
    ok_left = domain_findSplit_equalvolume(cpustart, ncpu_leftOfSplit, first, split - 1);
  else
    ok_left = 0;

  if((ncpu - ncpu_leftOfSplit) >= 2)
    ok_right =
      domain_findSplit_equalvolume(cpustart + ncpu_leftOfSplit, ncpu - ncpu_leftOfSplit, split, last);
  else
    ok_right = 0;

  if(ok_left == 0 && ok_right == 0)
    {
      /* found a viable split */

      if(ncpu_leftOfSplit == 1)
	{
	  for(i = first; i < split; i++)
	    tasktab[i] = cpustart;

	  list_load[cpustart] = load_leftOfSplit;
	  startlist[cpustart] = first;
	  endlist[cpustart] = split - 1;
	}

      if((ncpu - ncpu_leftOfSplit) == 1)
	{
	  for(i = split; i <= last; i++)
	    tasktab[i] = cpustart + ncpu_leftOfSplit;

	  list_load[cpustart + ncpu_leftOfSplit] = load - load_leftOfSplit;
	  startlist[cpustart + ncpu_leftOfSplit] = split;
	  endlist[cpustart + ncpu_leftOfSplit] = last;
	}

      for(task = 0; task < NTask; task++)
	list_work[task] = 0;

      for(i = 0; i < DOMAINGRID * DOMAINGRID * DOMAINGRID; i++)
	list_work[DomainTask[i]] += DomainWork[i];
      
      return 0;
    }

  /* we did not find a viable split */
  return -1;
}




void domain_shiftSplit(void)
{
  int i, task, iter = 0, moved;
  double maxw, newmaxw;

  for(task = 0; task < NTask; task++)
    list_work[task] = 0;

  for(i = 0; i < DOMAINGRID * DOMAINGRID * DOMAINGRID; i++)
    list_work[DomainTask[i]] += DomainWork[i];

  do
    {
      for(task = 0, moved = 0; task < NTask - 1; task++)
	{
	  maxw = dmax(list_work[task], list_work[task + 1]);

	  if(list_work[task] < list_work[task + 1])
	    {
	      newmaxw = dmax(list_work[task] + DomainWork[startlist[task + 1]],
			     list_work[task + 1] - DomainWork[startlist[task + 1]]);
	      if(newmaxw < maxw)
		{
		  if(list_load[task] + DomainCount[startlist[task + 1]] < maxload)
		    {
		      /* ok, we can move one domain cell from right to left */
		      list_work[task] += DomainWork[startlist[task + 1]];
		      list_load[task] += DomainCount[startlist[task + 1]];
		      list_work[task + 1] -= DomainWork[startlist[task + 1]];
		      list_load[task + 1] -= DomainCount[startlist[task + 1]];

		      DomainTask[startlist[task + 1]] = task;
		      startlist[task + 1] += 1;
		      endlist[task] += 1;

		      moved++;
		    }
		}
	    }
	  else
	    {
	      newmaxw = dmax(list_work[task] - DomainWork[endlist[task]],
			     list_work[task + 1] + DomainWork[endlist[task]]);
	      if(newmaxw < maxw)
		{
		  if(list_load[task + 1] + DomainCount[endlist[task]] < maxload)
		    {
		      /* ok, we can move one domain cell from left to right */
		      list_work[task] -= DomainWork[endlist[task]];
		      list_load[task] -= DomainCount[endlist[task]];
		      list_work[task + 1] += DomainWork[endlist[task]];
		      list_load[task + 1] += DomainCount[endlist[task]];

		      DomainTask[endlist[task]] = task + 1;
		      endlist[task] -= 1;
		      startlist[task + 1] -= 1;

		      moved++;
		    }
		}

	    }
	}

      iter++;
    }
  while(moved > 0 && iter < DOMAINGRID * DOMAINGRID * DOMAINGRID);
}



void domain_findExchangeNumbers(int task, int partner, int *send, int *recv)
{
  int numpartA, ntobesentA, maxsendA, maxsendA_old;
  int numpartB, ntobesentB, maxsendB, maxsendB_old;

  numpartA = list_NumPart[task];
  numpartB = list_NumPart[partner];

  ntobesentA = toGo[task * NTask + partner];
  ntobesentB = toGo[partner * NTask + task];

  maxsendA = imin(ntobesentA, All.BunchSizeDomain);
  maxsendB = imin(ntobesentB, All.BunchSizeDomain);

  do
    {
      maxsendA_old = maxsendA;
      maxsendB_old = maxsendB;

      maxsendA = imin(All.MaxPart - numpartB + ntobesentB, maxsendA);
      maxsendB = imin(All.MaxPart - numpartA + ntobesentA, maxsendB);
    }
  while((maxsendA != maxsendA_old) || (maxsendB != maxsendB_old));

  *send = maxsendA;
  *recv = maxsendB;
}






/*  Ok, we have settled on a split. This function tries to get rid of all its particles that are to the RIGHT
 *  of the split. It communicates with another processor which runs exchangeParticles_B(). A is sending first
 *  to B, then it is receiving from it.  The number of particles exchanged is negotiated such that the
 *  transfer is possible within the given memory constraints.
 */
void domain_exchangeParticles(int partner, int send_count, int recv_count)
{
  int i, j, k, n, count, rep, index;
  double xx, yy, zz;
  MPI_Status status;
  double scalefac;

  scalefac = 1.0 / All.BoxSize;

  for(n = 0, count = 0; count < send_count && n < NumPart; n++)
    {
      xx = P[n].Pos[0] * scalefac + 1.0;
      yy = P[n].Pos[1] * scalefac + 1.0;
      zz = P[n].Pos[2] * scalefac + 1.0;

      i = DOUBLE_to_DOMAINGRID(xx);
      j = DOUBLE_to_DOMAINGRID(yy);
      k = DOUBLE_to_DOMAINGRID(zz);

      index = DomainPeanoMap[(i * DOMAINGRID + j) * DOMAINGRID + k];


      if(DomainTask[index] == partner)
	{
	  DomainPartBuf[count] = P[n];	/* copy particle and collect in contiguous memory */
	  P[n] = P[NumPart - 1];

	  count++;
	  NumPart--;
	  n--;
	}
    }

  if(count != send_count)
    {
      printf("Houston, we got a problem...\n");
      printf("ThisTask=%d count=%d send_count=%d\n", ThisTask, count, send_count);
      endrun(888);
    }

  /* transmit */
  for(rep = 0; rep < 2; rep++)
    {
      if((rep == 0 && ThisTask < partner) || (rep == 1 && ThisTask > partner))
	{
	  if(send_count > 0)
	    {
	      MPI_Ssend(&DomainPartBuf[0], send_count * sizeof(struct particle_data), MPI_BYTE, partner,
			TAG_PDATA, MPI_COMM_WORLD);
	    }
	}

      if((rep == 1 && ThisTask < partner) || (rep == 0 && ThisTask > partner))
	{
	  if(recv_count > 0)
	    {

	      MPI_Recv(&P[NumPart], recv_count * sizeof(struct particle_data), MPI_BYTE, partner,
		       TAG_PDATA, MPI_COMM_WORLD, &status);

	      NumPart += recv_count;
	    }
	}
    }
}



void domain_countToGo(void)
{
  int n, i, j, k, index;
  double xx, yy, zz;
  double scalefac;

  scalefac = 1.0 / All.BoxSize;

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      xx = P[n].Pos[0] * scalefac + 1.0;
      yy = P[n].Pos[1] * scalefac + 1.0;
      zz = P[n].Pos[2] * scalefac + 1.0;

      i = DOUBLE_to_DOMAINGRID(xx);
      j = DOUBLE_to_DOMAINGRID(yy);
      k = DOUBLE_to_DOMAINGRID(zz);

      index = DomainPeanoMap[(i * DOMAINGRID + j) * DOMAINGRID + k];

      if(DomainTask[index] != ThisTask)
	local_toGo[DomainTask[index]] += 1;
    }

  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);
}

void domain_sumCost(void)
{
  int i, j, k, n, index;
  int local_DomainCount[DOMAINGRID * DOMAINGRID * DOMAINGRID];
  double local_DomainWork[DOMAINGRID * DOMAINGRID * DOMAINGRID];
  double xx, yy, zz;
  double scalefac;

  scalefac = 1.0 / All.BoxSize;

  for(i = 0; i < DOMAINGRID * DOMAINGRID * DOMAINGRID; i++)
    {
      local_DomainWork[i] = 0;
      local_DomainCount[i] = 0;
    }

  for(n = 0; n < NumPart; n++)
    {
      xx = P[n].Pos[0] * scalefac + 1.0;
      yy = P[n].Pos[1] * scalefac + 1.0;
      zz = P[n].Pos[2] * scalefac + 1.0;
      
      i = DOUBLE_to_DOMAINGRID(xx);
      j = DOUBLE_to_DOMAINGRID(yy);
      k = DOUBLE_to_DOMAINGRID(zz);
      
      index = DomainPeanoMap[(i * DOMAINGRID + j) * DOMAINGRID + k];

      local_DomainWork[index] += P[n].GravCost;
      local_DomainCount[index] += 1;
    }
  
  MPI_Allreduce(local_DomainWork, DomainWork, DOMAINGRID * DOMAINGRID * DOMAINGRID, MPI_DOUBLE, MPI_SUM,
		MPI_COMM_WORLD);
  MPI_Allreduce(local_DomainCount, DomainCount, DOMAINGRID * DOMAINGRID * DOMAINGRID, MPI_INT, MPI_SUM,
		MPI_COMM_WORLD);
}


/*! This routine finds the extent of the global domain grid.
 */
void domain_findExtent(void)
{
  int j;

  for(j = 0; j < 3; j++)
    DomainCenter[j] = 0.5 * All.BoxSize;

  DomainLen = All.BoxSize;
  DomainFac = DOMAINGRID / All.BoxSize;
}




void domain_allocatebuffer(void)
{
  size_t bytes;

  /* we use here the same amount of memory that we allocate for the BH
     tree, that's quite a bit of space, normally more than the
     communication buffer has */
  bytes = All.TreeAllocFactor * (All.MaxPart / All.PartAllocFactor) * sizeof(struct NODE);

  if(!(CommBuffer = mymalloc(bytes)))
    {
      printf("failed to allocate memory for `CommBuffer' in domain decomposition (%g MB).\n",
	     bytes / (1024.0 * 1024.0));
      endrun(2);
    }

  All.BunchSizeDomain = bytes / (sizeof(struct particle_data));

  DomainPartBuf = (struct particle_data *) CommBuffer;

  if(ThisTask == 0 && All.NumCurrentTiStep == 0)
    {
      printf("Communication buffer in domain decomposition has room for %d particles (%g MByte)\n",
	     All.BunchSizeDomain, bytes / (1024.0 * 1024.0));
      fflush(stdout);
    }
}

void domain_freebuffer(void)
{
  myfree(CommBuffer);
}
