#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

/*! \file twopoint.c
 *  \brief computes the two-point mass correlation function on the fly 
 */


#define BINS_TP  40		/* number of bins used */
#define ALPHA  -1.0		/* slope used in randomly selecting radii around target particles */

#ifndef FRACTION_TP
#define FRACTION_TP  0.2	
#endif                             /* fraction of particles selected for sphere
				   placement. Will be scaled with total
				   particle number so that a fixed value
				   should give roughly the same noise level
				   in the meaurement, indpendent of
				   simulation size */



/*! this macro maps a coordinate difference to the nearest periodic image
 */
#define NGB_PERIODIC(x) (((x)>BoxHalf)?((x)-BoxSize):(((x)<-BoxHalf)?((x)+BoxSize):(x)))
#define SQUARE_IT(x) ((x)*(x))



static long long Count[BINS_TP];
static long long CountSpheres[BINS_TP];
static double Xi[BINS_TP];
static double Rbin[BINS_TP];

static double R0, R1;		/* inner and outer radius for correlation function determination */

static double logR0;
static double binfac;

static double BoxHalf, BoxSize;


/*  This function computes the two-point function.
 */
void twopoint(void)
{
  int i, j, bin, n;
  double p, rs, vol, scaled_frac;
  int *noffset, *nbuffer, *nsend_local, *nsend, *ndonelist;
  int ndone, nexport;
  long long ntotleft, *countbuf;
  int maxfill, sendTask, recvTask, level, ngrp;
  double tstart, tend, t0, t1;
  void *state_buffer;
  MPI_Status status;

  if(ThisTask == 0)
    {
      printf("begin two-point correlation function...\n");
      fflush(stdout);
    }

  tstart = second();


  /* set inner and outer radius for the bins that are used for the correlation function estimate */
  R0 = All.Softening;
  R1 = All.BoxSize / 2;


  scaled_frac = FRACTION_TP * 1.0e7 / All.TotNumPart;

  logR0 = log(R0);
  binfac = BINS_TP / (log(R1) - log(R0));


  BoxSize = All.BoxSize;
  BoxHalf = 0.5 * All.BoxSize;


  for(i = 0; i < BINS_TP; i++)
    {
      Count[i] = 0;
      CountSpheres[i] = 0;
    }


  allocate_commbuffers();

  noffset = mymalloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = mymalloc(sizeof(int) * NTask);
  nsend_local = mymalloc(sizeof(int) * NTask);
  nsend = mymalloc(sizeof(int) * NTask * NTask);
  ndonelist = mymalloc(sizeof(int) * NTask);

  state_buffer = mymalloc(gsl_rng_size(random_generator));

  memcpy(state_buffer, gsl_rng_state(random_generator), gsl_rng_size(random_generator));

  gsl_rng_set(random_generator, P[0].ID + ThisTask);	/* seed things with first particle ID to make sure we are
							   different on each CPU */


  i = 0;			/* beginn with this index */

  ntotleft = All.TotNumPart;	/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */

      t0 = second();

      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeTwoPoint - NTask; i++, ndone++)
	if(gsl_rng_uniform(random_generator) < scaled_frac)
	  {
	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;

	    p = gsl_rng_uniform(random_generator);

	    rs = pow(pow(R0, ALPHA) + p * (pow(R1, ALPHA) - pow(R0, ALPHA)), 1 / ALPHA);

	    bin = (log(rs) - logR0) * binfac;

	    rs = exp((bin + 1) / binfac + logR0);

	    for(j = 0; j <= bin; j++)
	      CountSpheres[j]++;

	    count_local(i, 0, rs);

	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    TwoPointDataIn[nexport].Pos[0] = P[i].Pos[0];
		    TwoPointDataIn[nexport].Pos[1] = P[i].Pos[1];
		    TwoPointDataIn[nexport].Pos[2] = P[i].Pos[2];
		    TwoPointDataIn[nexport].Rs = rs;
		    TwoPointDataIn[nexport].Task = j;
		    nexport++;
		    nsend_local[j]++;
		  }
	      }
	  }

      t1 = second();
      if(ThisTask == 0)
	{
	  printf("have done %d local particles in %g on root-task\n", ndone, timediff(t0, t1));
	  fflush(stdout);
	}

      qsort(TwoPointDataIn, nexport, sizeof(struct twopointdata_in), twopoint_compare_key);

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);


      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
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
	      if(maxfill >= All.BunchSizeTwoPoint)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&TwoPointDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct twopointdata_in), MPI_BYTE,
				   recvTask, TAG_TWOPOINT,
				   &TwoPointDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct twopointdata_in),
				   MPI_BYTE, recvTask, TAG_TWOPOINT, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	  t0 = second();

	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    count_local(j, 1, 0);

	  t1 = second();
	  if(ThisTask == 0)
	    {
	      printf("have done %d imported  particles in %g on root-task\n", nbuffer[ThisTask], timediff(t0, t1));
	      fflush(stdout);
	    }

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  memcpy(gsl_rng_state(random_generator), state_buffer, gsl_rng_size(random_generator));
  myfree(state_buffer);

  myfree(ndonelist);
  myfree(nsend);
  myfree(nsend_local);
  myfree(nbuffer);
  myfree(noffset);


  free_commbuffers();




  /* Now compute the actual correlation function */

  countbuf = mymalloc(NTask * BINS_TP * sizeof(long long));

  MPI_Allgather(Count, BINS_TP * sizeof(long long), MPI_BYTE,
		countbuf, BINS_TP * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_TP; i++)
    {
      Count[i] = 0;
      for(n = 0; n < NTask; n++)
	Count[i] += countbuf[n * BINS_TP + i];
    }

  MPI_Allgather(CountSpheres, BINS_TP * sizeof(long long), MPI_BYTE,
		countbuf, BINS_TP * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_TP; i++)
    {
      CountSpheres[i] = 0;
      for(n = 0; n < NTask; n++)
	CountSpheres[i] += countbuf[n * BINS_TP + i];
    }

  myfree(countbuf);


  for(i = 0; i < BINS_TP; i++)
    {
      vol = 4 * M_PI / 3.0 * (pow(exp((i + 1.0) / binfac + logR0), 3)
			      - pow(exp((i + 0.0) / binfac + logR0), 3));

      if(CountSpheres[i] > 0)
	Xi[i] = -1 + Count[i] / ((double) CountSpheres[i]) / (All.TotNumPart / pow(All.BoxSize, 3)) / vol;
      else
	Xi[i] = 0;

      Rbin[i] = exp((i + 0.5) / binfac + logR0);
    }

  twopoint_save();

  tend = second();
  All.CPU_TwoPoint += timediff(tstart, tend);

  if(ThisTask == 0)
    {
      printf("end two-point: Took=%g seconds.\n", timediff(tstart, tend));
      fflush(stdout);
    }
}




void twopoint_save(void)
{
  FILE *fd;
  char buf[500];
  int i;

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/correl_%03d.txt", All.OutputDir, All.TwoPointFlag - 1);

      if(!(fd = fopen(buf, "w")))
	{
	  printf("can't open file `%s`\n", buf);
	  endrun(1323);
	}

      fprintf(fd, "%g\n", All.Time);
      i = BINS_TP;
      fprintf(fd, "%d\n", i);

      for(i = 0; i < BINS_TP; i++)
	fprintf(fd, "%g %g %g %g\n", Rbin[i], Xi[i], (double) Count[i], (double) CountSpheres[i]);

      fclose(fd);
    }
}




/*! This function counts the pairs in a sphere
 */
void count_local(int target, int mode, double rs)
{
  FLOAT pos[3];

  if(mode == 0)
    {
      pos[0] = P[target].Pos[0];
      pos[1] = P[target].Pos[1];
      pos[2] = P[target].Pos[2];
    }
  else
    {
      pos[0] = TwoPointDataGet[target].Pos[0];
      pos[1] = TwoPointDataGet[target].Pos[1];
      pos[2] = TwoPointDataGet[target].Pos[2];
      rs = TwoPointDataGet[target].Rs;
    }

  twopoint_ngb_treefind_variable(&pos[0], rs, mode);
}





/*! This function finds all particles within the radius "rsearch",
 *  and counts them in the bins used for the two-point correlation function.
 */
void twopoint_ngb_treefind_variable(FLOAT searchcenter[3], FLOAT rsearch, int mode)
{
  FLOAT searchmin[3], searchmax[3];
  double r2, r, ri, ro;
  int k, no, p, bin, bin2;
  struct NODE *nop;

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - rsearch;
      searchmax[k] = searchcenter[k] + rsearch;
    }


  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  r2 = SQUARE_IT(NGB_PERIODIC(P[p].Pos[0] - searchcenter[0])) +
	    SQUARE_IT(NGB_PERIODIC(P[p].Pos[1] - searchcenter[1])) +
	    SQUARE_IT(NGB_PERIODIC(P[p].Pos[2] - searchcenter[2]));


	  if(r2 >= R0 * R0 && r2 < R1 * R1)
	    {
	      if(r2 < rsearch * rsearch)
		{
		  bin = (log(sqrt(r2)) - logR0) * binfac;
		  if(bin < BINS_TP)
		    Count[bin]++;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;

	      no = DomainNextnode[no - (All.MaxPart + MaxNodes)];
	      continue;
	    }

	  nop = &Nodes[no];

	  if(mode == 1)
	    {
	      if((nop->u.d.cost & 3) == 1)	/* if it's a top-level node which does not contain local particles */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  no = nop->u.d.sibling;	/* make skipping the branch the default */

	  if((NGB_PERIODIC(nop->center[0] - searchcenter[0]) + 0.5 * nop->len) <
	     (searchmin[0] - searchcenter[0]))
	    continue;
	  if((NGB_PERIODIC(nop->center[0] - searchcenter[0]) - 0.5 * nop->len) >
	     (searchmax[0] - searchcenter[0]))
	    continue;
	  if((NGB_PERIODIC(nop->center[1] - searchcenter[1]) + 0.5 * nop->len) <
	     (searchmin[1] - searchcenter[1]))
	    continue;
	  if((NGB_PERIODIC(nop->center[1] - searchcenter[1]) - 0.5 * nop->len) >
	     (searchmax[1] - searchcenter[1]))
	    continue;
	  if((NGB_PERIODIC(nop->center[2] - searchcenter[2]) + 0.5 * nop->len) <
	     (searchmin[2] - searchcenter[2]))
	    continue;
	  if((NGB_PERIODIC(nop->center[2] - searchcenter[2]) - 0.5 * nop->len) >
	     (searchmax[2] - searchcenter[2]))
	    continue;


	  r2 = SQUARE_IT(NGB_PERIODIC(nop->center[0] - searchcenter[0])) +
	    SQUARE_IT(NGB_PERIODIC(nop->center[1] - searchcenter[1])) +
	    SQUARE_IT(NGB_PERIODIC(nop->center[2] - searchcenter[2]));

	  r = sqrt(r2);

	  ri = r - 0.5 * 1.732051 * nop->len;
	  ro = r + 0.5 * 1.732051 * nop->len;

	  if(ri >= R0 && ro < R1)
	    {
	      if(ro < rsearch)
		{
		  bin = (log(ri) - logR0) * binfac;
		  bin2 = (log(ro) - logR0) * binfac;
		  if(bin == bin2)
		    {
		      if(mode == 1)
			{
			  if((nop->u.d.cost & 1))	/* Bit 0 signals that this node belongs to top-level tree */
			    continue;
			}
		      Count[bin] += nop->u.d.mass;
		      continue;
		    }
		}
	    }

	  no = nop->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

}



/*! This routine is a comparison kernel used in a sort routine to group
 *  particles that are exported to the same processor.
 */
int twopoint_compare_key(const void *a, const void *b)
{
  if(((struct twopointdata_in *) a)->Task < (((struct twopointdata_in *) b)->Task))
    return -1;

  if(((struct twopointdata_in *) a)->Task > (((struct twopointdata_in *) b)->Task))
    return +1;

  return 0;
}
