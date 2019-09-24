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

/*! \file fof.c
 *  \brief parallel FoF group finder
 */

#define DOUBLE_to_DOMAINGRID(y)  ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - DOMAINLEVELS)))
#define DOUBLE_to_HASHBITS(y)    ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - HASHBITS)))

#define FOF_PERIODIC(x) (xtmp=(x),(xtmp>boxHalf)?(xtmp-boxSize):((xtmp<-boxHalf)?(xtmp+boxSize):xtmp))
#define FOF_PERIODIC_WRAP(x) (xtmp=(x),(xtmp>=boxSize)?(xtmp-boxSize):((xtmp<0)?(xtmp+boxSize):xtmp))

static double boxHalf, boxSize;


static struct export_fof
{
  FLOAT Pos[3];
  long long GrID;
}
 *ExportBuf, *ImportBuf;

static struct id_list
{
  long long GrID;
  long long ID;
}
 *local_ids, *export_ids;


static double LinkL;
static int Ngroups, TotNgroups, Nids;
static long long TotNids;
static int Nghosts, Nexport, Nlocal;
static int Grid;

static int *local_toGo, *toGo, *export_offset, *import_offset;
static int *Len, *Head, *Next, *Tail, *GridNext;
static int *GroupLen, *GroupOffset;
static int *GridFirst;
static int *CountBelowMinLen;
static long long *GrID;
static long long *GroupIDs;
static char *GridFlag;


void fof_fof(int num)
{
  int i, links, largestgroup;
  double t0, t1, tstart, tend;

  t0 = second();

#ifdef HPM
  hpmStart(7, "FoF");
#endif

  if(ThisTask == 0)
    {
      printf("\nBegin to compute FoF group catalogues...\n");
      fflush(stdout);
    }

  LinkL = LINKLENGTH * All.BoxSize / pow(All.TotNumPart, 1.0 / 3);

  if(ThisTask == 0)
    printf("Comoving linking length: %g\n", LinkL);

  boxHalf = 0.5 * All.BoxSize;
  boxSize = All.BoxSize;

  local_toGo = mymalloc(sizeof(int) * NTask);
  toGo = mymalloc(sizeof(int) * NTask * NTask);
  export_offset = mymalloc(sizeof(int) * NTask);
  import_offset = mymalloc(sizeof(int) * NTask);
  CountBelowMinLen = mymalloc(sizeof(int) * GROUP_MIN_LEN);

  fof_import_ghosts();

  Head = mymalloc((NumPart + Nghosts) * sizeof(int));
  GrID = mymalloc((NumPart + Nghosts) * sizeof(long long));

  Len = mymalloc((NumPart + Nghosts) * sizeof(int));
  Next = mymalloc((NumPart + Nghosts) * sizeof(int));
  Tail = mymalloc((NumPart + Nghosts) * sizeof(int));
  GridNext = mymalloc((NumPart + Nghosts) * sizeof(int));

  if(!(GridNext))
    {
      printf("failed to allocate memory\n");
      endrun(1231);
    }

  for(i = 0; i < NumPart + Nghosts; i++)
    {
      Head[i] = Tail[i] = i;
      Len[i] = 1;
      Next[i] = -1;

      if(i < NumPart)
	GrID[i] = P[i].ID + (((long long) ThisTask) << 48);
      else
	GrID[i] = ImportBuf[i - NumPart].GrID;
    }

  tstart= second();
  fof_course_binning();
  tend= second();

  if(ThisTask == 0)
    {
      printf("Coarse binning done. (%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }

  tstart= second();
  fof_find_groups();
  tend= second();

  if(ThisTask == 0)
    {
      printf("Local groups found. (%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }

  tstart= second();
  fof_find_minids();

  do
    {
      links = fof_link_accross();

      if(ThisTask == 0)
	{
	  printf("Made %d links accross domain boundaries...\n", links);
	  fflush(stdout);
	}
    }
  while(links > 0);
  tend= second();

  if(ThisTask == 0)
    {
      printf("Groups linked accross processor boundaries. (%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }

  myfree(GridFlag);
  myfree(GridFirst);
  myfree(GridNext);
  myfree(Tail);
  myfree(Next);
  myfree(Len);

  tend= second();
  fof_exchange_id_lists();

  fof_compile_catalogue();
  tend= second();

  if(ThisTask == 0)
    {
      printf("ID-lists exchanged and localized (%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }


  MPI_Allreduce(&Ngroups, &TotNgroups, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allgather(&Nids, 1, MPI_INT, local_toGo, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, TotNids = 0; i < NTask; i++)
    TotNids += local_toGo[i];

  MPI_Allreduce(&GroupLen[0], &largestgroup, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  /* GroupLen[0] always exists */

  t1 = second();

  All.CPU_FoF += timediff(t0, t1);

  if(ThisTask == 0)
    {
      printf("Total number of groups with at least %d particles: %d\n", GROUP_MIN_LEN, TotNgroups);
      if(TotNgroups > 0)
	{
	  printf("Group catalogues determined. (took in total %g sec)\n", timediff(t0, t1));
	  printf("Largest group has %d particles.\n", largestgroup);
	  printf("Total number of particles in groups: %d%09d\n",
		 (int) (TotNids / 1000000000), (int) (TotNids % 1000000000));
	}
      fflush(stdout);
    }

#ifdef HPM
  hpmStop(7);
#endif

  t0 = second();

  fof_save_groups(num);

  t1 = second();
  All.CPU_Snapshot += timediff(t0, t1);

  if(ThisTask == 0)
    {
      printf("Group catalogues saved. (I/O took %g sec)\n\n", timediff(t0, t1));
      fflush(stdout);
    }

  myfree(GroupOffset);
  myfree(GroupLen);
  myfree(local_ids);
  myfree(GrID);
  myfree(Head);
  myfree(ImportBuf);
  myfree(ExportBuf);
  myfree(CountBelowMinLen);
  myfree(import_offset);
  myfree(export_offset);
  myfree(toGo);
  myfree(local_toGo);
}



void fof_import_ghosts(void)
{
  int n, ix, iy, iz, i, j, k, index;
  int level, sendTask, recvTask;
  double x, y, z, xtmp;
  double scalefac;
  MPI_Status status;

  scalefac = 1.0 / All.BoxSize;

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      for(j = 0; j < NTask; j++)
	Exportflag[j] = 0;

      for(ix = -1; ix <= 1; ix++)
	for(iy = -1; iy <= 1; iy++)
	  for(iz = -1; iz <= 1; iz++)
	    if(ix != 0 || iy != 0 || iz != 0)
	      {
		x = FOF_PERIODIC_WRAP(P[n].Pos[0] + ix * LinkL) * scalefac + 1.0;
		y = FOF_PERIODIC_WRAP(P[n].Pos[1] + iy * LinkL) * scalefac + 1.0;
		z = FOF_PERIODIC_WRAP(P[n].Pos[2] + iz * LinkL) * scalefac + 1.0;

		i = DOUBLE_to_DOMAINGRID(x);
		j = DOUBLE_to_DOMAINGRID(y);
		k = DOUBLE_to_DOMAINGRID(z);

		index = DomainPeanoMap[(i * DOMAINGRID + j) * DOMAINGRID + k];

		if(DomainTask[index] != ThisTask)
		  Exportflag[DomainTask[index]] = 1;
	      }

      for(j = 0; j < NTask; j++)
	if(Exportflag[j])
	  local_toGo[j] += 1;
    }

  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, Nghosts = 0, Nexport = 0; j < NTask; j++)
    {
      Nghosts += toGo[j * NTask + ThisTask];
      Nexport += toGo[ThisTask * NTask + j];
    }

  printf("Task=%d needs to import %d ghost particles for FoF, and exports %d particles.\n", ThisTask, Nghosts,
	 Nexport);
  fflush(stdout);

  ExportBuf = mymalloc(1 + Nexport * sizeof(struct export_fof));
  ImportBuf = mymalloc(1 + Nghosts * sizeof(struct export_fof));

  if(!(ImportBuf) || !(ExportBuf))
    {
      printf("failed to allocate memory\n");
      endrun(12311);
    }

  export_offset[0] = 0;
  import_offset[0] = 0;

  for(n = 1; n < NTask; n++)
    {
      export_offset[n] = export_offset[n - 1] + toGo[ThisTask * NTask + (n - 1)];
      import_offset[n] = import_offset[n - 1] + toGo[(n - 1) * NTask + ThisTask];
    }


  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      for(j = 0; j < NTask; j++)
	Exportflag[j] = 0;

      for(ix = -1; ix <= 1; ix++)
	for(iy = -1; iy <= 1; iy++)
	  for(iz = -1; iz <= 1; iz++)
	    if(ix != 0 || iy != 0 || iz != 0)
	      {
		x = FOF_PERIODIC_WRAP(P[n].Pos[0] + ix * LinkL) * scalefac + 1.0;
		y = FOF_PERIODIC_WRAP(P[n].Pos[1] + iy * LinkL) * scalefac + 1.0;
		z = FOF_PERIODIC_WRAP(P[n].Pos[2] + iz * LinkL) * scalefac + 1.0;

		i = DOUBLE_to_DOMAINGRID(x);
		j = DOUBLE_to_DOMAINGRID(y);
		k = DOUBLE_to_DOMAINGRID(z);

		index = DomainPeanoMap[(i * DOMAINGRID + j) * DOMAINGRID + k];

		if(DomainTask[index] != ThisTask)
		  Exportflag[DomainTask[index]] = 1;
	      }

      for(j = 0; j < NTask; j++)
	if(Exportflag[j])
	  {
	    ExportBuf[export_offset[j] + local_toGo[j]].GrID = P[n].ID + (((long long) ThisTask) << 48);

	    for(k = 0; k < 3; k++)
	      ExportBuf[export_offset[j] + local_toGo[j]].Pos[k] = P[n].Pos[k];

	    local_toGo[j] += 1;
	  }
    }


  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&ExportBuf[export_offset[recvTask]],
		     toGo[ThisTask * NTask + recvTask] * sizeof(struct export_fof),
		     MPI_BYTE, recvTask, TAG_FOF_IMPORT,
		     &ImportBuf[import_offset[recvTask]],
		     toGo[recvTask * NTask + ThisTask] * sizeof(struct export_fof),
		     MPI_BYTE, recvTask, TAG_FOF_IMPORT, MPI_COMM_WORLD, &status);
    }

}


void fof_course_binning(void)
{
  int i, j, k, n;
  double fac;
  double xmin[3], xmax[3], len;


  /* determine extension */
  for(j = 0; j < 3; j++)
    {
      xmin[j] = MAX_REAL_NUMBER;
      xmax[j] = -MAX_REAL_NUMBER;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  if(xmin[j] > P[i].Pos[j])
	    xmin[j] = P[i].Pos[j];

	  if(xmax[j] < P[i].Pos[j])
	    xmax[j] = P[i].Pos[j];
	}
    }

  for(i = 0; i < Nghosts; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  if(ImportBuf[i].Pos[j] < xmin[j])
	    {
	      if((ImportBuf[i].Pos[j] + All.BoxSize - xmax[j]) < (xmin[j] - ImportBuf[i].Pos[j]))
		ImportBuf[i].Pos[j] += All.BoxSize;
	    }
	  else if(ImportBuf[i].Pos[j] > xmax[j])
	    {
	      if((xmin[j] - (ImportBuf[i].Pos[j] - All.BoxSize)) < (ImportBuf[i].Pos[j] - xmax[j]))
		ImportBuf[i].Pos[j] -= All.BoxSize;
	    }
	}
    }

  for(i = 0; i < Nghosts; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  if(xmin[j] > ImportBuf[i].Pos[j])
	    xmin[j] = ImportBuf[i].Pos[j];

	  if(xmax[j] < ImportBuf[i].Pos[j])
	    xmax[j] = ImportBuf[i].Pos[j];
	}
    }

  len = 0;
  for(j = 0; j < 3; j++)
    if(xmax[j] - xmin[j] > len)
      len = xmax[j] - xmin[j];

  /* determine grid size */

  Grid = pow(NumPart, 1.0 / 3);

  /* we need to make sure that one mesh cell is at least of thickness LinkL */

  if(len / Grid < LinkL)
    Grid = len / LinkL;

  fac = Grid / len;


  GridFirst = mymalloc(Grid * Grid * Grid * sizeof(int));
  GridFlag = mymalloc(Grid * Grid * Grid * sizeof(char));

  if(!(GridFirst) || !(GridFlag))
    {
      printf("failed to allocate memory\n");
      endrun(112311);
    }

  for(n = 0; n < Grid * Grid * Grid; n++)
    {
      GridFirst[n] = -1;
      GridFlag[n] = 0;
    }

  for(n = 0; n < NumPart + Nghosts; n++)
    GridNext[n] = -1;


  for(n = 0; n < NumPart + Nghosts; n++)
    {
      if(n < NumPart)
	{
	  i = (P[n].Pos[0] - xmin[0]) * fac;
	  j = (P[n].Pos[1] - xmin[1]) * fac;
	  k = (P[n].Pos[2] - xmin[2]) * fac;
	}
      else
	{
	  i = (ImportBuf[n - NumPart].Pos[0] - xmin[0]) * fac;
	  j = (ImportBuf[n - NumPart].Pos[1] - xmin[1]) * fac;
	  k = (ImportBuf[n - NumPart].Pos[2] - xmin[2]) * fac;
	}

      if(i >= Grid)
	i = Grid - 1;
      if(j >= Grid)
	j = Grid - 1;
      if(k >= Grid)
	k = Grid - 1;

      if(GridFirst[(i * Grid + j) * Grid + k] >= 0)
	{
	  GridNext[n] = GridFirst[(i * Grid + j) * Grid + k];
	  GridFirst[(i * Grid + j) * Grid + k] = n;
	}
      else
	{
	  GridFirst[(i * Grid + j) * Grid + k] = n;
	}
    }
}



void fof_find_groups(void)
{
  int i, j, k;
  int p;

  if(ThisTask == 0)
    {
      printf("linking...");
      fflush(stdout);
    }

  for(i = Grid - 1; i >= 0; i--)
    for(j = Grid - 1; j >= 0; j--)
      for(k = Grid - 1; k >= 0; k--)
	{
	  if((p = GridFirst[(i * Grid + j) * Grid + k]) >= 0)
	    {
	      do
		{
		  fof_check_cell(p, i + 1, j, k);
		  fof_check_cell(p, i + 1, j + 1, k);
		  fof_check_cell(p, i + 1, j, k + 1);
		  fof_check_cell(p, i + 1, j + 1, k + 1);
		  fof_check_cell(p, i, j + 1, k);
		  fof_check_cell(p, i, j + 1, k + 1);
		  fof_check_cell(p, i, j, k + 1);

		  fof_check_cell(p, i, j, k);

		  fof_check_cell(p, i + 1, j, k - 1);
		  fof_check_cell(p, i + 1, j - 1, k);
		  fof_check_cell(p, i, j - 1, k + 1);
		  fof_check_cell(p, i - 1, j + 1, k + 1);
		  fof_check_cell(p, i - 1, j - 1, k + 1);
		  fof_check_cell(p, i + 1, j - 1, k + 1);
		}
	      while((p = GridNext[p]) >= 0);
	    }
	}
}



void fof_check_cell(int p, int i, int j, int k)
{
  double r2, s2, dx, dy, dz, xtmp;
  int pp, ss;
  int s, flag;

  if(i < 0)
    i += Grid;
  if(j < 0)
    j += Grid;
  if(k < 0)
    k += Grid;

  if(i >= Grid)
    i -= Grid;
  if(j >= Grid)
    j -= Grid;
  if(k >= Grid)
    k -= Grid;


  if((s = GridFirst[(i * Grid + j) * Grid + k]) >= 0)
    {
      flag = Head[s];
      if(GridFlag[(i * Grid + j) * Grid + k])
	{
	  if(Head[s] == Head[p])
	    return;
	}
    }
  else
    flag = -1;

  while(s >= 0)
    {
      if(Head[s] != flag)
	flag = -1;

      if(Head[p] != Head[s])	/* only if not yet linked */
	{
	  if(p < NumPart)
	    {
	      dx = P[p].Pos[0];
	      dy = P[p].Pos[1];
	      dz = P[p].Pos[2];
	    }
	  else
	    {
	      dx = ImportBuf[p - NumPart].Pos[0];
	      dy = ImportBuf[p - NumPart].Pos[1];
	      dz = ImportBuf[p - NumPart].Pos[2];
	    }
	  if(s < NumPart)
	    {
	      dx -= P[s].Pos[0];
	      dy -= P[s].Pos[1];
	      dz -= P[s].Pos[2];
	    }
	  else
	    {
	      dx -= ImportBuf[s - NumPart].Pos[0];
	      dy -= ImportBuf[s - NumPart].Pos[1];
	      dz -= ImportBuf[s - NumPart].Pos[2];
	    }

	  dx = FOF_PERIODIC(dx);
	  dy = FOF_PERIODIC(dy);
	  dz = FOF_PERIODIC(dz);

	  r2 = dx * dx + dy * dy + dz * dz;

	  s2 = LinkL * LinkL;

	  if(r2 <= s2)
	    {
	      if(Len[Head[p]] > Len[Head[s]])	/* p group is longer */
		{
		  Next[Tail[Head[p]]] = Head[s];

		  Tail[Head[p]] = Tail[Head[s]];

		  Len[Head[p]] += Len[Head[s]];

		  ss = Head[s];
		  do
		    {
		      Head[ss] = Head[p];
		    }
		  while((ss = Next[ss]) >= 0);

		  flag = -1;
		}
	      else
		{
		  Next[Tail[Head[s]]] = Head[p];

		  Tail[Head[s]] = Tail[Head[p]];

		  Len[Head[s]] += Len[Head[p]];

		  pp = Head[p];
		  do
		    {
		      Head[pp] = Head[s];
		    }
		  while((pp = Next[pp]) >= 0);

		  flag = -1;
		}

	      if(GridFlag[(i * Grid + j) * Grid + k])
		return;
	    }
	}

      s = GridNext[s];
    }

  if(flag >= 0)
    GridFlag[(i * Grid + j) * Grid + k] = 1;
}



void fof_find_minids(void)
{
  int n, pp;
  long long minid, id;

  int grfound = 0;

  for(n = 0; n < NumPart + Nghosts; n++)
    {
      if(Head[n] == n)
	{
	  if(n < NumPart)
	    minid = P[n].ID + (((long long) ThisTask) << 48);
	  else
	    minid = ImportBuf[n - NumPart].GrID;

	  pp = n;
	  do
	    {
	      if(pp < NumPart)
		id = P[pp].ID + (((long long) ThisTask) << 48);
	      else
		id = ImportBuf[pp - NumPart].GrID;

	      if(minid > id)
		minid = id;
	    }
	  while((pp = Next[pp]) >= 0);

	  GrID[n] = minid;

	  grfound++;
	}
    }
}


int fof_link_accross(void)
{
  int n, index, i, j, k, ix, iy, iz, links, tot_links;
  int level, sendTask, recvTask;
  double x, y, z, xtmp;
  long long headid;
  double scalefac;
  MPI_Status status;

  scalefac = 1.0 / All.BoxSize;

  links = 0;

  for(n = 0; n < Nghosts; n++)
    ImportBuf[n].GrID = GrID[Head[n + NumPart]];

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&ImportBuf[import_offset[recvTask]],
		     toGo[recvTask * NTask + ThisTask] * sizeof(struct export_fof),
		     MPI_BYTE, recvTask, TAG_FOF_LINK_A,
		     &ExportBuf[export_offset[recvTask]],
		     toGo[ThisTask * NTask + recvTask] * sizeof(struct export_fof),
		     MPI_BYTE, recvTask, TAG_FOF_LINK_A, MPI_COMM_WORLD, &status);
    }

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      for(j = 0; j < NTask; j++)
	Exportflag[j] = 0;

      for(ix = -1; ix <= 1; ix++)
	for(iy = -1; iy <= 1; iy++)
	  for(iz = -1; iz <= 1; iz++)
	    if(ix != 0 || iy != 0 || iz != 0)
	      {
		x = FOF_PERIODIC_WRAP(P[n].Pos[0] + ix * LinkL) * scalefac + 1.0;
		y = FOF_PERIODIC_WRAP(P[n].Pos[1] + iy * LinkL) * scalefac + 1.0;
		z = FOF_PERIODIC_WRAP(P[n].Pos[2] + iz * LinkL) * scalefac + 1.0;

		i = DOUBLE_to_DOMAINGRID(x);
		j = DOUBLE_to_DOMAINGRID(y);
		k = DOUBLE_to_DOMAINGRID(z);

		index = DomainPeanoMap[(i * DOMAINGRID + j) * DOMAINGRID + k];

		if(DomainTask[index] != ThisTask)
		  Exportflag[DomainTask[index]] = 1;
	      }

      for(j = 0; j < NTask; j++)
	if(Exportflag[j])
	  {
	    headid = ExportBuf[export_offset[j] + local_toGo[j]].GrID;

	    if(headid < GrID[Head[n]])
	      {
		GrID[Head[n]] = headid;
		links++;
	      }
	    else
	      headid = GrID[Head[n]];

	    ExportBuf[export_offset[j] + local_toGo[j]].GrID = headid;

	    local_toGo[j] += 1;
	  }
    }


  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&ExportBuf[export_offset[recvTask]],
		     toGo[ThisTask * NTask + recvTask] * sizeof(struct export_fof),
		     MPI_BYTE, recvTask, TAG_FOF_LINK_B,
		     &ImportBuf[import_offset[recvTask]],
		     toGo[recvTask * NTask + ThisTask] * sizeof(struct export_fof),
		     MPI_BYTE, recvTask, TAG_FOF_LINK_B, MPI_COMM_WORLD, &status);
    }

  for(n = 0; n < Nghosts; n++)
    {
      headid = ImportBuf[n].GrID;

      if(headid < GrID[Head[n + NumPart]])
	{
	  GrID[Head[n + NumPart]] = headid;
	  links++;
	}
    }

  MPI_Allreduce(&links, &tot_links, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  return tot_links;
}



void fof_exchange_id_lists(void)
{
  int n;
  long long headid;
  int task, nsend, nget;
  int level, sendTask, recvTask;
  int i, j, k, hash_key;
  double scalefac;
  double xx, yy, zz;

  MPI_Status status;

  scalefac = 1.0 / All.BoxSize;

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      headid = GrID[Head[n]];
      task = (headid >> 48);

      if(task != ThisTask)
	local_toGo[task]++;
    }

  MPI_Allgather(local_toGo, NTask, MPI_INT, toGo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(j = 0, nsend = nget = 0; j < NTask; j++)
    {
      nsend += toGo[ThisTask * NTask + j];
      nget += toGo[j * NTask + ThisTask];
    }

  printf("Task=%d needs to send %d particles for FoF-grouplist and gets %d particles.\n", ThisTask, nsend,
	 nget);
  fflush(stdout);

  Nlocal = NumPart - nsend + nget;

  local_ids = mymalloc(1 + Nlocal * sizeof(struct id_list));
  export_ids = mymalloc(1 + nsend * sizeof(struct id_list));

  if(!(export_ids) || !(local_ids))
    {
      printf("failed to allocate memory\n");
      endrun(152311);
    }

  import_offset[0] = NumPart - nsend;
  export_offset[0] = 0;

  for(n = 1; n < NTask; n++)
    {
      export_offset[n] = export_offset[n - 1] + toGo[ThisTask * NTask + (n - 1)];
      import_offset[n] = import_offset[n - 1] + toGo[(n - 1) * NTask + ThisTask];
    }

  for(n = 0; n < NTask; n++)
    local_toGo[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      headid = GrID[Head[n]];
      task = headid >> 48;

      /* compute hash key */
      xx = P[n].Pos[0] * scalefac + 1.0;
      yy = P[n].Pos[1] * scalefac + 1.0;
      zz = P[n].Pos[2] * scalefac + 1.0;

      i = DOUBLE_to_HASHBITS(xx);
      j = DOUBLE_to_HASHBITS(yy);
      k = DOUBLE_to_HASHBITS(zz);

      hash_key = peano_hilbert_key(i, j, k, HASHBITS);

      if(task != ThisTask)
	{
	  export_ids[export_offset[task] + local_toGo[task]].GrID = headid;
	  export_ids[export_offset[task] + local_toGo[task]].ID = P[n].ID + (((long long) hash_key) << 34);
	}
      else
	{
	  local_ids[local_toGo[task]].GrID = headid;
	  local_ids[local_toGo[task]].ID = P[n].ID + (((long long) hash_key) << 34);
	}

      local_toGo[task]++;
    }

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&export_ids[export_offset[recvTask]],
		     toGo[ThisTask * NTask + recvTask] * sizeof(struct id_list),
		     MPI_BYTE, recvTask, TAG_FOF_EXCH,
		     &local_ids[import_offset[recvTask]],
		     toGo[recvTask * NTask + ThisTask] * sizeof(struct id_list),
		     MPI_BYTE, recvTask, TAG_FOF_EXCH, MPI_COMM_WORLD, &status);
    }

  myfree(export_ids);
}



int fof_grid_compare(const void *a, const void *b)
{
  if(((struct id_list *) a)->GrID > (((struct id_list *) b)->GrID))
    return -1;

  if(((struct id_list *) a)->GrID < (((struct id_list *) b)->GrID))
    return +1;

  if(((struct id_list *) a)->ID < (((struct id_list *) b)->ID))
    return -1;

  if(((struct id_list *) a)->ID > (((struct id_list *) b)->ID))
    return +1;

  return 0;
}




void fof_compile_catalogue(void)
{
  int i, j, start, groups;

  qsort(local_ids, Nlocal, sizeof(struct id_list), fof_grid_compare);

  Ngroups = 0;
  Nids = 0;
  groups = 0;

  for(i = 0; i < GROUP_MIN_LEN; i++)
    CountBelowMinLen[i] = 0;

  for(i = 0, start = 0; i < Nlocal; i++)
    {
      if(local_ids[i].GrID != local_ids[start].GrID)
	{
	  for(j = start; j < i; j++)
	    local_ids[j].GrID = groups + (((long long) (i - start)) << 32);

	  if((i - start) >= GROUP_MIN_LEN)
	    {
	      Ngroups++;
	      Nids += (i - start);
	    }
	  else
	    CountBelowMinLen[(i - start)]++;

	  groups++;
	  start = i;
	}
    }

  /* finish last group */
  for(j = start; j < Nlocal; j++)
    local_ids[j].GrID = groups + (((long long) (Nlocal - start)) << 32);

  if((Nlocal - start) >= GROUP_MIN_LEN)
    {
      Ngroups++;
      Nids += (Nlocal - start);
    }
  else
    CountBelowMinLen[(Nlocal - start)]++;


  /* now sort the groups by size */
  qsort(local_ids, Nlocal, sizeof(struct id_list), fof_grid_compare);

  GroupLen = mymalloc((1 + Ngroups) * sizeof(int));
  GroupOffset = mymalloc((1 + Ngroups) * sizeof(int));
  GroupIDs = (long long *) local_ids;

  GroupLen[0] = 0; /* such that we can later determine the largest group even if some tasks don't have one (yet) */

  if(!(GroupLen) || !(GroupOffset))
    {
      printf("failed to allocate memory\n");
      endrun(1523711);
    }



  for(i = 0, j = 0, GroupOffset[0] = 0; i < Ngroups; i++)
    {
      GroupLen[i] = (local_ids[j].GrID >> 32);

      j += GroupLen[i];
      if(i > 0)
	GroupOffset[i] = GroupOffset[i - 1] + GroupLen[i - 1];
    }

  for(i = 0; i < Nids; i++)
    GroupIDs[i] = local_ids[i].ID;
}



void fof_save_groups(int num)
{
  int nprocgroup, masterTask, groupTask;

  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(241931);
    }

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	fof_save_local_catalogue(num);

      MPI_Barrier(MPI_COMM_WORLD);	/* wait inside the group */
    }
}



void fof_save_local_catalogue(int num)
{
  FILE *fd;
  char buf[500];
  int count;

  sprintf(buf, "%ssnapdir_%03d/%s_%03d.%d", All.OutputDir, num, "group_tab", num, ThisTask);

  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1183);
    }

  fwrite(&Ngroups, sizeof(int), 1, fd);
  fwrite(&Nids, sizeof(int), 1, fd);
  fwrite(&TotNgroups, sizeof(int), 1, fd);
  fwrite(&NTask, sizeof(int), 1, fd);
  fwrite(GroupLen, sizeof(int), Ngroups, fd);
  fwrite(GroupOffset, sizeof(int), Ngroups, fd);
  count = GROUP_MIN_LEN;
  fwrite(&count, sizeof(int), 1, fd);
  fwrite(CountBelowMinLen, sizeof(int), GROUP_MIN_LEN, fd);
  fclose(fd);

  sprintf(buf, "%ssnapdir_%03d/%s_%03d.%d", All.OutputDir, num, "group_ids", num, ThisTask);

  if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s`\n", buf);
      endrun(1184);
    }

  fwrite(&Ngroups, sizeof(int), 1, fd);
  fwrite(&Nids, sizeof(int), 1, fd);
  fwrite(&TotNgroups, sizeof(int), 1, fd);
  fwrite(&NTask, sizeof(int), 1, fd);
  fwrite(GroupIDs, sizeof(long long), Nids, fd);
  fclose(fd);
}

