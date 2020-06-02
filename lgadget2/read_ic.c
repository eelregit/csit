#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file read_ic.c
 *  \brief code to read in initial conditions/snapshot files
 */



void read_ic(char *fname)
{
  char buf[100];
  int num_files, rest_files;
  int masterTask, lastTask, filenr;
  int ngroups, gr, groupMaster;
  double t0, t1;

#ifdef RESCALEVINI
 if(ThisTask == 0)
   {
     fprintf(stdout, "\nRescaling v_ini !\n\n");
     fflush(stdout);
   }
#endif

  num_files = find_files(fname);
  All.TotNumPart = header.npartTotal[1] + (((long long) header.npartTotal[2]) << 32);
  All.PartMass = header.mass[1];
  All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */
  if(ThisTask == 0)
    printf("Allocating memory for %d particles\n", All.MaxPart);
  allocate_memory();

  if(ThisTask == 0)
    {
      printf("\nPartMass= %g\n", All.PartMass);
      printf("TotNumPart= %d%09d\n\n",
	     (int) (All.TotNumPart / 1000000000), (int) (All.TotNumPart % 1000000000));
      fflush(stdout);
    }

  if(RestartFlag == 2)
    All.Time = All.TimeBegin = header.time;


  if(NTask < All.NumFilesWrittenInParallel)
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(24132131);
    }

  t0 = second();

  NumPart = 0;

  rest_files = num_files;

  while(rest_files >= NTask)
    {
      sprintf(buf, "%s.%d", fname, ThisTask + (rest_files - NTask));

      ngroups = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	ngroups++;
      groupMaster = (ThisTask / ngroups) * ngroups;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if(ThisTask == (groupMaster + gr))	/* ok, it's this processor's turn */
	    read_file(buf, ThisTask, ThisTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      rest_files -= NTask;
    }


  if(rest_files > 0)
    {
      distribute_file(rest_files, 0, 0, NTask - 1, &filenr, &masterTask, &lastTask);

      if(num_files > 1)
	sprintf(buf, "%s.%d", fname, filenr);
      else
	sprintf(buf, "%s", fname);

      ngroups = rest_files / All.NumFilesWrittenInParallel;
      if((rest_files % All.NumFilesWrittenInParallel))
	ngroups++;

      for(gr = 0; gr < ngroups; gr++)
	{
	  if((filenr / All.NumFilesWrittenInParallel) == gr)	/* ok, it's this processor's turn */
	    read_file(buf, masterTask, lastTask);
	  MPI_Barrier(MPI_COMM_WORLD);
	}
    }

  t1 = second();

  if(ThisTask == 0)
    {
      printf("\nReading of IC-files finished (took %g sec)\n\n", timediff(t0, t1));
      fflush(stdout);
    }
}

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last)
{
  int ntask, filesleft, filesright, tasksleft, tasksright;

  if(nfiles > 1)
    {
      ntask = lasttask - firsttask + 1;

      filesleft = (((double) (ntask / 2)) / ntask) * nfiles;
      if(filesleft <= 0)
	filesleft = 1;
      if(filesleft >= nfiles)
	filesleft = nfiles - 1;

      filesright = nfiles - filesleft;

      tasksleft = ntask / 2;
      tasksright = ntask - tasksleft;

      distribute_file(filesleft, firstfile, firsttask, firsttask + tasksleft - 1, filenr, master, last);
      distribute_file(filesright, firstfile + filesleft, firsttask + tasksleft, lasttask, filenr, master,
		      last);
    }
  else
    {
      if(ThisTask >= firsttask && ThisTask <= lasttask)
	{
	  *filenr = firstfile;
	  *master = firsttask;
	  *last = lasttask;
	}
    }
}


int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int dummy;

  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd = fopen(buf, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);

      return header.num_files;
    }

  if((fd = fopen(buf1, "r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header, sizeof(header), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      fclose(fd);
      header.num_files = 1;

      return header.num_files;
    }


  endrun(1);
  return 0;
}




void read_file(char *fname, int readTask, int lastTask)
{
  double *block;
  int *blockid;
  long long *blocklongid;
  int blockmaxlen, maxidlen, type_of_id;
  int n_in_file, n_for_this_task, ntask, n, k, pc, offset = 0, task;
  int4byte blksize;
  size_t bytes;
  MPI_Status status;
  FILE *fd;

#define SKIP  {my_fread(&blksize,sizeof(int4byte),1,fd);}

  if(!(block = mymalloc(bytes = 0.1 * All.MaxPart * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for `block' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(24);
    }

  blockmaxlen = bytes / (3 * sizeof(double));
  blockid = (int *) (block);
  blocklongid = (long long *) (block);


  if(ThisTask == readTask)
    {
      if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s' for reading initial conditions.\n", fname);
	  endrun(123);
	}

      SKIP;
      if(blksize != 256)
	{
	  printf("incorrect header format (2)\n");
	  fflush(stdout);
	  endrun(889);
	}
      my_fread(&header, sizeof(header), 1, fd);
      SKIP;

      n_in_file = header.npart[1];

      printf("reading file `%s' on task=%d (contains %d particles.)\n"
	     "distributing this file to tasks %d-%d\n\n", fname, ThisTask, n_in_file, readTask, lastTask);
      fflush(stdout);

      ntask = lastTask - readTask + 1;

      /* read coordinates */
      SKIP;
      for(task = readTask; task <= lastTask; task++)
	{
	  n_for_this_task = n_in_file / ntask;
	  if((task - readTask) < (n_in_file % ntask))
	    n_for_this_task++;

	  if(task != ThisTask)
	    MPI_Send(&n_for_this_task, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
	  else
	    {
	      if(NumPart + n_for_this_task > All.MaxPart)
		{
		  printf("too many particles\n");
		  endrun(1313);
		}
	    }

	  if(ThisTask == task)
	    offset = 0;

	  do
	    {
	      pc = n_for_this_task;

	      if(pc > blockmaxlen)
		pc = blockmaxlen;

	      my_fread(block, sizeof(double), 3 * pc, fd);

	      if(task != ThisTask)
		MPI_Ssend(block, sizeof(double) * 3 * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD);
	      else
		{
		  for(n = 0; n < pc; n++)
		    for(k = 0; k < 3; k++)
		      P[NumPart + offset + n].Pos[k] = block[3 * n + k];

		  offset += pc;
		}
	      n_for_this_task -= pc;
	    }
	  while(n_for_this_task > 0);
	}
      SKIP;



      /* read velocities */
      SKIP;
      for(task = readTask; task <= lastTask; task++)
	{
	  n_for_this_task = n_in_file / ntask;
	  if((task - readTask) < (n_in_file % ntask))
	    n_for_this_task++;

	  if(task != ThisTask)
	    MPI_Send(&n_for_this_task, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);

	  if(ThisTask == task)
	    offset = 0;

	  do
	    {
	      pc = n_for_this_task;

	      if(pc > blockmaxlen)
		pc = blockmaxlen;

	      my_fread(block, sizeof(double), 3 * pc, fd);

	      if(task != ThisTask)
		MPI_Ssend(block, sizeof(double) * 3 * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD);
	      else
		{
		  for(n = 0; n < pc; n++)
		    for(k = 0; k < 3; k++)
#ifdef RESCALEVINI
        /* scaling v to use same IC's for different cosmologies */
                      P[NumPart + offset + n].Vel[k] = (block[3 * n + k]) * All.VelIniScale;
#else
                      P[NumPart + offset + n].Vel[k] = block[3 * n + k];
#endif

		  offset += pc;
		}
	      n_for_this_task -= pc;
	    }
	  while(n_for_this_task > 0);
	}
      SKIP;



      /* read IDs */
      SKIP;
      if(blksize == sizeof(int) * header.npart[1])
	type_of_id = sizeof(int);
      else
	type_of_id = sizeof(long long);

      maxidlen = bytes / type_of_id;

      for(task = readTask; task <= lastTask; task++)
	{
	  n_for_this_task = n_in_file / ntask;
	  if((task - readTask) < (n_in_file % ntask))
	    n_for_this_task++;

	  if(task != ThisTask)
	    {
	      MPI_Send(&n_for_this_task, 1, MPI_INT, task, TAG_N, MPI_COMM_WORLD);
	      MPI_Send(&type_of_id, 1, MPI_INT, task, TAG_FLAG, MPI_COMM_WORLD);
	    }

	  if(ThisTask == task)
	    offset = 0;

	  do
	    {
	      pc = n_for_this_task;

	      if(pc > maxidlen)
		pc = maxidlen;

	      my_fread(blockid, type_of_id, pc, fd);

	      if(task != ThisTask)
		MPI_Ssend(blockid, type_of_id * pc, MPI_BYTE, task, TAG_PDATA, MPI_COMM_WORLD);
	      else
		{
		  if(type_of_id == sizeof(int))
		    {
		      for(n = 0; n < pc; n++)
			P[NumPart + offset + n].ID = blockid[n];
		    }
		  else
		    {
		      for(n = 0; n < pc; n++)
			P[NumPart + offset + n].ID = blocklongid[n];
		    }
		  offset += pc;
		}
	      n_for_this_task -= pc;
	    }
	  while(n_for_this_task > 0);
	}
      SKIP;

      fclose(fd);
    }
  else
    {
      /* receive coordinates */
      MPI_Recv(&n_for_this_task, 1, MPI_INT, readTask, TAG_N, MPI_COMM_WORLD, &status);

      if(NumPart + n_for_this_task > All.MaxPart)
	{
	  printf("too many particles\n");
	  endrun(1313);
	}

      offset = 0;
      do
	{
	  pc = n_for_this_task;

	  if(pc > blockmaxlen)
	    pc = blockmaxlen;

	  MPI_Recv(block, sizeof(double) * 3 * pc, MPI_BYTE, readTask, TAG_PDATA, MPI_COMM_WORLD, &status);

	  for(n = 0; n < pc; n++)
	    for(k = 0; k < 3; k++)
	      P[NumPart + offset + n].Pos[k] = block[3 * n + k];

	  offset += pc;

	  n_for_this_task -= pc;
	}
      while(n_for_this_task > 0);




      /* receive velocities */
      MPI_Recv(&n_for_this_task, 1, MPI_INT, readTask, TAG_N, MPI_COMM_WORLD, &status);

      offset = 0;
      do
	{
	  pc = n_for_this_task;

	  if(pc > blockmaxlen)
	    pc = blockmaxlen;

	  MPI_Recv(block, sizeof(double) * 3 * pc, MPI_BYTE, readTask, TAG_PDATA, MPI_COMM_WORLD, &status);

	  for(n = 0; n < pc; n++)
	    for(k = 0; k < 3; k++)
	      P[NumPart + offset + n].Vel[k] = block[3 * n + k];

	  offset += pc;

	  n_for_this_task -= pc;
	}
      while(n_for_this_task > 0);



      /* receive IDs */
      MPI_Recv(&n_for_this_task, 1, MPI_INT, readTask, TAG_N, MPI_COMM_WORLD, &status);
      MPI_Recv(&type_of_id, 1, MPI_INT, readTask, TAG_FLAG, MPI_COMM_WORLD, &status);
      maxidlen = bytes / type_of_id;

      offset = 0;
      do
	{
	  pc = n_for_this_task;

	  if(pc > maxidlen)
	    pc = maxidlen;

	  MPI_Recv(blockid, type_of_id * pc, MPI_BYTE, readTask, TAG_PDATA, MPI_COMM_WORLD, &status);

	  if(type_of_id == sizeof(int))
	    {
	      for(n = 0; n < pc; n++)
		P[NumPart + offset + n].ID = blockid[n];
	    }
	  else
	    {
	      for(n = 0; n < pc; n++)
		P[NumPart + offset + n].ID = blocklongid[n];
	    }

	  offset += pc;

	  n_for_this_task -= pc;
	}
      while(n_for_this_task > 0);

    }

  NumPart += offset;

  myfree(block);
}
