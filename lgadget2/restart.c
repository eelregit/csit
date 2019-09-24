#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/file.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

/*! \file restart.c
 *  \brief reads or writes restart files
 */


static int fd;
static void in(int *x, int modus);
static void byten(void *x, size_t n, int modus);



/* This function reads or writes the restart files.  Each processor writes its
 * own restart file, with the I/O being done in parallel. To avoid congestion
 * of the disks you can tell the program to restrict the number of files that
 * are simultaneously written to NumFilesWrittenInParallel.  On T3E systems,
 * you can also enable preallocation on different partitions to ensure high
 * I/O bandwidth.
 *
 * If modus>0  the restart()-routine reads, 
 * if modus==0 it writes a restart file. 
 */
void restart(int modus)
{
  char buf[200], buf_bak[200], buf_mv[500];
  double save_PartAllocFactor;
  double t0, t1;
  int nprocgroup, masterTask, groupTask;
  struct global_data_all_processes all_task0;

  t0 = second();
  if(ThisTask == 0)
    {
      printf("Starting restart-routine\n");
      fflush(stdout);
    }

#ifdef LIGHTCONE
  int LC;
  char fname_p[500], fname_v[500], fname_i[500], fname_a[500];
#endif
  
  sprintf(buf, "%s%s.%d", All.OutputDir, All.RestartFile, ThisTask);
  sprintf(buf_bak, "%s%s.%d.bak", All.OutputDir, All.RestartFile, ThisTask);
  sprintf(buf_mv, "mv %s %s", buf, buf_bak);


  if((NTask < All.NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(2131);
    }

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    {
      nprocgroup++;
    }

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{

#ifdef LIGHTCONE
	  //we are writing a restart file so dump light cone buffers here 
	  //- this will update NPart_Written to the proper count before it is put in restart file below
	  if(!modus)
	    {
	      for(LC=0;LC<NLCS;LC++)
		{
		  if(LCNum[LC] == 0)
		    continue;
		  dump_lc_buf(LC);
		}
	    }
#endif
	  
	  if(modus)
	    {
	      if((fd = open(buf, O_RDONLY)) < 0)
		{
		  printf("Restart file '%s' not found.\n", buf);
		  endrun(7870);
		}

	      printf("task=%d: reading restart-file '%s'\n", ThisTask, buf);
	      fflush(stdout);
	    }
	  else
	    {
	      system(buf_mv);	/* move old restart files to .bak files */

	      if((fd = open(buf, O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR)) < 0)
		{
		  printf("Restart file '%s' cannot be opened.\n", buf);
		  endrun(7878);
		}

	      printf("task=%d: writing restart-file '%s'\n", ThisTask, buf);
	      fflush(stdout);
	    }


	  save_PartAllocFactor = All.PartAllocFactor;

	  /* common data  */
	  byten(&All, sizeof(struct global_data_all_processes), modus);

	  if(ThisTask == 0 && modus > 0)
	    all_task0 = All;

	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }

	  if(modus)		/* read */
	    {
	      All.PartAllocFactor = save_PartAllocFactor;
	      All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);

	      if(all_task0.Time != All.Time)
		{
		  printf("The restart file on task=%d is not consistent with the one on task=0\n", ThisTask);
		  fflush(stdout);
		  endrun(16);
		}

	      allocate_memory();
	    }

	  in(&NumPart, modus);

	  if(NumPart > All.MaxPart)
	    {
	      printf
		("it seems you have reduced(!) 'PartAllocFactor' below the value of %g needed to load the restart file.\n",
		 ((double) (NumPart / (All.TotNumPart / NTask))));
	      printf("fatal error\n");
	      endrun(22);
	    }

	  /* Particle data  */
	  byten(&P[0], NumPart * sizeof(struct particle_data), modus);

	  /* write state of random number generator */
	  byten(gsl_rng_state(random_generator), gsl_rng_size(random_generator), modus);

#ifdef LIGHTCONE
	  /* counts total # of parts written by each task in each light cone 
	     - will be used to truncate files if restart from restart files being written now
	     - has been updated from dump of LC memory buffers above
	  */
	  byten(&NPart_Written[0], NLCS*sizeof(long long), modus);
#endif
	  
	  close(fd);
	  
#ifdef LIGHTCONE
	  if(modus) //we are reading a restart file - need to truncate light cone buffer files to correct size given in restart files
	    {
	      for(LC=0;LC<NLCS;LC++)
		{
		  if(LCNum[LC] == 0)
		    continue;
		  
		  sprintf(fname_p, "%slightcone_buff/lightcone%03d/p/%s_Lightcone.p.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
		  if(access(fname_p,F_OK) == 0)
		    {
		      if(truncate(fname_p,3*sizeof(float)*NPart_Written[LC]) != 0)
			{
			  printf("LIGHT CONE CORRUPTION: could not truncate '%s' for light cone!\n",fname_p);
			  endrun(666);
			}
		    }
		  
		  sprintf(fname_v, "%slightcone_buff/lightcone%03d/v/%s_Lightcone.v.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
		  if(access(fname_v,F_OK) == 0)
		    {
		      if(truncate(fname_v,3*sizeof(float)*NPart_Written[LC]) != 0)
			{
			  printf("LIGHT CONE CORRUPTION: could not truncate '%s' for light cone!\n",fname_v);
			  endrun(666);
			}
		    }
		  
		  sprintf(fname_i, "%slightcone_buff/lightcone%03d/i/%s_Lightcone.i.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
		  if(access(fname_i,F_OK) == 0)
		    {
		      if(truncate(fname_i,sizeof(long long)*NPart_Written[LC]) != 0)
			{
			  printf("LIGHT CONE CORRUPTION: could not truncate '%s' for light cone!\n",fname_i);
			  endrun(666);
			}
		    }
		  
		  sprintf(fname_a, "%slightcone_buff/lightcone%03d/a/%s_Lightcone.a.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
		  if(access(fname_a,F_OK) == 0)
		    {
		      if(truncate(fname_a,sizeof(float)*NPart_Written[LC]) != 0)
			{
			  printf("LIGHT CONE CORRUPTION: could not truncate '%s' for light cone!\n",fname_a);
			  endrun(666);
			}
		    }
		}
	    }
#endif
	}
      else			/* wait inside the group */
	{
	  if(modus > 0 && groupTask == 0)	/* read */
	    {
	      MPI_Bcast(&all_task0, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }


  MPI_Barrier(MPI_COMM_WORLD);

  t1 = second();
  if(ThisTask == 0)
    {
      printf("finished restart-routine  (took %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }
}



/* reads/writes n bytes 
 */
void byten(void *x, size_t n, int modus)
{
  if(modus)
    {
      if(read(fd, x, n * sizeof(char)) != n * sizeof(char))
	{
	  printf("read error on task=%d (restart file appears to be truncated)\n", ThisTask);
	  fflush(stdout);
	  endrun(7);
	}
    }
  else
    {
      if(write(fd, x, n * sizeof(char)) != n * sizeof(char))
	{
	  printf("write error on task=%d upon writing restart file\n", ThisTask);
	  fflush(stdout);
	  endrun(8);
	}
    }
}


/* reads/writes one int 
 */
void in(int *x, int modus)
{
  if(modus)
    {
      if(read(fd, x, 1 * sizeof(int)) != sizeof(int))
	{
	  printf("read error on task=%d (restart file appears to be truncated)\n", ThisTask);
	  fflush(stdout);
	  endrun(7);
	}
    }
  else
    {
      if(write(fd, x, 1 * sizeof(int)) != sizeof(int))
	{
	  printf("write error on task=%d upon writing restart file\n", ThisTask);
	  fflush(stdout);
	  endrun(8);
	}
    }
}
