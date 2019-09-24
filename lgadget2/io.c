#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "allvars.h"
#include "proto.h"


/*! \file io.c
 *  \brief routines to write snapshot files
 */


/* This macro converts a float in the range [0,1[ to an integer in the
 * range [0,2^HASHBITS-1]. The reason for doing it in this elaborate way
 * is basically because on the Intel platform, a correct IEEE float to int
 * conversion is very expensive because it would require to (re)set the 
 * default rounding mode, which can only be done by flushing the FPU state.
 * As a result, most optimizing compilers do a somewhat fuzzy float-to-int
 * conversion, which can sometimes lead to a round-up, at other times to 
 * a round-down... We here need however reproducible rounding, hence we 
 * do things ourselves by exploiting float-number storage directly.
 */
#define DOUBLE_to_HASHBITS(y) ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - HASHBITS)))



/* This wrapper function select the desired output 
 * routine for snapshot files.
 */
void savepositions(int num)
{
  double t0, t1;
  int nprocgroup, groupTask, masterTask;
  char buf[1000];

  if(ThisTask == 0)
    {
      printf("\nproducing snapshot file... \n");
      fflush(stdout);
    }

  DomainDecomposition(1);	/* try to generate equal volume decomposition */

  t0 = second();

  if((NTask < All.NumFilesWrittenInParallel))
    {
      printf
	("Fatal error.\nNumber of processors must be a smaller or equal than `NumFilesWrittenInParallel'.\n");
      endrun(24131);
    }

  if(ThisTask == 0)
    {
      sprintf(buf, "%ssnapdir_%03d", All.OutputDir, num);
      mkdir(buf, 02755);
    }
  MPI_Barrier(MPI_COMM_WORLD);	/* wait to make sure that directory has been created */

  /*
     if(num == 0)
     {
     All.TwoPointFlag = num + 1;
     All.PowerSpecFlag = num + 1;
     return;
     }
   */

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;


  masterTask = (ThisTask / nprocgroup) * nprocgroup;


  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	save_snapshot(num);

      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }

  All.TwoPointFlag = num + 1;
  All.PowerSpecFlag = num + 1;

  t1 = second();

  All.CPU_Snapshot += timediff(t0, t1);

  if(ThisTask == 0)
    {
      printf("done with writing snapshot.  (I/O took %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }



#ifdef FOF
  fof_fof(num);
#endif
}





/* This function writes a snapshot of the particle ditribution to one
 * or several files using Gadget's default file format.  If
 * NumFilesPerSnapshot>1, the snapshot is distributed into several
 * files, which are written simultaneously. Each file contains data
 * from a group of processors of size NTasl/NumFilesPerSnapshot.
 *
 * Each snapshot file contains a header first, then particle
 * positions, velocities and ID's.  Then particle masses are written
 * for those particle types with zero entry in MassTable.  After that,
 * first the internal energies u, and then the density is written for
 * the SPH particles.  Finally, if cooling is enabled, the mean
 * molecular weight is written for the gas particles.
 */
void save_snapshot(int num)
{
  size_t bytes;
  float *block;
  int *blockid;
  long long *blocklongid;
  int blockmaxlen, maxidlen, maxlongidlen;
  int4byte dummy;
  FILE *fd;
  char buf[300];
  int i, j, k, n, pc;
  double xx, yy, zz;
  double sqrta3inv;
  double scalefac;
  int prev_cell, hash_cell, first_hash_cell, last_hash_cell;



  sqrta3inv = sqrt(1 / (All.Time * All.Time * All.Time));

  sprintf(buf, "%ssnapdir_%03d/%s_%03d.%d", All.OutputDir, num, All.SnapshotFileBase, num, ThisTask);

  if(!(fd = fopen(buf, "w")))
    {
      printf("Error. Can't write in file '%s'\n", buf);
      endrun(10);
    }

  for(i = 0; i < 6; i++)
    {
      header.npart[i] = 0;
      header.npartTotal[i] = 0;
      header.mass[i] = 0;
    }

  header.npart[1] = NumPart;
  header.npartTotal[1] = All.TotNumPart;
  header.npartTotal[2] = (All.TotNumPart >> 32);
  header.mass[1] = All.PartMass;
  
  header.time = All.Time;
  header.redshift = 1.0 / All.Time - 1;

  header.flag_sfr = 0;
  header.flag_feedback = 0;
  header.flag_cooling = 0;
  header.flag_stellarage = 0;
  header.flag_metals = 0;

  header.num_files = NTask;

  header.BoxSize = All.BoxSize;
  header.Omega0 = All.Omega0;
  header.OmegaLambda = All.OmegaLambda;
  header.HubbleParam = All.HubbleParam;

  header.flag_stellarage = 0;
  header.flag_metals = 0;
  header.hashtabsize = (1 << (3 * HASHBITS));

  dummy = sizeof(header);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  my_fwrite(&header, sizeof(header), 1, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);


  if(!(block = mymalloc(bytes = 0.1 * All.MaxPart * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for `block' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(24);
    }

  blockmaxlen = bytes / (3 * sizeof(float));

  blockid = (int *) block;
  blocklongid = (long long *) block;
  maxidlen = bytes / (sizeof(int));
  maxlongidlen = bytes / (sizeof(long long));

  /* write coordinates */
  dummy = sizeof(float) * 3 * NumPart;
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = P[i].Pos[k];

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);



  /* write velocities */
  dummy = sizeof(float) * 3 * NumPart;
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      for(k = 0; k < 3; k++)
	block[3 * pc + k] = P[i].Vel[k] * sqrta3inv;

      pc++;

      if(pc == blockmaxlen)
	{
	  my_fwrite(block, sizeof(float), 3 * pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(block, sizeof(float), 3 * pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);


  /* write particle ID */
  dummy = sizeof(long long) * NumPart;
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  for(i = 0, pc = 0; i < NumPart; i++)
    {
      blocklongid[pc] = P[i].ID;

      pc++;

      if(pc == maxlongidlen)
	{
	  my_fwrite(blocklongid, sizeof(long long), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(blocklongid, sizeof(long long), pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);


  /* write hash-table */

  first_hash_cell = DomainMyStart * (1 << (3 * (HASHBITS - DOMAINLEVELS)));
  last_hash_cell = (DomainMyLast + 1) * (1 << (3 * (HASHBITS - DOMAINLEVELS))) - 1;


  /* let's do a check first */
  scalefac = 1.0 / All.BoxSize;

  for(n = 0, pc = 0, prev_cell = -1; n < NumPart; n++)
    {
      xx = P[n].Pos[0] * scalefac + 1.0;
      yy = P[n].Pos[1] * scalefac + 1.0;
      zz = P[n].Pos[2] * scalefac + 1.0;

      i = DOUBLE_to_HASHBITS(xx);
      j = DOUBLE_to_HASHBITS(yy);
      k = DOUBLE_to_HASHBITS(zz);

      hash_cell = peano_hilbert_key(i, j, k, HASHBITS);

      if(hash_cell < first_hash_cell)
	{
	  printf("Task=%d hash_cell=%d first_hash_cell=%d \n", ThisTask, hash_cell, first_hash_cell);

	  printf("Task=%d i,j,k= %d %d %d   %g %g %g  scalefac=%g\n",
		 ThisTask, i, j, k, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], scalefac);

	  fflush(stdout);

	  endrun(1313);
	}

      if(hash_cell > last_hash_cell)
	endrun(13134);

      if(hash_cell < prev_cell)
	endrun(13135);

      prev_cell = hash_cell;
    }

  dummy = sizeof(int) * 2;
  my_fwrite(&dummy, sizeof(dummy), 1, fd);
  my_fwrite(&first_hash_cell, sizeof(int), 1, fd);
  my_fwrite(&last_hash_cell, sizeof(int), 1, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  dummy = sizeof(int) * (last_hash_cell - first_hash_cell + 1);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);



  prev_cell = first_hash_cell - 1;

  for(n = 0, pc = 0; n < NumPart; n++)
    {
      xx = P[n].Pos[0] * scalefac + 1.0;
      yy = P[n].Pos[1] * scalefac + 1.0;
      zz = P[n].Pos[2] * scalefac + 1.0;

      i = DOUBLE_to_HASHBITS(xx);
      j = DOUBLE_to_HASHBITS(yy);
      k = DOUBLE_to_HASHBITS(zz);

      hash_cell = peano_hilbert_key(i, j, k, HASHBITS);

      if(hash_cell != prev_cell)
	{
	  for(i = prev_cell + 1; i <= hash_cell; i++)
	    {
	      blockid[pc] = n;
	      pc++;

	      if(pc == maxidlen)
		{
		  my_fwrite(blockid, sizeof(int), pc, fd);
		  pc = 0;
		}
	    }
	  prev_cell = hash_cell;
	}
    }

  for(i = prev_cell + 1; i <= last_hash_cell; i++)
    {
      blockid[pc] = NumPart;
      pc++;

      if(pc == maxidlen)
	{
	  my_fwrite(blockid, sizeof(int), pc, fd);
	  pc = 0;
	}
    }
  if(pc > 0)
    my_fwrite(blockid, sizeof(int), pc, fd);
  my_fwrite(&dummy, sizeof(dummy), 1, fd);

  myfree(block);

  fclose(fd);
}


/* This catches I/O errors occuring for my_fwrite(). In this case we better stop.
 */
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nwritten;

  if((nwritten = fwrite(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fwrite) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      endrun(777);
    }
  return nwritten;
}


/* This catches I/O errors occuring for fread(). In this case we better stop.
 */
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  size_t nread;

  if((nread = fread(ptr, size, nmemb, stream)) != nmemb)
    {
      printf("I/O error (fread) on task=%d has occured.\n", ThisTask);
      fflush(stdout);
      endrun(778);
    }
  return nread;
}
