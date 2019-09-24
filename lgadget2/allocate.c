#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"



/*! \file allocate.c
 *  \brief routines to allocate/free particle storage and comm-buffers
 */


/*! Allocates communication buffers.
 */
void allocate_commbuffers(void)
{
  size_t bytes;
  static int done_flag = 0;

  if(!(CommBuffer = mymalloc(bytes = All.BufferSize * 1024 * 1024)))
    {
      printf("LGadget Error: 1\n");
      printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(2);
    }

  All.BunchSizeForce =
    (All.BufferSize * 1024 * 1024) / (sizeof(struct gravdata_acctable) + sizeof(struct gravdata_index) +
				      2 * sizeof(struct gravdata_in));
  if(All.BunchSizeForce & 1)
    All.BunchSizeForce -= 1;	/* make sure that All.BunchSizeForce is even --> 8-byte alignment necessart for some 64bit processors */

  GravDataAccTable = (struct gravdata_acctable *) CommBuffer;
  GravDataIndexTable = (struct gravdata_index *) (GravDataAccTable + All.BunchSizeForce);
  GravDataIn = (struct gravdata_in *) (GravDataIndexTable + All.BunchSizeForce);
  GravDataGet = GravDataIn + All.BunchSizeForce;
  GravDataResult = GravDataGet;	/* this will overwrite the GravDataGet-Table */
  GravDataOut = GravDataIn;	/* this will overwrite the GravDataIn-Table */


  All.BunchSizeTwoPoint = (All.BufferSize * 1024 * 1024) / (2 * sizeof(struct twopointdata_in));
  TwoPointDataIn = (struct twopointdata_in *) CommBuffer;
  TwoPointDataGet = TwoPointDataIn + All.BunchSizeTwoPoint;


  if(ThisTask == 0 && done_flag == 0)
    {
      printf("Allocated %d MByte communication buffer per processor.\n", All.BufferSize);
      printf("Communication buffer has room for %d particles in gravity computation\n", All.BunchSizeForce);
      printf("Communication buffer has room for %d particles in two-point correlation function\n",
	     All.BunchSizeTwoPoint);
      fflush(stdout);
      done_flag = 1;
    }
}


void free_commbuffers(void)
{
  myfree(CommBuffer);
}


/* This routine allocates memory for 
 * particle storage, both the collisionless and the SPH particles.
 * The memory for the ordered binary tree of the timeline
 * is also allocated.
 */
void allocate_memory(void)
{
  size_t bytes, bytes_tot = 0;

  Exportflag = mymalloc(NTask * sizeof(char));
  DomainStartList = mymalloc(NTask * sizeof(int));
  DomainEndList = mymalloc(NTask * sizeof(int));

  if(!(P = mymalloc(bytes = All.MaxPart * sizeof(struct particle_data))))
    {
      printf("LGadget Error: 2");
      printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(ThisTask == 0)
    printf("\nAllocated %g MByte for particle storage.\n\n", bytes_tot / (1024.0 * 1024.0));
}

void reallocate_particle_memory_NumPart(void)
{
  size_t bytes;


  if(!(P = myrealloc(P, bytes = NumPart * sizeof(struct particle_data))))
    {
      printf("LGadget Error: 3");
      printf("Task=%d: failed to reallocate memory for `P' (%g MB).\n", ThisTask, bytes / (1024.0 * 1024.0));
      endrun(1);
    }

  if(ThisTask == 0)
    {
      printf("\nreduced particle storage to %g MByte\n", bytes / (1024.0 * 1024.0));
      fflush(stdout);
    }
}

void reallocate_particle_memory_MaxPart(void)
{
  size_t bytes;

  if(!(P = myrealloc(P, bytes = All.MaxPart * sizeof(struct particle_data))))
    {
      printf("Task=%d: failed to reallocate memory for `P' (%g MB).\n", ThisTask, bytes / (1024.0 * 1024.0));
      endrun(1);
    }

  if(ThisTask == 0)
    {
      printf("\nchanged particle storage to %g MByte\n", bytes / (1024.0 * 1024.0));
      fflush(stdout);
    }
}



/* This routine frees the memory for the particle storage
 * and for the communication buffers,
 * but we don't actually call it in the code. 
 * When the program terminats, the memory will be automatically
 * freed by the operating system.
 */
void free_memory(void)
{
  myfree(P);
  myfree(DomainEndList);
  myfree(DomainStartList);
  myfree(Exportflag);
}
