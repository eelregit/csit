#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifndef SYSMALLOC

#define MAXBLOCKS 128

static size_t TotBytes;
static size_t FreeBytes;

static unsigned long Nblocks;
static double Highmark;

static void *Base;
static void **Table;
static size_t *BlockSize;

void mymalloc_init(size_t n)
{
  BlockSize = malloc(MAXBLOCKS * sizeof(size_t));
  Table = malloc(MAXBLOCKS * sizeof(void *));

  if(!(Base = malloc(n)))
    {
      printf("Failed to allocate memory for `Base' (%u bytes).\n", (int) n);
      endrun(2);
    }

  TotBytes = FreeBytes = n;
  Nblocks = 0;
  Highmark = 0;
}



void *mymalloc(size_t n)
{
  if((n % 8) > 0)
    n = (n / 8 + 1) * 8;

  if(n > FreeBytes)
    {
      printf("Task=%d: Not enough memory in mymalloc(n=%g MB).  FreeBytes=%g MB\n",
	     ThisTask, n / (1024.0 * 1024.0), FreeBytes / (1024.0 * 1024.0));
      endrun(812);
    }

  if(Nblocks >= MAXBLOCKS)
    {
      printf("Task=%d: No blocks left in mymalloc().\n", ThisTask);
      endrun(813);
    }

  Table[Nblocks] = Base + (TotBytes - FreeBytes);
  BlockSize[Nblocks] = n;

  FreeBytes -= n;
  Nblocks += 1;

  if((TotBytes - FreeBytes) / (1024.0 * 1024.0) > Highmark)
    {
      Highmark = (TotBytes - FreeBytes) / (1024.0 * 1024.0);
      if(Highmark > 500.0)
	{
	  printf("Task=%d:   new highmark=%g MB   fractional=%g\n",
		 ThisTask, Highmark, Highmark / (TotBytes / (1024.0 * 1024.0)));
	  fflush(stdout);
	}
    }

  return Table[Nblocks - 1];
}


void myfree(void *p)
{
  if(Nblocks == 0)
    endrun(76878);

  if(p != Table[Nblocks - 1])
    {
      printf("Task=%d: Wrong call of myfree() - not the last allocated block!\n", ThisTask);
      fflush(stdout);
      endrun(814);
    }

  Nblocks -= 1;
  FreeBytes += BlockSize[Nblocks];
}


void *myrealloc(void *p, size_t n)
{
  myfree(p);
  return mymalloc(n);
}

#endif
