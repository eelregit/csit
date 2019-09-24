#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_heapsort.h>

#include "allvars.h"
#include "proto.h"
#include "domain.h"

#ifdef HPM
#include <libhpm.h>
#endif


/*! \file peano.c
 *  \brief constructs Peano-Hilbert ordering of particles
 */


static struct peano_hilbert_data
{
  peanokey key;
  int index;
}
 *mp;

static int *Id;


#define DOUBLE_to_PEANOGRID(y) ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - PEANO_BITS_PER_DIMENSION)))


void peano_hilbert_order(void)
{
  int i;
  int x, y, z;
  double t0, t1;
  double xx, yy, zz;
  double scalefac;

  scalefac = 1.0 / All.BoxSize;

  if(NumPart > 0)
    {
      if(ThisTask == 0)
	{
	  printf("Start Peano-Hilbert ordering.\n");
	  fflush(stdout);
	}

#ifdef HPM
      hpmStart(6, "Peano");
#endif
      mp = mymalloc(sizeof(struct peano_hilbert_data) * NumPart);
      Id = mymalloc(sizeof(int) * NumPart);

      t0 = second();
      for(i = 0; i < NumPart; i++)
	{
	  mp[i].index = i;

	  xx = P[i].Pos[0] * scalefac + 1.0;
	  yy = P[i].Pos[1] * scalefac + 1.0;
	  zz = P[i].Pos[2] * scalefac + 1.0;

	  x = DOUBLE_to_PEANOGRID(xx);
	  y = DOUBLE_to_PEANOGRID(yy);
	  z = DOUBLE_to_PEANOGRID(zz);

	  mp[i].key = peano_hilbert_key(x, y, z, PEANO_BITS_PER_DIMENSION);
	}
      t1 = second();

      if(ThisTask == 0)
	{
	  printf("keys generated (%g sec)\n", timediff(t0,t1));
	  fflush(stdout);
	}

      /* Note: In theory, Quicksort is slowest if the input array is already ordered!
       * This is basically the case for us because the order changes only
       * slowly. Hence depending on the detailed qsort() implementation, 
       * it may be faster to use the heapsort algorithm instead. This may be system 
       * andsize dependend.
       */

      t0 = second();
      qsort(mp, NumPart, sizeof(struct peano_hilbert_data), compare_key);
      /*
       * gsl_heapsort(mp, NumPart, sizeof(struct peano_hilbert_data), compare_key);
       */
      t1 = second();

      if(ThisTask == 0)
	{
	  printf("sort done. (%g sec)\n", timediff(t0,t1));
	  fflush(stdout);
	}

      t0 = second();

      for(i = 0; i < NumPart; i++)
	Id[mp[i].index] = i;
      reorder_particles();

      t1 = second();
      if(ThisTask == 0)
	{
	  printf("reordering done. (%g sec)\n", timediff(t0,t1));
	  fflush(stdout);
	}

      myfree(Id);
      myfree(mp);
#ifdef HPM
      hpmStop(6);
#endif
    }
}




int compare_key(const void *a, const void *b)
{
  if(((struct peano_hilbert_data *) a)->key < (((struct peano_hilbert_data *) b)->key))
    return -1;

  if(((struct peano_hilbert_data *) a)->key > (((struct peano_hilbert_data *) b)->key))
    return +1;

  return 0;
}


void reorder_particles(void)
{
  int i;
  struct particle_data Psave, Psource;
  int idsource, idsave, dest;

  for(i = 0; i < NumPart; i++)
    {
      if(Id[i] != i)
	{
	  Psource = P[i];
	  idsource = Id[i];

	  dest = Id[i];

	  do
	    {
	      Psave = P[dest];
	      idsave = Id[dest];

	      P[dest] = Psource;
	      Id[dest] = idsource;

	      if(dest == i)
		break;

	      Psource = Psave;
	      idsource = idsave;

	      dest = idsource;
	    }
	  while(1);
	}
    }
}






static char quadrants[24][2][2][2] = {
  /* rotx=0, roty=0-3 */
  {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
  {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
  {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
  {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
  /* rotx=1, roty=0-3 */
  {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
  {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
  {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
  {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
  /* rotx=2, roty=0-3 */
  {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
  {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
  {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
  {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
  /* rotx=3, roty=0-3 */
  {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
  {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
  {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
  {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
  /* rotx=4, roty=0-3 */
  {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
  {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
  {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
  {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
  /* rotx=5, roty=0-3 */
  {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
  {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
  {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
  {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
};


static char rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
  12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
};

static char rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
  11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
};

static char rotx_table[8] = { 3, 0, 0, 2, 2, 0, 0, 1 };
static char roty_table[8] = { 0, 1, 1, 2, 2, 3, 3, 0 };

static char sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };


peanokey peano_hilbert_key(int x, int y, int z, int bits)
{
  int i, bitx, bity, bitz, mask, quad, rotation;
  char sense, rotx, roty;
  peanokey key;

  mask = 1 << (bits - 1);
  key = 0;
  rotation = 0;
  sense = 1;

  for(i = 0; i < bits; i++, mask >>= 1)
    {
      bitx = (x & mask) ? 1 : 0;
      bity = (y & mask) ? 1 : 0;
      bitz = (z & mask) ? 1 : 0;

      quad = quadrants[rotation][bitx][bity][bitz];

      key <<= 3;
      key += (sense == 1) ? (quad) : (7 - quad);

      rotx = rotx_table[quad];
      roty = roty_table[quad];
      sense *= sense_table[quad];

      while(rotx > 0)
	{
	  rotation = rotxmap_table[rotation];
	  rotx--;
	}

      while(roty > 0)
	{
	  rotation = rotymap_table[rotation];
	  roty--;
	}
    }

  return key;
}
