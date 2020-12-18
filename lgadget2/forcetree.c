#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <gsl/gsl_sf_gamma.h>

#include "allvars.h"
#include "proto.h"


/*! \file forcetree.c
 *  \brief gravitational tree code for short-range forces
 */


/*
 *  This file contains the computation of the gravitational force by means of a tree. The type of tree
 *  implemented is a geometrical oct-tree, starting from a cube encompassing all particles. This cube is
 *  automatically found.  For particle types with different softenings, separate trees are always
 *  constructed. If 'COMBINETREES' is defined during compilation, different particle types with the same
 *  softening are put into the same tree.  The multipole moments of branches of the tree are dynamically
 *  updated during the tree walk, if needed.
 */


static int last;		/* auxialiary variable used to set-up non-recursive walk */

#define NTAB 10000
//static double shortrange_table[NTAB];
static double dI3_table[NTAB];
static double dI5_table[NTAB];
static double dI7_table[NTAB];
static double dI9_table[NTAB];
static double dI11_table[NTAB];
static double I7_table[NTAB];
static double I9_table[NTAB];
static double I11_table[NTAB];


#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))


#define DOUBLE_to_DOMAINGRID(y)  ((int)(((*((long long *) &y)) & 0xFFFFFFFFFFFFFllu) >> (52 - DOMAINLEVELS)))



/* index convention for accessing tree nodes:
   the indices   0...NumPart  reference single particles
   the indices   All.MaxPart.... All.MaxPart+nodes-1   reference tree nodes
    `Nodes_base' points to the first tree node
   `nodes' is shifted, such that nodes[All.MaxPart] gives the first tree node
*/



int force_treebuild(void)
{
  int i, j, subnode = 0, parent = -1, numnodes;
  int ii, jj, kk, index;
  int nfree, th, nn;
  double xx, yy, zz;
  double lenhalf;
  struct NODE *nfreep;
  double scalefac;

  scalefac = 1.0 / All.BoxSize;

  Nodes = Nodes_base - All.MaxPart;

  /* select first node */
  nfree = All.MaxPart;		/* index */
  nfreep = &Nodes[nfree];

  /* create an empty  root node  */
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];
  nfreep->len = DomainLen;

  for(i = 0; i < 8; i++)
    nfreep->u.suns[i] = -1;

  numnodes = 1;
  nfreep++;
  nfree++;


  /* create a set of empty nodes corresponding to the top-level domain grid  */
  force_create_empty_nodes(All.MaxPart, 0, 0, 0, 0, &numnodes, &nfree);

  nfreep = &Nodes[nfree];

  for(i = 0; i < NumPart; i++)	/* insert all particles */
    {
      xx = P[i].Pos[0] * scalefac + 1.0;
      yy = P[i].Pos[1] * scalefac + 1.0;
      zz = P[i].Pos[2] * scalefac + 1.0;

      ii = DOUBLE_to_DOMAINGRID(xx);
      jj = DOUBLE_to_DOMAINGRID(yy);
      kk = DOUBLE_to_DOMAINGRID(zz);

      index = DomainPeanoMap[(ii * DOMAINGRID + jj) * DOMAINGRID + kk];

      if(index < DomainMyStart || index > DomainMyLast)
	{
	  printf("Prob: Task=%d id=%d ii=%d jj=%d kk=%d   (%g|%g|%g) scalfac=%g %d %d %d\n",
		 ThisTask, (int) P[i].ID, ii, jj, kk, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], scalefac,
		 DomainMyStart, DomainMyLast, last);
	  fflush(stdout);
	  endrun(123131);
	}

      th = DomainNodeIndex[index];

      while(1)
	{
	  if(th >= All.MaxPart)	/* we are dealing with an internal node */
	    {
	      subnode = 0;
	      if(P[i].Pos[0] > Nodes[th].center[0])
		subnode += 1;
	      if(P[i].Pos[1] > Nodes[th].center[1])
		subnode += 2;
	      if(P[i].Pos[2] > Nodes[th].center[2])
		subnode += 4;

	      nn = Nodes[th].u.suns[subnode];

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = th;	/* note: subnode can still be used in the next step of the walk */
		  th = nn;
		}
	      else
		{
		  /* here we have found an empty slot where we can 
		   * attach the new particle as a leaf 
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* we try to insert into a leaf with a single particle
	       * need to generate a new internal node at this point 
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;

	      subnode = 0;
	      if(P[th].Pos[0] > nfreep->center[0])
		subnode += 1;
	      if(P[th].Pos[1] > nfreep->center[1])
		subnode += 2;
	      if(P[th].Pos[2] > nfreep->center[2])
		subnode += 4;
#ifndef NOTREERND
	      if(nfreep->len < 1.0e-3 * All.ComovSoftening)
		{
		  /* seems like we're dealing with particles   
		   * at identical locations. randomize 
		   * subnode index (well below gravitational softening scale). 
		   */
		  subnode = (int) (8.0 * gsl_rng_uniform(random_generator));
		  if(subnode >= 8)
		    subnode = 7;
		}
#endif
	      nfreep->u.suns[subnode] = th;

	      th = nfree;	/* resume trying to insert the new particle at 
				   the newly created internal node */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  /*
		     printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		     printf("for particle %d\n", i);
		     endrun(1);
		   */
		  return MaxNodes;	/* this will make the code try it once more with higher allocation factor */
		}
	    }
	}
    }


  force_insert_pseudo_particles();




  /* now compute the multipole moments recursively */
  last = -1;
  force_update_node_recursive(All.MaxPart, -1, -1);

  if(last >= All.MaxPart)
    {
      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
	DomainNextnode[last - (All.MaxPart + MaxNodes)] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  /* import top-level moments and update top-level tree */

  force_exchange_pseudodata();

  nfree = All.MaxPart + 1;
  force_treeupdate_toplevel(All.MaxPart, 0, &nfree);

  force_flag_localnodes();

  return numnodes;
}


/*! We create a set of empty tree nodes with a set of finest nodes that correspond to the domain grid. This is
 *  done to ensure that always internal nodes are created that correspond to the domain grid, which might not
 *  be the case for sparse particle population at the fringes of a particle set
*/
void force_create_empty_nodes(int no, int level, int nx, int ny, int nz, int *nodecount, int *nextfree)
{
  int count, i, j, k, n, index;

  if(level < DOMAINLEVELS)
    {
      level += 1;
      count = 0;

      for(k = -1; k <= 1; k += 2)
	for(j = -1; j <= 1; j += 2)
	  for(i = -1; i <= 1; i += 2)
	    {
	      Nodes[no].u.suns[count] = *nextfree;

	      Nodes[*nextfree].len = 0.5 * Nodes[no].len;
	      Nodes[*nextfree].center[0] = Nodes[no].center[0] + i * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[1] = Nodes[no].center[1] + j * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[2] = Nodes[no].center[2] + k * 0.25 * Nodes[no].len;
	      for(n = 0; n < 8; n++)
		Nodes[*nextfree].u.suns[n] = -1;

	      if(level == DOMAINLEVELS)
		{
		  index =
		    DomainPeanoMap[((nx * 2 + (i + 1) / 2) * DOMAINGRID + ny * 2 + (j + 1) / 2) * DOMAINGRID +
				   nz * 2 + (k + 1) / 2];

		  DomainNodeIndex[index] = *nextfree;
		}

	      *nextfree = *nextfree + 1;
	      *nodecount = *nodecount + 1;
	      count += 1;

	      if((*nodecount) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("in create empty nodes\n");
		  endrun(11);
		}

	      force_create_empty_nodes(*nextfree - 1, level,
				       nx * 2 + (i + 1) / 2, ny * 2 + (j + 1) / 2, nz * 2 + (k + 1) / 2,
				       nodecount, nextfree);
	    }
    }
}






/* now we insert some dummy pseudo-particles which represent the other domain grid cells */

void force_insert_pseudo_particles(void)
{
  int i, ii, jj, kk, index, subnode, nn, th;

  for(ii = 0; ii < DOMAINGRID; ii++)
    for(jj = 0; jj < DOMAINGRID; jj++)
      for(kk = 0; kk < DOMAINGRID; kk++)
	{
	  index = DomainPeanoMap[(ii * DOMAINGRID + jj) * DOMAINGRID + kk];

	  DomainMoment[index].mass = 0;
	  DomainMoment[index].s[0] = (ii + 0.5) / DomainFac;
	  DomainMoment[index].s[1] = (jj + 0.5) / DomainFac;
	  DomainMoment[index].s[2] = (kk + 0.5) / DomainFac;
	}

  for(ii = 0; ii < DOMAINGRID; ii++)
    for(jj = 0; jj < DOMAINGRID; jj++)
      for(kk = 0; kk < DOMAINGRID; kk++)
	{
	  i = DomainPeanoMap[(ii * DOMAINGRID + jj) * DOMAINGRID + kk];

	  if(i < DomainMyStart || i > DomainMyLast)
	    {
	      th = All.MaxPart;	/* select index of first node in tree */

	      while(1)
		{
		  if(th >= All.MaxPart)	/* we are dealing with an internal node */
		    {
		      if(th >= All.MaxPart + MaxNodes)
			{
			  endrun(888);	/* pseudo particle: this can't be */
			}

		      subnode = 0;
		      if(DomainMoment[i].s[0] > Nodes[th].center[0])
			subnode += 1;
		      if(DomainMoment[i].s[1] > Nodes[th].center[1])
			subnode += 2;
		      if(DomainMoment[i].s[2] > Nodes[th].center[2])
			subnode += 4;

		      nn = Nodes[th].u.suns[subnode];

		      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
			{
			  th = nn;
			}
		      else
			{
			  /* here we have found an empty slot where we can attach the pseudo particle as a
			   * leaf
			   */
			  Nodes[th].u.suns[subnode] = All.MaxPart + MaxNodes + i;

			  break;	/* done for this pseudo particle */
			}
		    }
		  else
		    {
		      endrun(889);	/* this can't be */
		    }
		}
	    }
	}
}






/* this routine computes the multipole moments for a given internal node and all its subnodes using a
 * recursive computation.  Note that the moments of the daughter nodes are already stored in single
 * precision. For very large particle numbers, loss of precision may results for certain particle
 * distributions
 */
void force_update_node_recursive(int no, int sib, int father)
{
  int j, jj, toplevel, p, pp = 0, nextsib, mass, suns[8];
  double s[3];



  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      if(Nodes[no].len > 0.9 / DomainFac)
	toplevel = 1;
      else
	toplevel = 0;

      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		DomainNextnode[last - (All.MaxPart + MaxNodes)] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no);

	      if(p >= All.MaxPart)	/* an internal node or pseudo particle */
		{
		  if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		    {
		      /* nothing to be done here */
		    }
		  else
		    {
		      mass += Nodes[p].u.d.mass;	/* we assume a fixed particle mass */
		      s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		      s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		      s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		    }
		}
	      else		/* a particle */
		{
		  mass += 1;
		  s[0] += P[p].Pos[0];
		  s[1] += P[p].Pos[1];
		  s[2] += P[p].Pos[2];
		}
	    }
	}

      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;

      if(toplevel)
	Nodes[no].u.d.cost = 1;
      else
	Nodes[no].u.d.cost = 0;

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else				/* single particle or pseudo particle */
    {
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		DomainNextnode[last - (All.MaxPart + MaxNodes)] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < All.MaxPart)	/* only set it for single particles */
	Father[no] = father;
    }
}





void force_exchange_pseudodata(void)
{
  int i, ni;
  MPI_Status status;
  int level, sendTask, recvTask;

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      ni = DomainNodeIndex[i];

      /* read out the multipole moments from the local base cells */
      DomainMoment[i].s[0] = Nodes[ni].u.d.s[0];
      DomainMoment[i].s[1] = Nodes[ni].u.d.s[1];
      DomainMoment[i].s[2] = Nodes[ni].u.d.s[2];
      DomainMoment[i].mass = Nodes[ni].u.d.mass;
    }

  /* share the pseudo-particle data accross CPUs */

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainMoment[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(struct DomainNODE),
		     MPI_BYTE, recvTask, TAG_PSEUDO,
		     &DomainMoment[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(struct DomainNODE),
		     MPI_BYTE, recvTask, TAG_PSEUDO, MPI_COMM_WORLD, &status);
    }


  for(i = 0; i < DOMAINGRID * DOMAINGRID * DOMAINGRID; i++)
    {
      if(i < DomainMyStart || i > DomainMyLast)
	{
	  ni = DomainNodeIndex[i];

	  /* write the multipole moments into the local pseudo cells */
	  Nodes[ni].u.d.s[0] = DomainMoment[i].s[0];
	  Nodes[ni].u.d.s[1] = DomainMoment[i].s[1];
	  Nodes[ni].u.d.s[2] = DomainMoment[i].s[2];
	  Nodes[ni].u.d.mass = DomainMoment[i].mass;
	}
    }
}


void force_treeupdate_toplevel(int no, int level, int *nextfree)
{
  int i, j, k, nn;
  long long mass;
  double s[3];

  if(level < DOMAINLEVELS)
    {
      level += 1;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;

      for(k = -1; k <= 1; k += 2)
	for(j = -1; j <= 1; j += 2)
	  for(i = -1; i <= 1; i += 2)
	    {
	      nn = *nextfree;

	      *nextfree = *nextfree + 1;

	      force_treeupdate_toplevel(*nextfree - 1, level, nextfree);

	      mass += Nodes[nn].u.d.mass;
	      s[0] += Nodes[nn].u.d.mass * Nodes[nn].u.d.s[0];
	      s[1] += Nodes[nn].u.d.mass * Nodes[nn].u.d.s[1];
	      s[2] += Nodes[nn].u.d.mass * Nodes[nn].u.d.s[2];
	    }

      if(mass < (1 << 30))
	Nodes[no].u.d.mass = mass;
      else
	Nodes[no].u.d.mass = (1 << 30);	/* to prevent wrap-around in millennium run */

      if(mass)
	{
	  Nodes[no].u.d.s[0] = s[0] / mass;
	  Nodes[no].u.d.s[1] = s[1] / mass;
	  Nodes[no].u.d.s[2] = s[2] / mass;
	}
    }
}


int force_treeevaluate_shortrange(int target, int mode, FLOAT * acc)
{
  struct NODE *nop = 0;
  int no, ninteractions, tabindex, tabindex_iso;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z, aold;
  double eff_dist;
  double rcut, asmthfac, rcut2, dist;
  double asmth2;
  double isofac;
  double boxsize, boxhalf;
  double anifacx = All.TimeAni[0] / All.Time;
  double anifacy = All.TimeAni[1] / All.Time;
  double anifacz = All.TimeAni[2] / All.Time;
  double Jacobi = 1.0/(anifacx*anifacy*anifacz);
 // double iso = pow(Jacobi, 1./3.);
  double danix = anifacx - 1.0;
  double daniy = anifacy - 1.0;
  double daniz = anifacz - 1.0;
  double fac0, fac1, fac11, fac12, fac21, fac22, fac23, fac24, fac25, fac26, fac2x, fac2y, fac2z;
  double dI3, dI5, dI7, dI9, dI11;
  double I7, I9, I11;
  double deci = 0.0;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;

  /*
  if(ThisTask == 0)
   {
     printf("\niso = %lf\n", iso);
     printf("\ndanix = %lf\n", danix);
     printf("\ndaniy = %lf\n", daniy);
     printf("\ndaniz = %lf\n", daniz);
   }
   */


  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      aold = All.ErrTolForceAcc * P[target].OldAcc;
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
      aold = All.ErrTolForceAcc * GravDataGet[target].OldAcc;
    }


  rcut = All.Rcut[0];
  rcut2 = rcut * rcut;
  asmthfac = 0.5 / All.Asmth[0] * (NTAB / 3.0);
//  asmthfac /= iso;

  h = 2.8 * All.ComovSoftening;
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;


  no = All.MaxPart;		/* first node */

  while(no >= 0)
    {
      if(no < All.MaxPart)
	{
	  /* the index of the node is the index of the particle */
	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;

	  dx = anifacx * NEAREST(dx);
	  dy = anifacy * NEAREST(dy);
	  dz = anifacz * NEAREST(dz);

	  r2 = dx * dx + dy * dy + dz * dz;

	  mass = 1;
	  /*
	     P[no].GravCost += 1.0;
	   */
	  no = Nextnode[no];
	}
      else			/* we have an  internal node */
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
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

	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;

	  dx = anifacx * NEAREST(dx);
	  dy = anifacy * NEAREST(dy);
	  dz = anifacz * NEAREST(dz);

	  r2 = dx * dx + dy * dy + dz * dz;

	  if(r2 > rcut2)
	    {
	      /* check whether we can stop walking along this branch */
	      eff_dist = rcut + 0.5 * nop->len;

	      dist = NEAREST(nop->center[0] - pos_x);
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}

	      dist = NEAREST(nop->center[1] - pos_y);
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}

	      dist = NEAREST(nop->center[2] - pos_z);
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  mass = nop->u.d.mass;

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* check in addition whether we lie inside the cell */

	      if( fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if( fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if( fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if((nop->u.d.cost & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }

	  nop->u.d.cost += 4;
	}

      r = sqrt(r2);

      if(r >= h)
	fac = mass / (r2 * r); //needs Jacobian?
      else
	{
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}

      tabindex = (int) (asmthfac * r);
      deci = asmthfac * r - tabindex;

    if(r >= h)
    {
	  fac1 = mass / (2.*sqrt(M_PI)); // needs Jacobian?
	  fac1 /= All.Asmth[0];
	} else 
	{
	  fac1 = 0.0;
	}

      if(tabindex < NTAB-1)
	{
	  //fac *= ( (1.0 - deci)*shortrange_table[tabindex] + deci*shortrange_table[tabindex+1] );

	  acc_x += dx * fac * anifacx;
	  acc_y += dy * fac * anifacy;
	  acc_z += dz * fac * anifacz;

	if( r >= h){
	  dI3 = ( (1.0 - deci) * dI3_table[tabindex] + deci*dI3_table[tabindex+1] ) * (2.0/r);

	  // 0th order correction of Delta_i (Jacobian)

	  fac0 = dI3 / r;

	  acc_x += dx * fac1 * anifacx * fac0;
	  acc_y += dy * fac1 * anifacy * fac0;
	  acc_z += dz * fac1 * anifacz * fac0;

	  if( ( fabs(danix) > 0.0001) || (fabs(daniy) > 0.0001) || (fabs(daniz) > 0.0001) )
	  {

	  	dI5 = ( (1.0 - deci) * dI5_table[tabindex] + deci*dI5_table[tabindex+1] ) * (2.0/r);

	  	dI7 = ( (1.0 - deci) * dI7_table[tabindex] + deci*dI7_table[tabindex+1] ) * (2.0/r);

	  	dI9 = ( (1.0 - deci) * dI9_table[tabindex] + deci*dI9_table[tabindex+1] ) * (2.0/r);

	  	dI11 = ( (1.0 - deci) * dI11_table[tabindex] + deci*dI11_table[tabindex+1] ) * (2.0/r);

	  	I7 =  (1.0 - deci) * I7_table[tabindex] + deci*I7_table[tabindex+1];

	  	I9 =  (1.0 - deci) * I9_table[tabindex] + deci*I9_table[tabindex+1];

	  	I11 =  (1.0 - deci) * I11_table[tabindex] + deci*I11_table[tabindex+1];

	
	  // 1st order correction of Delta_i 

	  asmth2 = All.Asmth[0] * All.Asmth[0];
	  fac11 = - dI5 * ( danix + daniy + daniz );
	  fac11 += dI7 * (danix*dx*dx + daniy*dy*dy + daniz*dz*dz) / (2.0 * asmth2);
	  fac11 /= r;
	  fac12 = I7 / asmth2;

	  acc_x += dx * fac1 * anifacx * ( fac11 + fac12 * danix );
	  acc_y += dy * fac1 * anifacy * ( fac11 + fac12 * daniy );
	  acc_z += dz * fac1 * anifacz * ( fac11 + fac12 * daniz );

	  // 2nd order correction of Delta_i
	  //Im derivative terms
	  // i=j terms
	  fac21 = ( danix*danix + daniy*daniy + daniz*daniz );
	  fac21 *= ( -dI5 + 3.0 * dI7 );

	  fac22 = dx*dx*danix*danix + dy*dy*daniy*daniy + dz*dz*daniz*daniz;
	  fac22 *= ( dI7/2.0 - 3.0 * dI9 ) ;
	  fac22 /= asmth2;

	  fac23 = danix*danix*dx*dx*dx*dx + daniy*daniy*dy*dy*dy*dy + daniz*daniz*dz*dz*dz*dz;
	  fac23 *= dI11;
	  fac23 /= (4.0 * asmth2*asmth2);

	  // i \neq j terms
	  fac24 =  danix*daniy * ( dI7 - (dx*dx + dy*dy)/(2.0*asmth2) * dI9 + dI11 * dx*dx*dy*dy/(4.0*asmth2*asmth2) );
	  fac24 += daniy*daniz * ( dI7 - (dy*dy + dz*dz)/(2.0*asmth2) * dI9 + dI11 * dy*dy*dz*dz/(4.0*asmth2*asmth2) );
	  fac24 += daniz*danix * ( dI7 - (dz*dz + dx*dx)/(2.0*asmth2) * dI9 + dI11 * dz*dz*dx*dx/(4.0*asmth2*asmth2) );
	  fac24 *= 2.0;

	  //f derivative terms
	  // i=j terms
	  fac25 = ( I7/2.0 - 3.0 * I9  ) / asmth2; 
	  fac26 = I11 / (2.0*asmth2*asmth2);

	  // i \neq j terms
	  fac2x = ( -I9 * (daniy + daniz) / (2.0*asmth2) + I11 * ( daniy*dy*dy + daniz*dz*dz ) / (4.0*asmth2*asmth2) );
	  fac2y = ( -I9 * (daniz + danix) / (2.0*asmth2) + I11 * ( daniz*dz*dz + danix*dx*dx ) / (4.0*asmth2*asmth2) );
	  fac2z = ( -I9 * (danix + daniy) / (2.0*asmth2) + I11 * ( danix*dx*dx + daniy*dy*dy ) / (4.0*asmth2*asmth2) );

	  acc_x += dx * fac1 * anifacx * ( fac21 + fac22 + fac23 + fac24  ) / (2.0 * r);
	  acc_x += dx * fac1 * anifacx * danix * ( danix * ( fac25 + fac26 * dx*dx )  + fac2x );
	  acc_y += dy * fac1 * anifacy * ( fac21 + fac22 + fac23 + fac24  ) / (2.0 * r);
	  acc_y += dy * fac1 * anifacy * daniy * ( daniy * ( fac25 + fac26 * dy*dy )  + fac2y );
	  acc_z += dz * fac1 * anifacz * ( fac21 + fac22 + fac23 + fac24  ) / (2.0 * r);
	  acc_z += dz * fac1 * anifacz * daniz * ( daniz * ( fac25 + fac26 * dz*dz )  + fac2z );

	  }

	}
	  ninteractions++;
	}
    }

  /* store result at the proper place */

  if(mode == 0)
    {
      acc[0] = acc_x;
      acc[1] = acc_y;
      acc[2] = acc_z;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
    }

  return ninteractions;
}







void force_costevaluate(void)
{
  int i, no;

  for(i = 0; i < NumPart; i++)
    {
      no = Father[i];
      while(no >= All.MaxPart)
	{
	  P[i].GravCost += (Nodes[no].u.d.cost >> 2) / Nodes[no].u.d.mass;

	  no = Nodes[no].u.d.father;
	}
    }
}


void force_flag_localnodes(void)
{
  int no, i;

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      while(no >= 0)
	{
	  if((Nodes[no].u.d.cost & 2))
	    break;

	  Nodes[no].u.d.cost |= 2;

	  no = Nodes[no].u.d.father;
	}
    }
}


double dIm_func(int m, double u)
{
	double gammafac, gamma;
    if( u > 0. )
    {
        gammafac = (m-2) / 2.0 ;
        gammafac *= pow(u, 2-m);
        gamma = gsl_sf_gamma_inc((m-2)/2.0, 0) - gsl_sf_gamma_inc((m-2)/2.0, u*u);
        return exp(-u * u) - gammafac * gamma;
    } else {
        return 0.0;
    }
}

double Im_func(int m, double u)
{
    if( u > 0. )
    {
        double prefac, gamma;
        prefac = pow(1./u, m-2);
        gamma = gsl_sf_gamma_inc((m-2)/2.0, 0) - gsl_sf_gamma_inc((m-2)/2.0, u*u);
        return prefac * gamma;
    } else {
        return 2.0 / (m-2);
    }
}


void force_treeinit(void)
{
  int i;
  double u;

  for(i = 0; i < NTAB; i++)
    {
      u = 3.0 / NTAB * i;
      //shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
      dI3_table[i] = dIm_func(3, u);
      dI5_table[i] = dIm_func(5, u);
      dI7_table[i] = dIm_func(7, u);
      dI9_table[i] = dIm_func(9, u);
      dI11_table[i] = dIm_func(11, u);
      I7_table[i]  = Im_func(7, u);
      I9_table[i]  = Im_func(9, u);
      I11_table[i] = Im_func(11, u);
    }
    if(ThisTask == 0)
    {
      printf("\nforce_treeinit finished.\n");
    }
}


/* this function allocates memory used for storage of the tree
 * and auxiliary arrays for tree-walk and link-lists.
 */
size_t force_treeallocate(int maxnodes, int maxpart)	/* usually maxnodes=0.7*maxpart is sufficient */
{
  size_t bytes, allbytes = 0;

  MaxNodes = maxnodes;

  if(!(Nodes_base = mymalloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  if(!(Nextnode = mymalloc(bytes = (maxpart) * sizeof(int4byte))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n", maxpart,
	     bytes / (1024.0 * 1024.0));
      endrun(4);
    }
  allbytes += bytes;

  if(!(Father = mymalloc(bytes = (maxpart) * sizeof(int4byte))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      endrun(5);
    }
  allbytes += bytes;

  return allbytes;
}


/* free the allocated memory
 */
void force_treefree(void)
{
  myfree(Father);
  myfree(Nextnode);
  myfree(Nodes_base);
}
