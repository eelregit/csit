#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#ifdef HPM
#include <libhpm.h>
#endif


/*! \file pm_periodic.c
 *  \brief routines for periodic PM force computation
 */

#ifdef NOTYPEPREFIX_FFTW
#include	<rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif

#include "allvars.h"
#include "proto.h"

#define  PMGRID2   (2*(PMGRID/2 + 1))


static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;

static int slab_to_task[PMGRID];
static int *slabs_per_task;
static int *first_slab_of_task;
static int *meshmin_list, *meshmax_list;

static int slabstart_x, nslab_x, slabstart_y, nslab_y;

static int fftsize, maxfftsize, dimprodmax;

static fftw_complex *fft_of_rhogrid;
static fftw_real *rhogrid, *forcegrid, *workspace;

static double to_slab_fac;


void pm_init_periodic(void)
{
  int i;
  int slab_to_task_local[PMGRID];


  All.Asmth[0] = ASMTH * All.BoxSize / PMGRID;
  All.Rcut[0] = RCUT * All.Asmth[0];

  /* Set up the FFTW plan files. */

  fft_forward_plan =
    rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID, FFTW_REAL_TO_COMPLEX,
			    FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan =
    rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID, FFTW_COMPLEX_TO_REAL,
			    FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  slabs_per_task = mymalloc(NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = mymalloc(NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  meshmin_list = mymalloc(3 * NTask * sizeof(int));
  meshmax_list = mymalloc(3 * NTask * sizeof(int));


  to_slab_fac = PMGRID / All.BoxSize;

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

}





/* Calculate the force given the density field using the PM method.  
 * The force is Gaussian filtered with asmth, given in box units. 
 */
void pmforce_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz, acc_dim;
  double fx, fy, fz, ff;
  double asmth2, fac, fac2;
  double dt_gravkick;
  double atime, fac1, hubble_a;
  double mbyte, mbytesum;
  double t0, t1;
  int ti_step, tend, tstart;
  int i, j, slab, level, sendTask, recvTask;
  int x, y, z, xl, yl, zl, xr, yr, zr, xll, yll, zll, xrr, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int meshmin[3], meshmax[3], sendmin, sendmax, recvmin, recvmax;
  int dimx, dimy, dimz, recv_dimx, recv_dimy, recv_dimz;
  int rep, ncont, cont_sendmin[2], cont_sendmax[2], cont_recvmin[2], cont_recvmax[2];
  int dimprod;
  size_t bytes_tot = 0, bytes;
  MPI_Status status;
  double anifac = (All.Time*All.Time*All.Time) / (All.TimeAni[0]*All.TimeAni[1]*All.TimeAni[2]);
  double anifacx = (All.Time*All.Time) / (All.TimeAni[0]*All.TimeAni[0]);
  double anifacy = (All.Time*All.Time) / (All.TimeAni[1]*All.TimeAni[1]);
  double anifacz = (All.Time*All.Time) / (All.TimeAni[2]*All.TimeAni[2]);

  if(ThisTask == 0)
    {
      printf("Starting periodic PM calculation.\n");
      fflush(stdout);
    }

  ti_step = TIMEBASE;
  while(ti_step > (All.Dt_displacement / All.Timebase_interval))
    ti_step >>= 1;

  if(ti_step > (All.PM_Ti_endstep - All.PM_Ti_begstep))	/* PM-timestep wants to increase */
    {
      /* we only increase if an integer number of steps will bring us to the end */
      if(((TIMEBASE - All.PM_Ti_endstep) % ti_step) > 0)
	ti_step = All.PM_Ti_endstep - All.PM_Ti_begstep;	/* leave at old step */
    }

  if(All.Ti_Current == TIMEBASE)	/* we here finish the last timestep. */
    ti_step = 0;

  tstart = (All.PM_Ti_begstep + All.PM_Ti_endstep) / 2;
  tend = All.PM_Ti_endstep + ti_step / 2;


  dt_gravkick = get_gravkick_factor(tstart, tend);

  All.PM_Ti_begstep = All.PM_Ti_endstep;
  All.PM_Ti_endstep = All.PM_Ti_begstep + ti_step;

  fac1 = 1 / (All.Time * All.Time);
  atime = All.Time;
  hubble_a = hubble_function(All.Time);


  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = 1.0 / (M_PI * All.BoxSize);	/* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */
#ifdef MAKEGLASS
  dt_gravkick = 1;
#endif
  fac2 = dt_gravkick * All.G * All.PartMass;



  /* first, establish the extension of the local patch in the PMGRID  */

  for(j = 0; j < 3; j++)
    {
      meshmin[j] = PMGRID;
      meshmax[j] = 0;
    }

  for(i = 0; i < NumPart; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  slab = to_slab_fac * P[i].Pos[j];
	  if(slab >= PMGRID)
	    slab = PMGRID - 1;

	  if(slab < meshmin[j])
	    meshmin[j] = slab;

	  if(slab > meshmax[j])
	    meshmax[j] = slab;
	}
    }

  MPI_Allgather(meshmin, 3, MPI_INT, meshmin_list, 3, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(meshmax, 3, MPI_INT, meshmax_list, 3, MPI_INT, MPI_COMM_WORLD);

  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  dimprod = (dimx + 4) * (dimy + 4) * (dimz + 4);
  MPI_Allreduce(&dimprod, &dimprodmax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  /* now allocate memory to hold the FFT fields */

  if(!(rhogrid = (fftw_real *) mymalloc(bytes = fftsize * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  fft_of_rhogrid = (fftw_complex *) & rhogrid[0];

  if(!(forcegrid = (fftw_real *) mymalloc(bytes = imax(fftsize, dimprodmax) * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-forcegrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(workspace = (fftw_real *) mymalloc(bytes = imax(maxfftsize, dimprodmax) * sizeof(fftw_real))))
    {
      printf("failed to allocate memory for `FFT-workspace' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  mbyte = bytes_tot / (1024.0 * 1024.0);
  MPI_Reduce(&mbyte, &mbytesum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Use of %g MByte for PM calculation.\n", mbytesum / NTask);
      fflush(stdout);
    }

#ifdef POWERSPEC
  if(All.PowerSpecFlag)
    {
      foldonitself();
      powerspec(0);
    }
#endif

  for(i = 0; i < dimx * dimy * dimz; i++)
    workspace[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      slab_x = to_slab_fac * P[i].Pos[0];
      if(slab_x >= PMGRID)
	slab_x = PMGRID - 1;
      dx = to_slab_fac * P[i].Pos[0] - slab_x;
      slab_x -= meshmin[0];
      slab_xx = slab_x + 1;

      slab_y = to_slab_fac * P[i].Pos[1];
      if(slab_y >= PMGRID)
	slab_y = PMGRID - 1;
      dy = to_slab_fac * P[i].Pos[1] - slab_y;
      slab_y -= meshmin[1];
      slab_yy = slab_y + 1;

      slab_z = to_slab_fac * P[i].Pos[2];
      if(slab_z >= PMGRID)
	slab_z = PMGRID - 1;
      dz = to_slab_fac * P[i].Pos[2] - slab_z;
      slab_z -= meshmin[2];
      slab_zz = slab_z + 1;

      workspace[(slab_x * dimy + slab_y) * dimz + slab_z] += (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_z] += (1.0 - dx) * dy * (1.0 - dz);
      workspace[(slab_x * dimy + slab_y) * dimz + slab_zz] += (1.0 - dx) * (1.0 - dy) * dz;
      workspace[(slab_x * dimy + slab_yy) * dimz + slab_zz] += (1.0 - dx) * dy * dz;

      workspace[(slab_xx * dimy + slab_y) * dimz + slab_z] += (dx) * (1.0 - dy) * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_z] += (dx) * dy * (1.0 - dz);
      workspace[(slab_xx * dimy + slab_y) * dimz + slab_zz] += (dx) * (1.0 - dy) * dz;
      workspace[(slab_xx * dimy + slab_yy) * dimz + slab_zz] += (dx) * dy * dz;
    }

  if(ThisTask == 0)
    {
      printf("start sharing density field.\n");
      fflush(stdout);
    }
  t0 = second();

  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;
      if(recvTask < NTask)
	{
	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -1;
	  for(slab_x = meshmin[0]; slab_x < meshmax[0] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == recvTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -1)
	    sendmin = 0;

	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -1;
	  for(slab_x = meshmin_list[3 * recvTask]; slab_x < meshmax_list[3 * recvTask] + 2; slab_x++)
	    if(slab_to_task[slab_x % PMGRID] == sendTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -1)
	    recvmin = 0;


	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 2;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 2;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 2;

	      if(level > 0)
		{
		  MPI_Sendrecv(workspace + (sendmin - meshmin[0]) * dimy * dimz,
			       (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE, recvTask,
			       TAG_PM_A, forcegrid,
			       (recvmax - recvmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_PM_A, MPI_COMM_WORLD, &status);
		}
	      else
		{
		  memcpy(forcegrid, workspace + (sendmin - meshmin[0]) * dimy * dimz,
			 (sendmax - sendmin + 1) * dimy * dimz * sizeof(fftw_real));
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  slab_xx = (slab_x % PMGRID) - first_slab_of_task[ThisTask];

		  if(slab_xx >= 0 && slab_xx < slabs_per_task[ThisTask])
		    {
		      for(slab_y = meshmin_list[3 * recvTask + 1];
			  slab_y <= meshmax_list[3 * recvTask + 1] + 1; slab_y++)
			{
			  slab_yy = slab_y;
			  if(slab_yy >= PMGRID)
			    slab_yy -= PMGRID;

			  for(slab_z = meshmin_list[3 * recvTask + 2];
			      slab_z <= meshmax_list[3 * recvTask + 2] + 1; slab_z++)
			    {
			      slab_zz = slab_z;
			      if(slab_zz >= PMGRID)
				slab_zz -= PMGRID;

			      rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz] +=
				forcegrid[((slab_x - recvmin) * recv_dimy +
					   (slab_y - meshmin_list[3 * recvTask + 1])) * recv_dimz +
					  (slab_z - meshmin_list[3 * recvTask + 2])];
			    }
			}
		    }
		}
	    }
	}
    }

  t1 = second();
  if(ThisTask == 0)
    {
      printf("finished sharing density field. (took %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }


  /* Do the FFT of the density field */

  t0 = second();

#ifdef HPM
  hpmStart(4, "FFT");
#endif
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#ifdef HPM
  hpmStop(4);
#endif

  t1 = second();
  if(ThisTask == 0)
    {
      printf("finished forward FFT (took %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }

#ifdef POWERSPEC
  if(All.PowerSpecFlag)
    powerspec(1);
#endif

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = anifacx * kx*kx + anifacy * ky*ky + anifacz * kz*kz;

	  if(k2 > 0)
	    {
	      smth = - anifac * exp(- (kx*kx + ky*ky + kz*kz) * asmth2) / k2;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;

	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      fft_of_rhogrid[ip].re *= smth;
	      fft_of_rhogrid[ip].im *= smth;
	    }
	}

  if(slabstart_y == 0)
    fft_of_rhogrid[0].re = fft_of_rhogrid[0].im = 0.0;

  /* Do the FFT to get the potential */


  t0 = second();

#ifdef HPM
  hpmStart(4, "FFT");
#endif
  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#ifdef HPM
  hpmStop(4);
#endif


  t1 = second();
  if(ThisTask == 0)
    {
      printf("finished backwards FFT (took %g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }

#ifdef POWERSPEC
  if(All.PowerSpecFlag)
    {
#ifdef OUTPUT_LONGRANGE_POTENTIAL
      dump_potential();
#endif
      All.PowerSpecFlag = 0;
    }
#endif

  /* Now rhogrid holds the potential */
  /* construct the potential for the local patch */

  t0 = second();

  dimx = meshmax[0] - meshmin[0] + 6;
  dimy = meshmax[1] - meshmin[1] + 6;
  dimz = meshmax[2] - meshmin[2] + 6;

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{

	  /* check how much we have to send */
	  sendmin = 2 * PMGRID;
	  sendmax = -PMGRID;
	  for(slab_x = meshmin_list[3 * recvTask] - 2; slab_x < meshmax_list[3 * recvTask] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == sendTask)
	      {
		if(slab_x < sendmin)
		  sendmin = slab_x;
		if(slab_x > sendmax)
		  sendmax = slab_x;
	      }
	  if(sendmax == -PMGRID)
	    sendmin = sendmax + 1;


	  /* check how much we have to receive */
	  recvmin = 2 * PMGRID;
	  recvmax = -PMGRID;
	  for(slab_x = meshmin[0] - 2; slab_x < meshmax[0] + 4; slab_x++)
	    if(slab_to_task[(slab_x + PMGRID) % PMGRID] == recvTask)
	      {
		if(slab_x < recvmin)
		  recvmin = slab_x;
		if(slab_x > recvmax)
		  recvmax = slab_x;
	      }
	  if(recvmax == -PMGRID)
	    recvmin = recvmax + 1;

	  if((recvmax - recvmin) >= 0 || (sendmax - sendmin) >= 0)	/* ok, we have a contribution to the slab */
	    {
	      recv_dimx = meshmax_list[3 * recvTask + 0] - meshmin_list[3 * recvTask + 0] + 6;
	      recv_dimy = meshmax_list[3 * recvTask + 1] - meshmin_list[3 * recvTask + 1] + 6;
	      recv_dimz = meshmax_list[3 * recvTask + 2] - meshmin_list[3 * recvTask + 2] + 6;

	      ncont = 1;
	      cont_sendmin[0] = sendmin;
	      cont_sendmax[0] = sendmax;
	      cont_sendmin[1] = sendmax + 1;
	      cont_sendmax[1] = sendmax;

	      cont_recvmin[0] = recvmin;
	      cont_recvmax[0] = recvmax;
	      cont_recvmin[1] = recvmax + 1;
	      cont_recvmax[1] = recvmax;

	      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
		    {
		      /* non-contiguous */
		      cont_sendmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != ThisTask)
			slab_x++;
		      cont_sendmin[1] = slab_x;
		      ncont++;
		    }
		}

	      for(slab_x = recvmin; slab_x <= recvmax; slab_x++)
		{
		  if(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
		    {
		      /* non-contiguous */
		      cont_recvmax[0] = slab_x - 1;
		      while(slab_to_task[(slab_x + PMGRID) % PMGRID] != recvTask)
			slab_x++;
		      cont_recvmin[1] = slab_x;
		      if(ncont == 1)
			ncont++;
		    }
		}


	      for(rep = 0; rep < ncont; rep++)
		{
		  sendmin = cont_sendmin[rep];
		  sendmax = cont_sendmax[rep];
		  recvmin = cont_recvmin[rep];
		  recvmax = cont_recvmax[rep];

		  /* prepare what we want to send */
		  if(sendmax - sendmin >= 0)
		    {
		      for(slab_x = sendmin; slab_x <= sendmax; slab_x++)
			{
			  slab_xx = ((slab_x + PMGRID) % PMGRID) - first_slab_of_task[ThisTask];

			  for(slab_y = meshmin_list[3 * recvTask + 1] - 2;
			      slab_y < meshmax_list[3 * recvTask + 1] + 4; slab_y++)
			    {
			      slab_yy = (slab_y + PMGRID) % PMGRID;

			      for(slab_z = meshmin_list[3 * recvTask + 2] - 2;
				  slab_z < meshmax_list[3 * recvTask + 2] + 4; slab_z++)
				{
				  slab_zz = (slab_z + PMGRID) % PMGRID;

				  forcegrid[((slab_x - sendmin) * recv_dimy +
					     (slab_y - (meshmin_list[3 * recvTask + 1] - 2))) * recv_dimz +
					    slab_z - (meshmin_list[3 * recvTask + 2] - 2)] =
				    rhogrid[PMGRID * PMGRID2 * slab_xx + PMGRID2 * slab_yy + slab_zz];
				}
			    }
			}
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(forcegrid,
				   (sendmax - sendmin + 1) * recv_dimy * recv_dimz * sizeof(fftw_real),
				   MPI_BYTE, recvTask, TAG_PM_B,
				   workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
				   (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PM_B, MPI_COMM_WORLD, &status);
		    }
		  else
		    {
		      memcpy(workspace + (recvmin - (meshmin[0] - 2)) * dimy * dimz,
			     forcegrid, (recvmax - recvmin + 1) * dimy * dimz * sizeof(fftw_real));
		    }
		}
	    }
	}
    }

  t1 = second();
  if(ThisTask == 0)
    {
      printf("finished assembling local potential (%g sec)\n", timediff(t0, t1));
      fflush(stdout);
    }

  dimx = meshmax[0] - meshmin[0] + 2;
  dimy = meshmax[1] - meshmin[1] + 2;
  dimz = meshmax[2] - meshmin[2] + 2;

  recv_dimx = meshmax[0] - meshmin[0] + 6;
  recv_dimy = meshmax[1] - meshmin[1] + 6;
  recv_dimz = meshmax[2] - meshmin[2] + 6;

  for(i = 0; i < NumPart; i++)
    P[i].GravCost = 0;

  for(dim = 0; dim < 3; dim++)	/* Calculate each component of the force. */
    {
      /* get the force component by finite differencing the potential */
      /* note: "workspace" now contains the potential for the local patch, plus a suffiently large buffer region */

      for(x = 0; x < meshmax[0] - meshmin[0] + 2; x++)
	for(y = 0; y < meshmax[1] - meshmin[1] + 2; y++)
	  for(z = 0; z < meshmax[2] - meshmin[2] + 2; z++)
	    {
	      xrr = xll = xr = xl = x;
	      yrr = yll = yr = yl = y;
	      zrr = zll = zr = zl = z;

	      switch (dim)
		{
		case 0:
		  xr = x + 1;
		  xrr = x + 2;
		  xl = x - 1;
		  xll = x - 2;
		  break;
		case 1:
		  yr = y + 1;
		  yl = y - 1;
		  yrr = y + 2;
		  yll = y - 2;
		  break;
		case 2:
		  zr = z + 1;
		  zl = z - 1;
		  zrr = z + 2;
		  zll = z - 2;
		  break;
		}

	      forcegrid[(x * dimy + y) * dimz + z]
		=
		fac * ((4.0 / 3) *
		       (workspace[((xl + 2) * recv_dimy + (yl + 2)) * recv_dimz + (zl + 2)]
			- workspace[((xr + 2) * recv_dimy + (yr + 2)) * recv_dimz + (zr + 2)]) -
		       (1.0 / 6) *
		       (workspace[((xll + 2) * recv_dimy + (yll + 2)) * recv_dimz + (zll + 2)] -
			workspace[((xrr + 2) * recv_dimy + (yrr + 2)) * recv_dimz + (zrr + 2)]));
	    }


      /* read out the forces and do the PM kick */

      for(i = 0; i < NumPart; i++)
	{
	  slab_x = to_slab_fac * P[i].Pos[0];
	  if(slab_x >= PMGRID)
	    slab_x = PMGRID - 1;
	  dx = to_slab_fac * P[i].Pos[0] - slab_x;
	  slab_x -= meshmin[0];
	  slab_xx = slab_x + 1;

	  slab_y = to_slab_fac * P[i].Pos[1];
	  if(slab_y >= PMGRID)
	    slab_y = PMGRID - 1;
	  dy = to_slab_fac * P[i].Pos[1] - slab_y;
	  slab_y -= meshmin[1];
	  slab_yy = slab_y + 1;

	  slab_z = to_slab_fac * P[i].Pos[2];
	  if(slab_z >= PMGRID)
	    slab_z = PMGRID - 1;
	  dz = to_slab_fac * P[i].Pos[2] - slab_z;
	  slab_z -= meshmin[2];
	  slab_zz = slab_z + 1;

	  acc_dim =
	    forcegrid[(slab_x * dimy + slab_y) * dimz + slab_z] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_z] * (1.0 - dx) * dy * (1.0 - dz);
	  acc_dim += forcegrid[(slab_x * dimy + slab_y) * dimz + slab_zz] * (1.0 - dx) * (1.0 - dy) * dz;
	  acc_dim += forcegrid[(slab_x * dimy + slab_yy) * dimz + slab_zz] * (1.0 - dx) * dy * dz;

	  acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_z] * (dx) * (1.0 - dy) * (1.0 - dz);
	  acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_z] * (dx) * dy * (1.0 - dz);
	  acc_dim += forcegrid[(slab_xx * dimy + slab_y) * dimz + slab_zz] * (dx) * (1.0 - dy) * dz;
	  acc_dim += forcegrid[(slab_xx * dimy + slab_yy) * dimz + slab_zz] * (dx) * dy * dz;

	  P[i].Vel[dim] += fac2 * acc_dim;

	  P[i].GravCost += acc_dim * acc_dim;	/* we use GravCost here to temporarily store the
						   magnitude of the PM acceleration */
	}
    }


  myfree(workspace);
  myfree(forcegrid);
  myfree(rhogrid);

  if(ThisTask == 0)
    {
      printf("done PM.\n");
      fflush(stdout);
    }
}



/* Here comes code for the power-sepctrum computation.
 */

#define BINS_PS  2000		/* number of bins for power spectrum computation */
#define POWERSPEC_FOLDFAC 32


static long long CountModes[2][BINS_PS];
static double SumPower[2][BINS_PS];
static double SumPowerUncorrected[2][BINS_PS];	/* without binning correction (as for shot noise) */
static double Power[2][BINS_PS];
static double PowerUncorrected[2][BINS_PS];	/* without binning correction */
static double Delta[2][BINS_PS];
static double DeltaUncorrected[2][BINS_PS];	/* without binning correction */
static double ShotLimit[2][BINS_PS];
static double Kbin[BINS_PS];
static double K0, K1;
static double binfac;


//#define RAYLEIGH_BINS 100

#ifdef RAYLEIGH_BINS
static int RayleighCountModes[BINS_PS][RAYLEIGH_BINS];
static int RayleighCountAbove[BINS_PS];
static double RayleighFactor[BINS_PS];
static float RayleighMax[BINS_PS];
#endif

void powerspec(int flag)
{
  int i, n, x, y, z, kx, ky, kz, bin, ip, rep, zz;
  double k, k2, po, ponorm, smth, fac;
  double fx, fy, fz, ff;
  double *powerbuf;
  long long *countbuf;
  double tstart, tend;

#ifdef RAYLEIGH_BINS
  double ratio;
#endif

  if(ThisTask == 0)
    {
      printf("begin power spectrum. (step=%d)\n", flag);
      fflush(stdout);
    }

  tstart = second();

  fac = 1.0 / All.TotNumPart;	/* to convert to density fluctuation,
				   we need the factor
				   PMGRID^3/TotNumPart. However, the
				   FFT is defined with the factor
				   1/PMGRID^3 */

  K0 = 2 * M_PI / All.BoxSize;	/* minimum k */
  K1 = K0 * All.BoxSize / All.Softening;	/* maximum k */
  binfac = BINS_PS / (log(K1) - log(K0));

  if(flag == 0)
    {
      for(rep = 0; rep < 2; rep++)
	for(i = 0; i < BINS_PS; i++)
	  {
	    SumPower[rep][i] = 0;
	    SumPowerUncorrected[rep][i] = 0;
	    CountModes[rep][i] = 0;
	  }

#ifdef RAYLEIGH_BINS
      for(i = 0; i < BINS_PS; i++)
	{
	  RayleighCountAbove[i] = 0;
	  RayleighMax[i] = 0;

	  for(n = 0; n < RAYLEIGH_BINS; n++)
	    RayleighCountModes[i][n] = 0;
	}
      prepare_rayleigh();
#endif
    }

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID; z++)
	{
	  zz = z;
	  if(z >= PMGRID / 2 + 1)
	    zz = PMGRID - z;

	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      if(k2 < (PMGRID / 2.0) * (PMGRID / 2.0))
		{
		  /* do deconvolution */

		  fx = fy = fz = 1;
		  if(kx != 0)
		    {
		      fx = (M_PI * kx) / PMGRID;
		      fx = sin(fx) / fx;
		    }
		  if(ky != 0)
		    {
		      fy = (M_PI * ky) / PMGRID;
		      fy = sin(fy) / fy;
		    }
		  if(kz != 0)
		    {
		      fz = (M_PI * kz) / PMGRID;
		      fz = sin(fz) / fz;
		    }
		  ff = 1 / (fx * fy * fz);
		  smth = ff * ff * ff * ff;

		  /* end deconvolution */

		  ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + zz;

		  po = (fft_of_rhogrid[ip].re * fft_of_rhogrid[ip].re
			+ fft_of_rhogrid[ip].im * fft_of_rhogrid[ip].im);

		  po *= fac * fac * smth;

		  k = sqrt(k2) * 2 * M_PI / All.BoxSize;

		  if(flag == 0)
		    k *= POWERSPEC_FOLDFAC;

		  if(k >= K0 && k < K1)
		    {
		      bin = log(k / K0) * binfac;

		      ponorm = po / PowerSpec_Efstathiou(k);

		      SumPower[flag][bin] += ponorm;
		      SumPowerUncorrected[flag][bin] += po;
		      CountModes[flag][bin] += 1;

#ifdef RAYLEIGH_BINS
		      if(flag == 1)	/* only for course grid */
			{
			  ratio = sqrt(ponorm * RayleighFactor[bin]);

			  if(ratio > RayleighMax[bin])
			    RayleighMax[bin] = ratio;

			  if(ratio >= 10.0)
			    RayleighCountAbove[bin]++;
			  else
			    {
			      i = RAYLEIGH_BINS * ratio / 10.0;
			      RayleighCountModes[bin][i]++;
			    }
			}
#endif
		    }
		}
	    }
	}

  /* Now compute the power spectrum */

  countbuf = mymalloc(NTask * BINS_PS * sizeof(long long));
  powerbuf = mymalloc(NTask * BINS_PS * sizeof(double));

  MPI_Allgather(CountModes[flag], BINS_PS * sizeof(long long), MPI_BYTE,
		countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      CountModes[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	CountModes[flag][i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPower[flag], BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	SumPower[flag][i] += powerbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPowerUncorrected[flag], BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPowerUncorrected[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	SumPowerUncorrected[flag][i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);


  for(i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[flag][i] > 0)
	{
	  Power[flag][i] = PowerSpec_Efstathiou(Kbin[i]) * SumPower[flag][i] / CountModes[flag][i];
	  PowerUncorrected[flag][i] = SumPowerUncorrected[flag][i] / CountModes[flag][i];
	}
      else
	{
	  Power[flag][i] = 0;
	  PowerUncorrected[flag][i] = 0;
	}

      Delta[flag][i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * Power[flag][i];

      DeltaUncorrected[flag][i] = 4 * M_PI * pow(Kbin[i], 3) /
	pow(2 * M_PI / All.BoxSize, 3) * PowerUncorrected[flag][i];

      ShotLimit[flag][i] = 4 * M_PI * pow(Kbin[i], 3) /
	pow(2 * M_PI / All.BoxSize, 3) * (1.0 / All.TotNumPart);
    }

  if(flag == 1)
    {
      powerspec_save();
#ifdef RAYLEIGH_BINS
      rayleigh_save();
#endif
    }

  tend = second();
  All.CPU_PowerSpec += timediff(tstart, tend);

  if(ThisTask == 0)
    {
      printf("end power spectrum. (step=%d) took %g seconds\n", flag, timediff(tstart, tend));
      fflush(stdout);
    }
}


double PowerSpec_Efstathiou(double k)
{
  double AA, BB, CC, nu, ShapeGamma;

  ShapeGamma = 0.21;
  AA = 6.4 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  BB = 3.0 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  CC = 1.7 / ShapeGamma * (3.085678e24 / All.UnitLength_in_cm);
  nu = 1.13;


  return k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}




void powerspec_save(void)
{
  FILE *fd;
  char buf[500];
  int i, flag;

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/powerspec_%03d.txt", All.OutputDir, All.PowerSpecFlag - 1);

      if(!(fd = fopen(buf, "w")))
	{
	  printf("can't open file `%s`\n", buf);
	  endrun(1324);
	}
      for(flag = 0; flag < 2; flag++)
	{
	  fprintf(fd, "%g\n", All.Time);
	  i = BINS_PS;
	  fprintf(fd, "%d\n", i);

	  for(i = 0; i < BINS_PS; i++)
	    {
	      fprintf(fd, "%g %g %g %g %g %g %g %g %g %g\n", Kbin[i], Delta[flag][i], ShotLimit[flag][i],
		      Power[flag][i], (double) CountModes[flag][i], DeltaUncorrected[flag][i],
		      PowerUncorrected[flag][i], PowerSpec_Efstathiou(Kbin[i]), SumPower[flag][i],
		      4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3));
	    }
	}
      fclose(fd);
    }
}


#ifdef RAYLEIGH_BINS


struct pow_table
{
  double logk, logD;
}
 *PowerTable;

int compare_logk(const void *a, const void *b)
{
  if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk))
    return -1;

  if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk))
    return +1;

  return 0;
}

void prepare_rayleigh(void)
{
  int i;
  double kbin;
  FILE *fd;
  char buf[500];
  double k, p, dummy, ddd;
  int NPowerTable;
  double logk, logD, u, dlogk, Delta2;
  int binlow, binhigh, binmid;



  sprintf(buf, "inputspec_ics_millennium.txt");

  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
      endrun(17);
    }

  fscanf(fd, " %lg %lg ", &dummy, &dummy);


  NPowerTable = 0;
  do
    {
      if(fscanf(fd, " %lg %lg %lg %lg ", &k, &p, &dummy, &dummy) == 4)
	{
	  if(p>0)
	    NPowerTable++;
	}
      else
	break;
    }
  while(1);

  fclose(fd);


  if(ThisTask == 0)
    {
      printf("found %d pairs of values in input spectrum table\n", NPowerTable);
      fflush(stdout);
    }


  PowerTable = malloc(NPowerTable * sizeof(struct pow_table));


  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
      endrun(18);
    }

  fscanf(fd, " %lg %lg ", &dummy, &dummy);

  NPowerTable = 0;
  do
    {
      if(fscanf(fd, " %lg %lg %lg %lg ", &k, &p, &dummy, &dummy) == 4)
	{
	  if(p>0)
	    {
	      PowerTable[NPowerTable].logk = log10(k);
	      PowerTable[NPowerTable].logD = log10(p);
	      NPowerTable++;
	      
	      if(ThisTask==0)
		{
		  printf("k=%g  p=%g\n", k, p);
		  fflush(stdout);
		}
	    }
	}
      else
	break;
    }
  while(1);

  fclose(fd);

  qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);


  ddd = GrowthFactor(All.Time, 1.0);

  if(ThisTask == 0)
    {
      printf("growthfactor=%g\n", ddd);
      fflush(stdout);
    }


  for(i = 0; i < BINS_PS; i++)
    {
      kbin = exp((i + 0.5) / binfac + log(K0));

      RayleighFactor[i] = PowerSpec_Efstathiou(kbin);

      /* convert to h/Mpc */
      /*
         kbin *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);
       */


      logk = log10(kbin);

      if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
	{
	  if(ThisTask==0)
	    {
	      printf("Warning: logk=%g   %g %g\n", logk, PowerTable[0].logk, PowerTable[NPowerTable - 1].logk);
	    }
	  Delta2 = 1;
	}
      else
	{
	  binlow = 0;
	  binhigh = NPowerTable - 1;

	  while(binhigh - binlow > 1)
	    {
	      binmid = (binhigh + binlow) / 2;
	      if(logk < PowerTable[binmid].logk)
		binhigh = binmid;
	      else
		binlow = binmid;
	    }

	  dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

	  if(dlogk == 0)
	    endrun(777);

	  u = (logk - PowerTable[binlow].logk) / dlogk;

	  logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;

	  Delta2 = pow(10.0, logD);
	}

      Delta2 /= ddd * ddd;

      RayleighFactor[i] /= Delta2;

      RayleighFactor[i] *= 4 * M_PI * pow(kbin, 3) / pow(2 * M_PI / All.BoxSize, 3);
    }

  free(PowerTable);
}


double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}

double growth(double a)
{
#define WORKSIZE 100000
  double hubble_a, result, abserr;
  gsl_function F;
  gsl_integration_workspace *workspace;


  hubble_a = hubble_function(a);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  F.function = &growth_int;

  gsl_integration_qag(&F, 0, a, hubble_a,	/* note: absolute error just a dummy */
		      1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}

double growth_int(double a, void *param)
{
  return pow(a / (a * a * a * pow(hubble_function(a),2.0)), 1.5);
}

void rayleigh_save(void)
{
  FILE *fd;
  char buf[500];
  int i, rep;
  int *sum_RayleighCountModes;
  int *sum_RayleighCountAbove;
  float *max_RayleighMax;
  long long count_modes;


  sum_RayleighCountModes = mymalloc(BINS_PS * RAYLEIGH_BINS * sizeof(int));
  sum_RayleighCountAbove = mymalloc(BINS_PS * sizeof(int));
  max_RayleighMax = mymalloc(BINS_PS * sizeof(float));


  MPI_Allreduce(RayleighCountModes, sum_RayleighCountModes, BINS_PS * RAYLEIGH_BINS,
		MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(RayleighCountAbove, sum_RayleighCountAbove, BINS_PS, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(RayleighMax, max_RayleighMax, BINS_PS, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);


  if(ThisTask == 0)
    {
      for(rep = 0; rep < BINS_PS; rep++)
	{
	  for(i = 0, count_modes = sum_RayleighCountAbove[rep]; i < RAYLEIGH_BINS; i++)
	    count_modes += sum_RayleighCountModes[rep * RAYLEIGH_BINS + i];

	  if(count_modes)
	    {
	      sprintf(buf, "%s/power_rayleigh_%03d_%d.txt", All.OutputDir, All.PowerSpecFlag - 1, rep);

	      if(!(fd = fopen(buf, "w")))
		{
		  printf("can't open file `%s`\n", buf);
		  endrun(1326);
		}

	      fprintf(fd, "%g\n", All.Time);
	      fprintf(fd, "%g\n", max_RayleighMax[rep]);
	      fprintf(fd, "%d\n", (int) count_modes);
	      fprintf(fd, "%d\n", (int) sum_RayleighCountAbove[rep]);
	      fprintf(fd, "%d\n", (int) RAYLEIGH_BINS);
	      fprintf(fd, "%g\n", Kbin[rep]);

	      for(i = 0; i < RAYLEIGH_BINS; i++)
		fprintf(fd, "%d\n", sum_RayleighCountModes[rep * RAYLEIGH_BINS + i]);

	      fclose(fd);
	    }
	}
    }

  myfree(max_RayleighMax);
  myfree(sum_RayleighCountAbove);
  myfree(sum_RayleighCountModes);
}

#endif



void foldonitself(void)
{
  int i, j, level, sendTask, recvTask, istart, nbuf, n, rest, iter = 0;
  long int slab_x, slab_xx, slab_y, slab_yy, slab_z, slab_zz;
  int *nsend_local, *nsend_offset, *nsend, count, buf_capacity;
  double to_slab_fac_folded, dx, dy, dz;
  double tstart0, tstart, tend, t0, t1;
  FLOAT *pos_sendbuf, *pos_recvbuf, *pos;
  MPI_Status status;


  if(ThisTask == 0)
    {
      printf("begin folding for power spectrum estimation...\n");
      fflush(stdout);
    }

  tstart0 = tstart = second();

  nsend_local = mymalloc(NTask * sizeof(int));
  nsend_offset = mymalloc(NTask * sizeof(int));
  nsend = mymalloc(NTask * NTask * sizeof(int));

  pos_sendbuf = (FLOAT *) workspace;
  pos_recvbuf = (FLOAT *) forcegrid;

  buf_capacity = (dimprodmax * sizeof(fftw_real)) / (3 * sizeof(FLOAT));

  if(ThisTask == 0)
    {
      printf("buf_capacity = %d\n", buf_capacity);
      fflush(stdout);
    }


  to_slab_fac_folded = to_slab_fac * POWERSPEC_FOLDFAC;

  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;



  istart = 0;

  do
    {
      t0 = second();

      for(i = 0; i < NTask; i++)
	nsend_local[i] = 0;

      for(i = istart, nbuf = 0; i < NumPart; i++)
	{
	  if(nbuf + 1 >= buf_capacity)
	    break;

	  slab_x = to_slab_fac_folded * P[i].Pos[0];
	  slab_xx = slab_x + 1;
	  slab_x %= PMGRID;
	  slab_xx %= PMGRID;

	  nsend_local[slab_to_task[slab_x]]++;
	  nbuf++;

	  if(slab_to_task[slab_x] != slab_to_task[slab_xx])
	    {
	      nsend_local[slab_to_task[slab_xx]]++;
	      nbuf++;
	    }
	}

      for(i = 1, nsend_offset[0] = 0; i < NTask; i++)
	nsend_offset[i] = nsend_offset[i - 1] + nsend_local[i - 1];

      for(i = 0; i < NTask; i++)
	nsend_local[i] = 0;

      for(i = istart, nbuf = 0; i < NumPart; i++)
	{
	  if(nbuf + 1 >= buf_capacity)
	    break;

	  slab_x = to_slab_fac_folded * P[i].Pos[0];
	  slab_xx = slab_x + 1;
	  slab_x %= PMGRID;
	  slab_xx %= PMGRID;

	  for(j = 0; j < 3; j++)
	    pos_sendbuf[3 * (nsend_offset[slab_to_task[slab_x]] + nsend_local[slab_to_task[slab_x]]) + j] =
	      P[i].Pos[j];
	  nsend_local[slab_to_task[slab_x]]++;
	  nbuf++;

	  if(slab_to_task[slab_x] != slab_to_task[slab_xx])
	    {
	      for(j = 0; j < 3; j++)
		pos_sendbuf[3 * (nsend_offset[slab_to_task[slab_xx]] + nsend_local[slab_to_task[slab_xx]]) +
			    j] = P[i].Pos[j];
	      nsend_local[slab_to_task[slab_xx]]++;
	      nbuf++;
	    }
	}

      istart = i;


      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      t1 = second();
      if(ThisTask == 0)
	{
	  printf("buffer filled (took %g sec)\n", timediff(t0, t1));
	  fflush(stdout);
	}

      t0 = second();
      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(recvTask != sendTask)
		{
		  MPI_Sendrecv(&pos_sendbuf[3 * nsend_offset[recvTask]],
			       3 * nsend_local[recvTask] * sizeof(FLOAT), MPI_BYTE,
			       recvTask, TAG_PM_FOLD,
			       &pos_recvbuf[0],
			       3 * nsend[recvTask * NTask + ThisTask] * sizeof(FLOAT), MPI_BYTE,
			       recvTask, TAG_PM_FOLD, MPI_COMM_WORLD, &status);

		  pos = &pos_recvbuf[0];
		  count = nsend[recvTask * NTask + ThisTask];
		}
	      else
		{
		  pos = &pos_sendbuf[3 * nsend_offset[ThisTask]];
		  count = nsend_local[ThisTask];
		}

	      for(n = 0; n < count; n++, pos += 3)
		{
		  slab_x = to_slab_fac_folded * pos[0];
		  dx = to_slab_fac_folded * pos[0] - slab_x;
		  slab_xx = slab_x + 1;
		  slab_x %= PMGRID;
		  slab_xx %= PMGRID;

		  slab_y = to_slab_fac_folded * pos[1];
		  dy = to_slab_fac_folded * pos[1] - slab_y;
		  slab_yy = slab_y + 1;
		  slab_y %= PMGRID;
		  slab_yy %= PMGRID;

		  slab_z = to_slab_fac_folded * pos[2];
		  dz = to_slab_fac_folded * pos[2] - slab_z;
		  slab_zz = slab_z + 1;
		  slab_z %= PMGRID;
		  slab_zz %= PMGRID;

		  if(slab_to_task[slab_x] == ThisTask)
		    {
		      slab_x -= first_slab_of_task[ThisTask];

		      rhogrid[(slab_x * PMGRID + slab_y) * PMGRID2 + slab_z] +=
			(1.0 - dx) * (1.0 - dy) * (1.0 - dz);
		      rhogrid[(slab_x * PMGRID + slab_yy) * PMGRID2 + slab_z] += (1.0 - dx) * dy * (1.0 - dz);
		      rhogrid[(slab_x * PMGRID + slab_y) * PMGRID2 + slab_zz] += (1.0 - dx) * (1.0 - dy) * dz;
		      rhogrid[(slab_x * PMGRID + slab_yy) * PMGRID2 + slab_zz] += (1.0 - dx) * dy * dz;
		    }

		  if(slab_to_task[slab_xx] == ThisTask)
		    {
		      slab_xx -= first_slab_of_task[ThisTask];

		      rhogrid[(slab_xx * PMGRID + slab_y) * PMGRID2 + slab_z] +=
			(dx) * (1.0 - dy) * (1.0 - dz);
		      rhogrid[(slab_xx * PMGRID + slab_yy) * PMGRID2 + slab_z] += (dx) * dy * (1.0 - dz);
		      rhogrid[(slab_xx * PMGRID + slab_y) * PMGRID2 + slab_zz] += (dx) * (1.0 - dy) * dz;
		      rhogrid[(slab_xx * PMGRID + slab_yy) * PMGRID2 + slab_zz] += (dx) * dy * dz;
		    }
		}
	    }
	}

      count = NumPart - istart;	/* local remaining particles */
      MPI_Allreduce(&count, &rest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      iter++;

      t1 = second();
      if(ThisTask == 0)
	{
	  printf("particles exchanged and binned. (took %g sec) max-rest=%d\n", timediff(t0, t1), rest);
	  fflush(stdout);
	}
    }
  while(rest > 0);

  tend = second();

  if(ThisTask == 0)
    {
      printf("folded density field assembled (took %g seconds, iter=%d)\n", timediff(tstart, tend), iter);
      fflush(stdout);
    }

  tstart = second();

  /* Do the FFT of the self-folded density field */
#ifdef HPM
  hpmStart(4, "FFT");
#endif
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#ifdef HPM
  hpmStop(4);
#endif

  tend = second();

  if(ThisTask == 0)
    {
      printf("FFT for folded density done (took %g seconds)\n", timediff(tstart, tend));
      fflush(stdout);
    }

  myfree(nsend);
  myfree(nsend_offset);
  myfree(nsend_local);

  tend = second();
  All.CPU_PowerSpec += timediff(tstart0, tend);
}


#ifdef OUTPUT_LONGRANGE_POTENTIAL
void dump_potential(void)
{
  char buf[1000];
  int nprocgroup, masterTask, groupTask, n, i, j, k;
  double asmth, fac, box, tstart, tend;
  float *potential;
  FILE *fd;

  tstart = second();

  if(ThisTask == 0)
    {
      printf("Start dumping potential\n");
      fflush(stdout);
    }

  sprintf(buf, "%s/snapdir_%03d/potential_%03d.%d", All.OutputDir, All.PowerSpecFlag - 1,
	  All.PowerSpecFlag - 1, ThisTask);

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;

  masterTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(!(fd = fopen(buf, "w")))
	    {
	      printf("Error. Can't write in file '%s'\n", buf);
	      endrun(11);
	    }

	  n = PMGRID;
	  fwrite(&n, sizeof(int), 1, fd);

	  n = sizeof(float);
	  fwrite(&n, sizeof(int), 1, fd);

	  fwrite(&slabs_per_task[ThisTask], sizeof(int), 1, fd);
	  fwrite(&first_slab_of_task[ThisTask], sizeof(int), 1, fd);

	  box = All.BoxSize;
	  asmth = All.Asmth[0];

	  fwrite(&box, sizeof(double), 1, fd);
	  fwrite(&asmth, sizeof(double), 1, fd);

	  fac = All.G * All.PartMass / (M_PI * All.BoxSize);

	  potential = (float *) forcegrid;

	  for(i = 0; i < slabs_per_task[ThisTask]; i++)
	    for(j = 0; j < PMGRID; j++)
	      for(k = 0; k < PMGRID; k++)
		*potential++ = fac * rhogrid[(i * PMGRID + j) * PMGRID2 + k];

	  potential = (float *) forcegrid;

	  fwrite(potential, sizeof(float), PMGRID * PMGRID * slabs_per_task[ThisTask], fd);

	  fclose(fd);
	}

      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }


  MPI_Barrier(MPI_COMM_WORLD);

  tend = second();

  if(ThisTask == 0)
    {
      printf("finished writing potential (took=%g sec)\n", timediff(tstart, tend));
      fflush(stdout);
    }
}
#endif
