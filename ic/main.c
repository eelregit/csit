#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <drfftw_mpi.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"


#define ASSERT_ALLOC(cond) {                                                                                  \
   if(cond)                                                                                                   \
    {                                                                                                         \
      if(ThisTask == 0)                                                                                       \
	printf("\nallocated %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                     \
    }                                                                                                         \
  else                                                                                                        \
    {                                                                                                         \
      printf("failed to allocate %g Mbyte on Task %d\n", bytes / (1024.0 * 1024.0), ThisTask);                \
      printf("bailing out.\n");                                                                               \
      FatalError(1);                                                                                          \
    }                                                                                                         \
}





int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  if(argc < 2)
    {
      if(ThisTask == 0)
	{
	  fprintf(stdout, "\nParameters are missing.\n");
	  fprintf(stdout, "Call with <ParameterFile>\n\n");
	}
      MPI_Finalize();
      exit(0);
    }

  read_parameterfile(argv[1]);

  set_units();

  initialize_powerspectrum();

  initialize_ffts();

  read_glass(GlassFile);

  displacement_fields();

  write_particle_data();

  if(NumPart)
    free(P);

  free_ffts();


  if(ThisTask == 0)
    {
      printf("\nIC's generated.\n\n");
      printf("Initial scale factor = %g\n", InitTime);
      printf("\n");
    }

  MPI_Barrier(MPI_COMM_WORLD);
  print_spec();

  MPI_Finalize();		/* clean up & finalize MPI */
  exit(0);
}





void displacement_fields(void)
{
  MPI_Request request;
  MPI_Status status;
  gsl_rng *random_generator;
  int i, j, k, ii, jj, kk, axes;
  int n;
  int sendTask, recvTask;
  double fac, vel_prefac, vel_prefac2;
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl, hubble_a;
  double tides = 0.; /* tidal correction to 1st displacement*/
  double tides2 = 0.; /* tidal correction to 2nd displacement */
  double trace = 0.; /* det = 1 + trace */
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double dis, dis2, maxdisp, max_disp_glob;
  double eps, eps2; /* 1st and 2nd tidal correction for displacement */
  unsigned int *seedtable;

  double lambda[3]; /*Tidal correction @z=initial */
  lambda[0] = Lambda_x/Dplus;
  lambda[1] = Lambda_y/Dplus;
  lambda[2] = Lambda_z/Dplus;
  if(ThisTask == 0)
  {
    printf("\nlambda[0] = %lf, lambda[1] = %lf, lambda[2] = %lf\n",lambda[0], lambda[1], lambda[2]);
  }

  unsigned long int bytes, nmesh3;
  int coord;
  fftw_complex *(cdisp[3]), *(cdisp2[3]) ; /* ZA and 2nd order displacements */
  fftw_real *(disp[3]), *(disp2[3]) ;

  fftw_complex *(cepsi[3]), *(cepsi2[3]) ; /* 1st and 2nd order epsilons */
  fftw_real *(epsi[3]), *(epsi2[3]) ;

  fftw_complex *(cdigrad[6]);
  fftw_real *(digrad[6]);

  fftw_complex *(cepgrad[6]); /* gradient of epsilon */
  fftw_real *(epgrad[6]);

#ifdef CORRECT_CIC
  double fx, fy, fz, ff, smth;
#endif

  if(ThisTask == 0)
    {
      printf("\nTidal Field @ z=0 : (%lf, %lf, %lf) \n", Lambda_x, Lambda_y, Lambda_z);
      printf("\nstart computing displacement fields...\n");
      fflush(stdout);
    }

  hubble_a =
    Hubble * sqrt(Omega / pow(InitTime, 3) + (1 - Omega - OmegaLambda) / pow(InitTime, 2) + OmegaLambda);

  vel_prefac = InitTime * hubble_a * F_Omega(InitTime);
  vel_prefac2 = InitTime * hubble_a * F2_Omega(InitTime);

  vel_prefac /= sqrt(InitTime);	/* converts to Gadget velocity */
  vel_prefac2 /= sqrt(InitTime);

  if(ThisTask == 0)
    printf("vel_prefac= %g, vel_prefac2= %g,  hubble_a=%g fom=%g \n", vel_prefac, vel_prefac2,
                                                                      hubble_a, F_Omega(InitTime));

  fac = pow(2 * PI / Box, 1.5);

  maxdisp = 0;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, Seed);

  if(!(seedtable = malloc(Nmesh * Nmesh * sizeof(unsigned int))))
    FatalError(4);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }



  for(axes=0,bytes=0; axes < 3; axes++)
    {
      cdisp[axes] = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      disp[axes] = (fftw_real *) cdisp[axes];
    }

  ASSERT_ALLOC(cdisp[0] && cdisp[1] && cdisp[2]);

  for(axes=0,bytes=0; axes < 3; axes++)
    {
      cepsi[axes] = (fftw_complex *) malloc(bytes += sizeof(fftw_real) * TotalSizePlusAdditional);
      epsi[axes] = (fftw_real *) cepsi[axes];
    }

  ASSERT_ALLOC(cepsi[0] && cepsi[1] && cepsi[2]);


#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
  for(Type = MinType; Type <= MaxType; Type++)
#endif
    {
      if(ThisTask == 0)
	{
	  printf("\nstarting axes=%d...\n", axes);
	  fflush(stdout);
	}

      /* first, clean the array */
      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    for(axes = 0; axes < 3; axes++)
	      {
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		cdisp[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
		cepsi[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].re = 0;
		cepsi[axes][(i * Nmesh + j) * (Nmesh / 2 + 1) + k].im = 0;
	      }

      for(i = 0; i < Nmesh; i++)
	{
	  ii = Nmesh - i;
	  if(ii == Nmesh)
	    ii = 0;
	  if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
	     (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
	    {
	      for(j = 0; j < Nmesh; j++)
		{
		  gsl_rng_set(random_generator, seedtable[i * Nmesh + j]);

		  for(k = 0; k < Nmesh / 2; k++)
		    {
		      phase = gsl_rng_uniform(random_generator) * 2 * PI;
		      do
			ampl = gsl_rng_uniform(random_generator);
		      while(ampl == 0);

		      if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
			continue;
		      if(i == 0 && j == 0 && k == 0)
			continue;

		      if(i < Nmesh / 2)
			kvec[0] = i * 2 * PI / Box;
		      else
			kvec[0] = -(Nmesh - i) * 2 * PI / Box;

		      if(j < Nmesh / 2)
			kvec[1] = j * 2 * PI / Box;
		      else
			kvec[1] = -(Nmesh - j) * 2 * PI / Box;

		      if(k < Nmesh / 2)
			kvec[2] = k * 2 * PI / Box;
		      else
			kvec[2] = -(Nmesh - k) * 2 * PI / Box;

		      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
		      kmag = sqrt(kmag2);

		      if(SphereMode == 1)
			{
			  if(kmag * Box / (2 * PI) > Nsample / 2)	/* select a sphere in k-space */
			    continue;
			}
		      else
			{
			  if(fabs(kvec[0]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[1]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			  if(fabs(kvec[2]) * Box / (2 * PI) > Nsample / 2)
			    continue;
			}

		      p_of_k = PowerSpec(kmag);

		      p_of_k *= -log(ampl);

		      delta = fac * sqrt(p_of_k) / Dplus;	/* scale back to starting redshift */

		      tides = (4./7.)*(kvec[0]*kvec[0]*lambda[0] + kvec[1]*kvec[1]*lambda[1] + kvec[2]*kvec[2]*lambda[2])/kmag2
		      			+3./7.*(lambda[0]+lambda[1]+lambda[2]); /* lambda[] are already scaled back.*/

		      if(k > 0)
			{
			  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
			    for(axes = 0; axes < 3; axes++)
			      {
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				  -kvec[axes] / kmag2 * delta * sin(phase);
				cepsi[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				  -kvec[axes] / kmag2 * delta * sin(phase) * tides;
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				  kvec[axes] / kmag2 * delta * cos(phase);
				cepsi[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				  kvec[axes] / kmag2 * delta * cos(phase) * tides;
			      }
			}
		      else	/* k=0 plane needs special treatment */
			{
			  if(i == 0)
			    {
			      if(j >= Nmesh / 2)
				continue;
			      else
				{
				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    {
				      jj = Nmesh - j;	/* note: j!=0 surely holds at this point */

				      for(axes = 0; axes < 3; axes++)
					{
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cepsi[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase) * tides;
					  cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase);
					  cepsi[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					    kvec[axes] / kmag2 * delta * cos(phase) * tides;

					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase);
					  cepsi[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
					    -kvec[axes] / kmag2 * delta * sin(phase) * tides;
					  cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					    -kvec[axes] / kmag2 * delta * cos(phase);
					  cepsi[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
					    -kvec[axes] / kmag2 * delta * cos(phase) * tides;
					}
				    }
				}
			    }
			  else	/* here comes i!=0 : conjugate can be on other processor! */
			    {
			      if(i >= Nmesh / 2)
				continue;
			      else
				{
				  ii = Nmesh - i;
				  if(ii == Nmesh)
				    ii = 0;
				  jj = Nmesh - j;
				  if(jj == Nmesh)
				    jj = 0;

				  if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					  -kvec[axes] / kmag2 * delta * sin(phase);
					cepsi[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
					  -kvec[axes] / kmag2 * delta * sin(phase) * tides;
					cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					  kvec[axes] / kmag2 * delta * cos(phase);
					cepsi[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
					  kvec[axes] / kmag2 * delta * cos(phase) * tides;
				      }

				  if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
				    for(axes = 0; axes < 3; axes++)
				      {
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].re = -kvec[axes] / kmag2 * delta * sin(phase);
					cepsi[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].re = -kvec[axes] / kmag2 * delta * sin(phase) * tides;
					cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].im = -kvec[axes] / kmag2 * delta * cos(phase);
					cepsi[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
					      k].im = -kvec[axes] / kmag2 * delta * cos(phase) * tides;
				      }
				}
			    }
			}
		    }
		}
	    }
	}



      /* At this point, cdisp[axes] contains the complex Zeldovich displacement */

       if(ThisTask == 0) printf("Done Zeldovich.\n");

      /* Compute displacement gradient */

      for(i = 0; i < 6; i++)
	{
	  cdigrad[i] = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  digrad[i] = (fftw_real *) cdigrad[i];
	  ASSERT_ALLOC(cdigrad[i]);
	}

      for(i = 0; i < 6; i++)
	{
	  cepgrad[i] = (fftw_complex *) malloc(bytes = sizeof(fftw_real) * TotalSizePlusAdditional);
	  epgrad[i] = (fftw_real *) cepgrad[i];
	  ASSERT_ALLOC(cepgrad[i]);
	}

      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;

	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;

	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;

	      /* Derivatives of ZA displacement  */
	      /* d(dis_i)/d(q_j)  -> sqrt(-1) k_j dis_i */
	      cdigrad[0][coord].re = -cdisp[0][coord].im * kvec[0]; /* disp0,0 */
	      cdigrad[0][coord].im = cdisp[0][coord].re * kvec[0];

	      cdigrad[1][coord].re = -cdisp[0][coord].im * kvec[1]; /* disp0,1 */
	      cdigrad[1][coord].im = cdisp[0][coord].re * kvec[1];

	      cdigrad[2][coord].re = -cdisp[0][coord].im * kvec[2]; /* disp0,2 */
	      cdigrad[2][coord].im = cdisp[0][coord].re * kvec[2];

	      cdigrad[3][coord].re = -cdisp[1][coord].im * kvec[1]; /* disp1,1 */
	      cdigrad[3][coord].im = cdisp[1][coord].re * kvec[1];

	      cdigrad[4][coord].re = -cdisp[1][coord].im * kvec[2]; /* disp1,2 */
	      cdigrad[4][coord].im = cdisp[1][coord].re * kvec[2];

	      cdigrad[5][coord].re = -cdisp[2][coord].im * kvec[2]; /* disp2,2 */
	      cdigrad[5][coord].im = cdisp[2][coord].re * kvec[2];

	      /* Derivatives of epsilon field  */
	      /* d(eps_i)/d(q_j) => sqrt(-1) k_j eps_i */
	      cepgrad[0][coord].re = -cepsi[0][coord].im * kvec[0]; /* eps0,0 */
	      cepgrad[0][coord].im = cepsi[0][coord].re * kvec[0];

	      cepgrad[1][coord].re = -cepsi[0][coord].im * kvec[1]; /* eps0,1 */
	      cepgrad[1][coord].im = cepsi[0][coord].re * kvec[1];

	      cepgrad[2][coord].re = -cepsi[0][coord].im * kvec[2]; /* eps0,2 */
	      cepgrad[2][coord].im = cepsi[0][coord].re * kvec[2];

	      cepgrad[3][coord].re = -cepsi[1][coord].im * kvec[1]; /* eps1,1 */
	      cepgrad[3][coord].im = cepsi[1][coord].re * kvec[1];

	      cepgrad[4][coord].re = -cepsi[1][coord].im * kvec[2]; /* eps1,2 */
	      cepgrad[4][coord].im = cepsi[1][coord].re * kvec[2];

	      cepgrad[5][coord].re = -cepsi[2][coord].im * kvec[2]; /* eps2,2 */
	      cepgrad[5][coord].im = cepsi[2][coord].re * kvec[2];
	    }


      if(ThisTask == 0) printf("Fourier transforming displacement gradient...");
      for(i = 0; i < 6; i++) rfftwnd_mpi(Inverse_plan, 1, digrad[i], Workspace, FFTW_NORMAL_ORDER);
      for(i = 0; i < 6; i++) rfftwnd_mpi(Inverse_plan, 1, epgrad[i], Workspace, FFTW_NORMAL_ORDER);
      if(ThisTask == 0) printf("Done.\n");

      /* Compute second order source and store it in digrad[3]*/
  	trace = lambda[0] + lambda[1] + lambda[2];


      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k < Nmesh; k++)
	    {
	      coord = (i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k;

	      epgrad[3][coord] =

	    -2./9.*lambda[0]*(digrad[0][coord]*digrad[0][coord] + digrad[1][coord]*digrad[1][coord] + digrad[2][coord]*digrad[2][coord] )
	    -2./9.*lambda[1]*(digrad[1][coord]*digrad[1][coord] + digrad[3][coord]*digrad[3][coord] + digrad[4][coord]*digrad[4][coord] )
	    -2./9.*lambda[2]*(digrad[2][coord]*digrad[2][coord] + digrad[4][coord]*digrad[4][coord] + digrad[5][coord]*digrad[5][coord] )
	    -trace/6.*(
	    	 digrad[0][coord]*digrad[0][coord] + digrad[3][coord]*digrad[3][coord] + digrad[5][coord]*digrad[5][coord]
	    	+digrad[1][coord]*digrad[1][coord] + digrad[2][coord]*digrad[2][coord] + digrad[4][coord]*digrad[4][coord]
	    	+digrad[0][coord]*digrad[3][coord] + digrad[3][coord]*digrad[5][coord] + digrad[5][coord]*digrad[0][coord] )
	    +1./18.*(digrad[0][coord]*(7.*epgrad[0][coord]-3.*epgrad[3][coord]-3.*epgrad[5][coord]) /* psi_xx, psi_yy, psi_zz part*/
	    		+digrad[3][coord]*(-3.*epgrad[0][coord]+7.*epgrad[3][coord]-3.*epgrad[5][coord])
	    		+digrad[5][coord]*(-3.*epgrad[0][coord]-3.*epgrad[3][coord]+7.*epgrad[5][coord]) )
	    +10./9.*(digrad[1][coord]*epgrad[1][coord] + digrad[2][coord]*epgrad[2][coord] + digrad[4][coord]*epgrad[4][coord]); /*psi_xy, psi_xz, psi_yz part*/

	      digrad[3][coord] =

		digrad[0][coord]*(digrad[3][coord]+digrad[5][coord])+digrad[3][coord]*digrad[5][coord]
                -digrad[1][coord]*digrad[1][coord]-digrad[2][coord]*digrad[2][coord]-digrad[4][coord]*digrad[4][coord];
	    }

      if(ThisTask == 0) printf("Fourier transforming second order source...");
      rfftwnd_mpi(Forward_plan, 1, digrad[3], Workspace, FFTW_NORMAL_ORDER);
      rfftwnd_mpi(Forward_plan, 1, epgrad[3], Workspace, FFTW_NORMAL_ORDER);
      if(ThisTask == 0) printf("Done.\n");

      /* The memory allocated for cdigrad[0], [1], and [2] will be used for 2nd order displacements */
      /* Freeing the rest. cdigrad[3] still has 2nd order displacement source, free later */

      for(axes = 0; axes < 3; axes++)
	{
	  cdisp2[axes] = cdigrad[axes];
	  disp2[axes] = (fftw_real *) cdisp2[axes];
	}

      free(cdigrad[4]); free(cdigrad[5]);

      for(axes = 0; axes < 3; axes++)
	{
	  cepsi2[axes] = cepgrad[axes];
	  epsi2[axes] = (fftw_real *) cepsi2[axes];
	}

      free(cepgrad[4]); free(cepgrad[5]);

      /* Solve Poisson eq. and calculate 2nd order displacements */

      for(i = 0; i < Local_nx; i++)
	for(j = 0; j < Nmesh; j++)
	  for(k = 0; k <= Nmesh / 2; k++)
	    {
	      coord = (i * Nmesh + j) * (Nmesh / 2 + 1) + k;
	      if((i + Local_x_start) < Nmesh / 2)
		kvec[0] = (i + Local_x_start) * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - (i + Local_x_start)) * 2 * PI / Box;

	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;

	      if(k < Nmesh / 2)
		kvec[2] = k * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - k) * 2 * PI / Box;

	      kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
#ifdef CORRECT_CIC
	      /* calculate smooth factor for deconvolution of CIC interpolation */
	      fx = fy = fz = 1;
	      if(kvec[0] != 0)
		{
		  fx = (kvec[0] * Box / 2) / Nmesh;
		  fx = sin(fx) / fx;
		}
	      if(kvec[1] != 0)
		{
		  fy = (kvec[1] * Box / 2) / Nmesh;
		  fy = sin(fy) / fy;
		}
	      if(kvec[2] != 0)
		{
		  fz = (kvec[2] * Box / 2) / Nmesh;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth = ff * ff;
	      /*  */
#endif

	      /* cdisp2 = source * k / (sqrt(-1) k^2) */
	    tides2 = (4./9.) * (kvec[0]*kvec[0]*lambda[0] + kvec[1]*kvec[1]*lambda[1] + kvec[2]*kvec[2]*lambda[2]) / kmag2;
	    trace = lambda[0] + lambda[1] + lambda[2];
	      for(axes = 0; axes < 3; axes++)
		{
		  if(kmag2 > 0.0)
		    {
		      cdisp2[axes][coord].re = cdigrad[3][coord].im * kvec[axes] / kmag2;
		      cdisp2[axes][coord].im = -cdigrad[3][coord].re * kvec[axes] / kmag2;
		      cepsi2[axes][coord].re = cepgrad[3][coord].im * kvec[axes] / kmag2
		      								- (3./7.)* (tides2 + trace/6. )* cdisp2[axes][coord].re;
		      cepsi2[axes][coord].im = -cepgrad[3][coord].re * kvec[axes] / kmag2
		      								- (3./7.)* (tides2 + trace/6. )* cdisp2[axes][coord].im;
		    }
		  else
		    {
		    	cdisp2[axes][coord].re = cdisp2[axes][coord].im = 0.0;
		  		cepsi2[axes][coord].re = cepsi2[axes][coord].im = 0.0;
		    }
#ifdef CORRECT_CIC
		  cdisp[axes][coord].re *= smth;   cdisp[axes][coord].im *= smth;
		  cdisp2[axes][coord].re *= smth;  cdisp2[axes][coord].im *= smth;
#endif
		}
	    }

      /* Free cdigrad[3] */
      free(cdigrad[3]);
      free(cepgrad[3]);


      /* Now, both cdisp, and cdisp2 have the ZA and 2nd order displacements */

      for(axes = 0; axes < 3; axes++)
	{
          if(ThisTask == 0) printf("Fourier transforming displacements, axis %d.\n",axes);

	  rfftwnd_mpi(Inverse_plan, 1, disp[axes], Workspace, FFTW_NORMAL_ORDER);
	  rfftwnd_mpi(Inverse_plan, 1, epsi[axes], Workspace, FFTW_NORMAL_ORDER);
	  rfftwnd_mpi(Inverse_plan, 1, disp2[axes], Workspace, FFTW_NORMAL_ORDER);
	  rfftwnd_mpi(Inverse_plan, 1, epsi2[axes], Workspace, FFTW_NORMAL_ORDER);

	  /* now get the plane on the right side from neighbour on the right,
	     and send the left plane */

	  recvTask = ThisTask;
	  do
	    {
	      recvTask--;
	      if(recvTask < 0)
		recvTask = NTask - 1;
	    }
	  while(Local_nx_table[recvTask] == 0);

	  sendTask = ThisTask;
	  do
	    {
	      sendTask++;
	      if(sendTask >= NTask)
		sendTask = 0;
	    }
	  while(Local_nx_table[sendTask] == 0);

	  /* use non-blocking send */

	  if(Local_nx > 0)
	    {
	      /* send ZA disp */
	      MPI_Isend(&(disp[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);

	      MPI_Recv(&(disp[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);

	      MPI_Wait(&request, &status);

	      /* send epsilon */
	      MPI_Isend(&(epsi[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);

	      MPI_Recv(&(epsi[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);

	      MPI_Wait(&request, &status);

	      /* send 2nd order disp */
	      MPI_Isend(&(disp2[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);

	      MPI_Recv(&(disp2[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);

	      MPI_Wait(&request, &status);

	      /* send 2nd order epsilon */
	      MPI_Isend(&(epsi2[axes][0]),
			sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
			MPI_BYTE, recvTask, 10, MPI_COMM_WORLD, &request);

	      MPI_Recv(&(epsi2[axes][(Local_nx * Nmesh) * (2 * (Nmesh / 2 + 1))]),
		       sizeof(fftw_real) * Nmesh * (2 * (Nmesh / 2 + 1)),
		       MPI_BYTE, sendTask, 10, MPI_COMM_WORLD, &status);

	      MPI_Wait(&request, &status);
	    }
	}

      /* read-out displacements */

	  nmesh3 = ((unsigned int ) Nmesh ) * ((unsigned int) Nmesh) *  ((unsigned int) Nmesh);

      for(n = 0; n < NumPart; n++)
	{
#if defined(MULTICOMPONENTGLASSFILE) && defined(DIFFERENT_TRANSFER_FUNC)
	  if(P[n].Type == Type)
#endif
	    {
	      u = P[n].Pos[0] / Box * Nmesh;
	      v = P[n].Pos[1] / Box * Nmesh;
	      w = P[n].Pos[2] / Box * Nmesh;

	      i = (int) u;
	      j = (int) v;
	      k = (int) w;

	      if(i == (Local_x_start + Local_nx))
		i = (Local_x_start + Local_nx) - 1;
	      if(i < Local_x_start)
		i = Local_x_start;
	      if(j == Nmesh)
		j = Nmesh - 1;
	      if(k == Nmesh)
		k = Nmesh - 1;

	      u -= i;
	      v -= j;
	      w -= k;

	      i -= Local_x_start;
	      ii = i + 1;
	      jj = j + 1;
	      kk = k + 1;

	      if(jj >= Nmesh)
		jj -= Nmesh;
	      if(kk >= Nmesh)
		kk -= Nmesh;

	      f1 = (1 - u) * (1 - v) * (1 - w);
	      f2 = (1 - u) * (1 - v) * (w);
	      f3 = (1 - u) * (v) * (1 - w);
	      f4 = (1 - u) * (v) * (w);
	      f5 = (u) * (1 - v) * (1 - w);
	      f6 = (u) * (1 - v) * (w);
	      f7 = (u) * (v) * (1 - w);
	      f8 = (u) * (v) * (w);

	      for(axes = 0; axes < 3; axes++)
		{
		  dis = disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

		  eps = epsi[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    epsi[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    epsi[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    epsi[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    epsi[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    epsi[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    epsi[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    epsi[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;

		  dis2 = disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    disp2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    disp2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    disp2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    disp2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;
		  dis2 /= (double) nmesh3;

		  eps2 = epsi2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f1 +
		    epsi2[axes][(i * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f2 +
		    epsi2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f3 +
		    epsi2[axes][(i * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f4 +
		    epsi2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + k] * f5 +
		    epsi2[axes][(ii * Nmesh + j) * (2 * (Nmesh / 2 + 1)) + kk] * f6 +
		    epsi2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + k] * f7 +
		    epsi2[axes][(ii * Nmesh + jj) * (2 * (Nmesh / 2 + 1)) + kk] * f8;
		  eps2 /= (double) nmesh3;

#ifdef ONLY_ZA
		  P[n].Pos[axes] += dis + eps;
		  P[n].Vel[axes] = dis * vel_prefac * sqrt(1.-lambda[axes]) + eps * vel_prefac2 * sqrt(1.-lambda[axes]);
#else
		  P[n].Pos[axes] += dis + eps - 3./7. * dis2 + eps2;
		  P[n].Vel[axes] = dis * vel_prefac * sqrt(1.-lambda[axes]) + eps * vel_prefac2 * sqrt(1.-lambda[axes])
		  					 + ( - 3./7. * dis2 + eps2 ) * vel_prefac2 * sqrt(1.-lambda[axes]) + eps2 * vel_prefac * sqrt(1.-lambda[axes]);
#endif

		  P[n].Pos[axes] = periodic_wrap(P[n].Pos[axes]);

		  if(dis + eps - 3./7. * dis2 + eps2 > maxdisp)
		    maxdisp = dis;
		}
	    }
	}
    }


  for(axes = 0; axes < 3; axes++) free(cdisp[axes]);
  for(axes = 0; axes < 3; axes++) free(cepsi[axes]);
  for(axes = 0; axes < 3; axes++) free(cdisp2[axes]);
  for(axes = 0; axes < 3; axes++) free(cepsi2[axes]);

  gsl_rng_free(random_generator);

  MPI_Reduce(&maxdisp, &max_disp_glob, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("\nMaximum displacement: %g kpc/h, in units of the part-spacing= %g\n",
	     max_disp_glob, max_disp_glob / (Box / Nmesh));
    }
}

double periodic_wrap(double x)
{
  while(x >= Box)
    x -= Box;

  while(x < 0)
    x += Box;

  return x;
}


void set_units(void)		/* ... set some units */
{
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;

  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  Hubble = HUBBLE * UnitTime_in_s;
}



void initialize_ffts(void)
{
  int total_size, i, additional;
  int local_ny_after_transpose, local_y_start_after_transpose;
  int *slab_to_task_local;
  size_t bytes;


  Inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);

  Forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD,
					 Nmesh, Nmesh, Nmesh, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

  rfftwnd_mpi_local_sizes(Forward_plan, &Local_nx, &Local_x_start,
			  &local_ny_after_transpose, &local_y_start_after_transpose, &total_size);

  Local_nx_table = malloc(sizeof(int) * NTask);
  MPI_Allgather(&Local_nx, 1, MPI_INT, Local_nx_table, 1, MPI_INT, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      for(i = 0; i < NTask; i++)
	printf("Task=%d Local_nx=%d\n", i, Local_nx_table[i]);
      fflush(stdout);
    }


  Slab_to_task = malloc(sizeof(int) * Nmesh);
  slab_to_task_local = malloc(sizeof(int) * Nmesh);

  for(i = 0; i < Nmesh; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < Local_nx; i++)
    slab_to_task_local[Local_x_start + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, Slab_to_task, Nmesh, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  free(slab_to_task_local);



  additional = (Nmesh) * (2 * (Nmesh / 2 + 1));	/* additional plane on the right side */

  TotalSizePlusAdditional = total_size + additional;

  //Disp = (fftw_real *) malloc(bytes = sizeof(fftw_real) * (total_size + additional));

  Workspace = (fftw_real *) malloc(bytes = sizeof(fftw_real) * total_size);

  ASSERT_ALLOC(Workspace)

  //Cdata = (fftw_complex *) Disp;	/* transformed array */
}



void free_ffts(void)
{
  free(Workspace);
  //free(Disp);
  free(Slab_to_task);
  rfftwnd_mpi_destroy_plan(Inverse_plan);
  rfftwnd_mpi_destroy_plan(Forward_plan);
}


int FatalError(int errnum)
{
  printf("FatalError called with number=%d\n", errnum);
  fflush(stdout);
  MPI_Abort(MPI_COMM_WORLD, errnum);
  exit(0);
}




static double A, B, alpha, beta, V, gf;

double fnl(double x)		/* Peacock & Dodds formula */
{
  return x * pow((1 + B * beta * x + pow(A * x, alpha * beta)) /
		 (1 + pow(pow(A * x, alpha) * gf * gf * gf / (V * sqrt(x)), beta)), 1 / beta);
}

void print_spec(void)
{
  double k, knl, po, dl, dnl, neff, kf, kstart, kend, po2, po1, DDD;
  char buf[1000];
  FILE *fd;

  if(ThisTask == 0)
    {
      sprintf(buf, "%s/inputspec_%s.txt", OutputDir, FileBase);

      fd = fopen(buf, "w");

      gf = GrowthFactor(0.001, 1.0) / (1.0 / 0.001);

      DDD = GrowthFactor(1.0 / (Redshift + 1), 1.0);

      fprintf(fd, "%12g %12g\n", Redshift, DDD);	/* print actual starting redshift and
							   linear growth factor for this cosmology */

      kstart = 2 * PI / (1000.0 * (3.085678e24 / UnitLength_in_cm));	/* 1000 Mpc/h */
      kend = 2 * PI / (0.001 * (3.085678e24 / UnitLength_in_cm));	/* 0.001 Mpc/h */

      for(k = kstart; k < kend; k *= 1.025)
	{
	  po = PowerSpec(k);
	  dl = 4.0 * PI * k * k * k * po;

	  kf = 0.5;

	  po2 = PowerSpec(1.001 * k * kf);
	  po1 = PowerSpec(k * kf);

	  if(po != 0 && po1 != 0 && po2 != 0)
	    {
	      neff = (log(po2) - log(po1)) / (log(1.001 * k * kf) - log(k * kf));

	      if(1 + neff / 3 > 0)
		{
		  A = 0.482 * pow(1 + neff / 3, -0.947);
		  B = 0.226 * pow(1 + neff / 3, -1.778);
		  alpha = 3.310 * pow(1 + neff / 3, -0.244);
		  beta = 0.862 * pow(1 + neff / 3, -0.287);
		  V = 11.55 * pow(1 + neff / 3, -0.423) * 1.2;

		  dnl = fnl(dl);
		  knl = k * pow(1 + dnl, 1.0 / 3);
		}
	      else
		{
		  dnl = 0;
		  knl = 0;
		}
	    }
	  else
	    {
	      dnl = 0;
	      knl = 0;
	    }

	  fprintf(fd, "%12g %12g    %12g %12g\n", k, dl, knl, dnl);
	}
      fclose(fd);
    }
}
