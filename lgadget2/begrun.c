#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"


/*! \file begrun.c
 *  \brief does the initial set-up of a simulation run
 */


/*!
 *  Does the initial set-up of the simulation.
 *  Reading the parameterfile, setting units, 
 *  getting IC's/restart files, etc.
 */
void begrun(void)
{
  int pmgrid_size;
  struct global_data_all_processes all;

  pmgrid_size = PMGRID;
  if(ThisTask == 0)
    {
      printf("\nThis is L-Gadget, version `%s'.\n", GADGETVERSION);
      printf("\nRunning on %d processors.\n", NTask);
      printf("\nCompiled with PMGRID = %d\n", pmgrid_size);
    }

#ifdef LIGHTCONE
  setup_lightcone_indexing();
#endif

  read_parameter_file(ParameterFile);	/* ... read in parameters for this run */

#ifndef SYSMALLOC
  mymalloc_init(All.MaxMemSize * 1024 * 1024);
#endif
  
  write_pid_file();

  set_units();

  open_outputfiles();

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  gsl_rng_set(random_generator, 42);	/* start-up seed */

  init_peano_map();

  All.TimeLastRestartFile = CPUThisRun;

  pm_init_periodic();

  if(RestartFlag == 0 || RestartFlag == 2)
    {
      init();			/* ... read in initial model */

      set_random_numbers();
    }
  else
    {
      all = All;		/* save global variables. (will be read from restart file) */

      restart(RestartFlag);	/* ... read restart file. Note: This also resets all variables in the struct
				   `All'.  However, during the run, some variables in the parameter file are
				   allowed to be changed, if desired. These need to copied in the way below.
				   Note: All.PartAllocFactor is treated in restart() separately.
				 */

      All.MinSizeTimestep = all.MinSizeTimestep;
      All.MaxSizeTimestep = all.MaxSizeTimestep;
      All.TreeAllocFactor = all.TreeAllocFactor;
      All.BufferSize = all.BufferSize;
      All.MaxMemSize = all.MaxMemSize;
      All.BunchSizeForce = all.BunchSizeForce;
      All.BunchSizeDomain = all.BunchSizeDomain;
      All.TimeLimitCPU = all.TimeLimitCPU;
      All.ResubmitOn = all.ResubmitOn;
      All.TimeBetSnapshot = all.TimeBetSnapshot;
      All.TimeBetStatistics = all.TimeBetStatistics;
      All.CpuTimeBetRestartFile = all.CpuTimeBetRestartFile;
      All.ErrTolIntAccuracy = all.ErrTolIntAccuracy;
      All.ErrTolForceAcc = all.ErrTolForceAcc;
      All.NumFilesWrittenInParallel = all.NumFilesWrittenInParallel;
      All.OutputListOn = all.OutputListOn;
      All.OutputListLength = all.OutputListLength;
#ifdef DARKENERGY
      All.DarkEnergyParam = all.DarkEnergyParam;
#endif

      memcpy(All.OutputListTimes, all.OutputListTimes, sizeof(double) * All.OutputListLength);

      strcpy(All.ResubmitCommand, all.ResubmitCommand);
      strcpy(All.OutputListFilename, all.OutputListFilename);
      strcpy(All.OutputDir, all.OutputDir);
      strcpy(All.RestartFile, all.RestartFile);
      strcpy(All.EnergyFile, all.EnergyFile);
      strcpy(All.InfoFile, all.InfoFile);
      strcpy(All.CpuFile, all.CpuFile);
      strcpy(All.TimingsFile, all.TimingsFile);
      strcpy(All.SnapshotFileBase, all.SnapshotFileBase);

      if(All.TimeMax != all.TimeMax)
	readjust_timebase(All.TimeMax, all.TimeMax);
    }



  set_softenings();

  force_treeinit();

  init_drift_table();

  if(RestartFlag == 2)
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current + 1);
  else
    All.Ti_nextoutput = find_next_outputtime(All.Ti_Current);

  All.TimeLastRestartFile = CPUThisRun;

#ifdef LIGHTCONE
  setup_lightcone();
  
  /*a couple of output diagnostics*/
  //float dorigin = sqrt(P[0].Pos[0]*P[0].Pos[0] + P[0].Pos[1]*P[0].Pos[1] + P[0].Pos[2]*P[0].Pos[2]);
  //printf("Position of first particle on task %d is %f, %f, %f -> z = %f, a = %f\n", ThisTask, 
  //P[0].Pos[0], P[0].Pos[1], P[0].Pos[2], ZofR(dorigin), 1./(1+ZofR(dorigin)));
#endif

#ifdef DARKENERGY
#ifdef TIMEDEPDE
  fwa_init();
#endif
#endif


}



void init_peano_map(void)
{
  int i, j, k, n;

  n = 0;

  for(i = 0; i < DOMAINGRID; i++)
    for(j = 0; j < DOMAINGRID; j++)
      for(k = 0; k < DOMAINGRID; k++)
	{
	  DomainPeanoMap[n] = peano_hilbert_key(i, j, k, DOMAINLEVELS);

	  n++;
	}
}

/*
 *  Compute conversion factors between internal code units
 *  and the cgs-system.  
 */
void set_units(void)
{
  All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
  All.UnitDensity_in_cgs = All.UnitMass_in_g / pow(All.UnitLength_in_cm, 3);
  All.UnitEnergy_in_cgs = All.UnitMass_in_g * pow(All.UnitLength_in_cm, 2) / pow(All.UnitTime_in_s, 2);

  /* convert some physical input parameters to internal units */

  All.Hubble = HUBBLE * All.UnitTime_in_s;

  if(All.GravityConstantInternal == 0)
    All.G = GRAVITY / pow(All.UnitLength_in_cm, 3) * All.UnitMass_in_g * pow(All.UnitTime_in_s, 2);
  else
    All.G = All.GravityConstantInternal;


  if(ThisTask == 0)
    {
      printf("\nHubble (internal units) = %g\n", All.Hubble);
      printf("G (internal units) = %g\n", All.G);
      printf("UnitMass_in_g = %g \n", All.UnitMass_in_g);
      printf("UnitTime_in_s = %g \n", All.UnitTime_in_s);
      printf("UnitVelocity_in_cm_per_s = %g \n", All.UnitVelocity_in_cm_per_s);
      printf("UnitDensity_in_cgs = %g \n", All.UnitDensity_in_cgs);
      printf("UnitEnergy_in_cgs = %g \n", All.UnitEnergy_in_cgs);
      printf("\n");
    }
}




/*
 *  Opens various log-files. On restart, the code 
 *  will append to the files.
 */
void open_outputfiles(void)
{
  char mode[2], buf[200];

  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;


  if(RestartFlag == 0)
    strcpy(mode, "w");
  else
    strcpy(mode, "a");


  sprintf(buf, "%s%s", All.OutputDir, All.CpuFile);
  if(!(FdCPU = fopen(buf, mode)))
  {
    printf("LGadget Error: 4\n");
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  }
  else
  {
    printf("LGadget Printing File: '%s'\n", buf);
    printf("outputcpu=%s\n", buf);
  }

  sprintf(buf, "%s%s", All.OutputDir, All.InfoFile);
  if(!(FdInfo = fopen(buf, mode)))
  {
    printf("LGadget Error: 5\n");
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  } 
  else
  {
    printf("LGadget Printing File: '%s'\n", buf);
    printf("outputinfo=%s\n", buf);
  }

  sprintf(buf, "%s%s", All.OutputDir, All.EnergyFile);
  if(!(FdEnergy = fopen(buf, mode)))
  {
    printf("LGadget Error: 6\n");
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  }
  else
  {
    printf("LGadget Printing File: '%s'\n", buf);
    printf("outputenergy=%s\n", buf);
  }

  sprintf(buf, "%s%s", All.OutputDir, All.TimingsFile);
  if(!(FdTimings = fopen(buf, mode)))
  {
    printf("LGadget Error: 7\n");
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  }
  else
  {
    printf("LGadget Printing File: '%s'\n", buf);
    printf("outputtimings=%s\n", buf);
  }

#ifdef DARKENERGY
  sprintf(buf, "%s%s", All.OutputDir, "darkenergy.txt");
  if(!(FdDE = fopen(buf, mode)))
  {
    printf("LGadget Error: 8\n");
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  }
  else
  {
    printf("LGadget Printing File: '%s'\n", buf);
    printf("outputdarkenergy=%s\n", buf);
    if(RestartFlag == 0)
    {
      fprintf(FdDE, "nstep time H(a) ");
#ifndef TIMEDEPDE
      fprintf(FdDE, "w0 Omega_L ");
#else
      fprintf(FdDE, "w(a) Omega_L ");
#endif
      fprintf(FdDE, "\n");
      fflush(FdDE);
    }
  }
#endif

#ifdef LIGHTCONE
  sprintf(buf, "%s%s", All.OutputDir, "lightcone.txt");
  if(!(FdLC = fopen(buf, mode)))
  {
    printf("LGadget Error: 9\n");
    printf("error in opening file '%s'\n", buf);
    endrun(1);
  }
  else
  {
    printf("LGadget Printing File: '%s'\n", buf);
    printf("outputlightcone=%s\n", buf);
  }
#endif

}


void close_outputfiles(void)
{
  if(ThisTask != 0)		/* only the root processors writes to the log files */
    return;

  fclose(FdCPU);
  fclose(FdInfo);
  fclose(FdEnergy);
  fclose(FdTimings);

#ifdef DARKENERGY
  fclose(FdDE);
#endif

#ifdef LICHTCONE
  fclose(FdLC);
#endif
  
}




/*
 *  This function parses the parameterfile in a simple way.
 *  Each paramater is defined by a keyword (`tag'), and can be
 *  either of type douple, int, or character string.
 *  The routine makes sure that each parameter appears 
 *  exactly once in the parameterfile.
 */
void read_parameter_file(char *fname)
{
#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd, *fdout;
  char buf[200], buf1[200], buf2[200], buf3[400];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int pnum, errorFlag = 0;

   if(sizeof(long long) != 8)
   {
     if(ThisTask == 0)
     {
       printf("LGadget Error: 10\n");
       printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
       endrun(0);
     }
   }
   if(sizeof(int) != 4)
   {
     if(ThisTask == 0)
     {
       printf("LGadget Error: 11\n");
       printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
       endrun(0);
     }
   }
   if(sizeof(float) != 4)
   {
     if(ThisTask == 0)
     {
       printf("LGadget Error: 12\n");
       printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
       endrun(0);
     }
   }
 
   if(sizeof(double) != 8)
   {
     if(ThisTask == 0)
     {
       printf("LGadget Error: 13\n");
       printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
       endrun(0);
     }
   }


  if(ThisTask == 0)		/* read parameter file on process 0 */
    {
      nt = 0;

      strcpy(tag[nt], "InitCondFile");
      addr[nt] = All.InitCondFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputDir");
      addr[nt] = All.OutputDir;
      id[nt++] = STRING;

      strcpy(tag[nt], "SnapshotFileBase");
      addr[nt] = All.SnapshotFileBase;
      id[nt++] = STRING;

      strcpy(tag[nt], "EnergyFile");
      addr[nt] = All.EnergyFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CpuFile");
      addr[nt] = All.CpuFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "InfoFile");
      addr[nt] = All.InfoFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "TimingsFile");
      addr[nt] = All.TimingsFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "RestartFile");
      addr[nt] = All.RestartFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "ResubmitCommand");
      addr[nt] = All.ResubmitCommand;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListFilename");
      addr[nt] = All.OutputListFilename;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputListOn");
      addr[nt] = &All.OutputListOn;
      id[nt++] = INT;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &All.Omega0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &All.OmegaBaryon;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &All.OmegaLambda;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "HubbleParam");
      addr[nt] = &All.HubbleParam;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &All.BoxSize;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeOfFirstSnapshot");
      addr[nt] = &All.TimeOfFirstSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "CpuTimeBetRestartFile");
      addr[nt] = &All.CpuTimeBetRestartFile;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetStatistics");
      addr[nt] = &All.TimeBetStatistics;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBegin");
      addr[nt] = &All.TimeBegin;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeMax");
      addr[nt] = &All.TimeMax;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TimeBetSnapshot");
      addr[nt] = &All.TimeBetSnapshot;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
      addr[nt] = &All.UnitVelocity_in_cm_per_s;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitLength_in_cm");
      addr[nt] = &All.UnitLength_in_cm;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "UnitMass_in_g");
      addr[nt] = &All.UnitMass_in_g;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolIntAccuracy");
      addr[nt] = &All.ErrTolIntAccuracy;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolTheta");
      addr[nt] = &All.ErrTolTheta;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ErrTolForceAcc");
      addr[nt] = &All.ErrTolForceAcc;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxSizeTimestep");
      addr[nt] = &All.MaxSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinSizeTimestep");
      addr[nt] = &All.MinSizeTimestep;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MaxRMSDisplacementFac");
      addr[nt] = &All.MaxRMSDisplacementFac;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "ICFormat");
      addr[nt] = &All.ICFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "SnapFormat");
      addr[nt] = &All.SnapFormat;
      id[nt++] = INT;

      strcpy(tag[nt], "NumFilesWrittenInParallel");
      addr[nt] = &All.NumFilesWrittenInParallel;
      id[nt++] = INT;

      strcpy(tag[nt], "ResubmitOn");
      addr[nt] = &All.ResubmitOn;
      id[nt++] = INT;

      strcpy(tag[nt], "TypeOfOpeningCriterion");
      addr[nt] = &All.TypeOfOpeningCriterion;
      id[nt++] = INT;

      strcpy(tag[nt], "TimeLimitCPU");
      addr[nt] = &All.TimeLimitCPU;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "Softening");
      addr[nt] = &All.Softening;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "SofteningMaxPhys");
      addr[nt] = &All.SofteningMaxPhys;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BufferSize");
      addr[nt] = &All.BufferSize;
      id[nt++] = INT;


      strcpy(tag[nt], "MaxMemSize");
      addr[nt] = &All.MaxMemSize;
      id[nt++] = INT;


      strcpy(tag[nt], "PartAllocFactor");
      addr[nt] = &All.PartAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TreeAllocFactor");
      addr[nt] = &All.TreeAllocFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "GravityConstantInternal");
      addr[nt] = &All.GravityConstantInternal;
      id[nt++] = DOUBLE;

#ifdef DARKENERGY
#ifdef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyFile");
      addr[nt] = All.DarkEnergyFile;
      id[nt++] = STRING;
#endif

#ifdef RESCALEVINI
      strcpy(tag[nt], "VelIniScale");
      addr[nt] = &All.VelIniScale;
      id[nt++] = DOUBLE;
#endif

#ifndef TIMEDEPDE
      strcpy(tag[nt], "DarkEnergyParam");
      addr[nt] = &All.DarkEnergyParam;
      id[nt++] = DOUBLE;
#endif
#endif



      if((fd = fopen(fname, "r")))
	{
	  sprintf(buf, "%s%s", fname, "-usedvalues");
	  if(!(fdout = fopen(buf, "w")))
	    {
	      printf("error opening file '%s' \n", buf);
	      errorFlag = 1;
	    }
	  else
	    {
	      while(!feof(fd))
		{
		  *buf = 0;
		  fgets(buf, 200, fd);
		  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
		    continue;

		  if(buf1[0] == '%')
		    continue;

		  for(i = 0, j = -1; i < nt; i++)
		    if(strcmp(buf1, tag[i]) == 0)
		      {
			j = i;
			tag[i][0] = 0;
			break;
		      }

		  if(j >= 0)
		    {
		      switch (id[j])
			{
			case DOUBLE:
			  *((double *) addr[j]) = atof(buf2);
			  fprintf(fdout, "%-35s%g\n", buf1, *((double *) addr[j]));
			  break;
			case STRING:
			  strcpy(addr[j], buf2);
			  fprintf(fdout, "%-35s%s\n", buf1, buf2);
			  break;
			case INT:
			  *((int *) addr[j]) = atoi(buf2);
			  fprintf(fdout, "%-35s%d\n", buf1, *((int *) addr[j]));
			  break;
			}
		    }
		  else
		    {
		      fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n",
			      fname, buf1);
		      errorFlag = 1;
		    }
		}
	      fclose(fd);
	      fclose(fdout);

	      i = strlen(All.OutputDir);
	      if(i > 0)
		if(All.OutputDir[i - 1] != '/')
		  strcat(All.OutputDir, "/");

	      sprintf(buf1, "%s%s", fname, "-usedvalues");
	      sprintf(buf2, "%s%s", All.OutputDir, "parameters-usedvalues");
	      sprintf(buf3, "cp %s %s", buf1, buf2);
	      system(buf3);
	    }
	}
      else
	{
	  printf("Parameter file %s not found.\n", fname);
	  errorFlag = 1;
	}


      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	      errorFlag = 1;
	    }
	}

      if(All.OutputListOn && errorFlag == 0)
	errorFlag += read_outputlist(All.OutputListFilename);
      else
	All.OutputListLength = 0;
    }

  MPI_Bcast(&errorFlag, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(errorFlag)
    {
      MPI_Finalize();
      exit(0);
    }

  /* now communicate the relevant parameters to the other processes */
  MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);



  for(pnum = 0; All.NumFilesWrittenInParallel > (1 << pnum); pnum++);

  if(All.NumFilesWrittenInParallel != (1 << pnum))
    {
      if(ThisTask == 0)
      {
        printf("LGadget Error: 14\n");
	printf("NumFilesWrittenInParallel MUST be a power of 2\n");
      }
      endrun(0);
    }

  if(All.NumFilesWrittenInParallel > NTask)
    {
      if(ThisTask == 0)
      {
        printf("LGadget Error: 15\n");
	printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
      }
      endrun(0);
    }

#undef DOUBLE
#undef STRING
#undef INT
#undef MAXTAGS

#ifdef TIMEDEPDE
#ifndef DARKENERGY
  if(ThisTask == 0)
  {
    printf("LGadget Error: 16\n");
    fprintf(stdout, "Code was compiled with TIMEDEPDE, but not with DARKENERGY.\n")     ;
    fprintf(stdout, "This is not allowed.\n");
  }
  endrun(0);
#endif
#endif

}


/* this function reads a table with a list of desired output
 * times. The table does not have to be ordered in any way,
 * but may not contain more than MAXLEN_OUTPUTLIST entries.
 */
int read_outputlist(char *fname)
{
  FILE *fd;

  if(!(fd = fopen(fname, "r")))
    {
      printf("can't read output list in file '%s'\n", fname);
      return 1;
    }

  All.OutputListLength = 0;
  do
    {
      if(fscanf(fd, " %lg ", &All.OutputListTimes[All.OutputListLength]) == 1)
	All.OutputListLength++;
      else
	break;
    }
  while(All.OutputListLength < MAXLEN_OUTPUTLIST);

  fclose(fd);

  printf("\nfound %d times in output-list.\n", All.OutputListLength);

  return 0;
}



void readjust_timebase(double TimeMax_old, double TimeMax_new)
{
  int i;
  long long ti_end;

  if(sizeof(long long) != 8)
    {
      if(ThisTask == 0)
      {
        printf("LGadget Error: 17\n");
	printf("\nType 'long long' is not 64 bit on this platform\n\n");
      }
      endrun(555);
    }

  if(ThisTask == 0)
    {
      printf("\nAll.TimeMax has been changed in the parameterfile\n");
      printf("Need to adjust integer timeline\n\n\n");
    }

  if(TimeMax_new < TimeMax_old)
    {
      if(ThisTask == 0)
      {
        printf("LGadget Error: 18\n");
	printf("\nIt is not allowed to reduce All.TimeMax\n\n");
      }
      endrun(556);
    }


  ti_end = log(TimeMax_new / All.TimeBegin) / All.Timebase_interval;

  while(ti_end > TIMEBASE)
    {
      All.Timebase_interval *= 2.0;

      ti_end /= 2;
      All.Ti_Current /= 2;
      All.PM_Ti_endstep /= 2;
      All.PM_Ti_begstep /= 2;

      for(i = 0; i < NumPart; i++)
	{
	  P[i].Ti_begstep /= 2;
	  P[i].Ti_endstep /= 2;
	}
    }

  All.TimeMax = TimeMax_new;
}
