#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>

#include "allvars.h"
#include "proto.h"

#ifdef WALLCLOCK
#include <mpi.h>
#endif

/*! \file system.c
 *  \brief  contains a few low-level routines
 */



/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
#ifdef WALLCLOCK
  return MPI_Wtime();
#else
  return ((double) clock()) / CLOCKS_PER_SEC;
#endif

  /* note: on AIX and presumably many other 32bit systems, 
   * clock() has only a resolution of 10ms=0.01sec 
   */
}

/* returns the time difference between two measurements 
 * obtained with second(). The routine takes care of the 
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0, double t1)
{
  double dt;

  dt = t1 - t0;

  if(dt < 0)			/* overflow has occured (for systems with 32bit tick counter) */
    {
#ifdef WALLCLOCK
      dt = 0;
#else
      dt = t1 + pow(2, 32) / CLOCKS_PER_SEC - t0;
#endif
    }

  return dt;
}

double get_random_number(int id)
{
  return RndTable[(id % RNDTABLE)];
}

void set_random_numbers(void)
{
  int i;

  for(i = 0; i < RNDTABLE; i++)
    RndTable[i] = gsl_rng_uniform(random_generator);
}



/* returns the maximum of two double
 */
double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}

/* returns the minimum of two double
 */
double dmin(double x, double y)
{
  if(x < y)
    return x;
  else
    return y;
}



/* returns the maximum of two integers
 */
int imax(int x, int y)
{
  if(x > y)
    return x;
  else
    return y;
}

/* returns the minimum of two integers
 */
int imin(int x, int y)
{
  if(x < y)
    return x;
  else
    return y;
}

#ifdef DEBUG
#include <fenv.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <signal.h>
void enable_core_dumps_and_fpu_exceptions(void)
{
  struct rlimit rlim;
  extern int feenableexcept(int __excepts);

  /* enable floating point exceptions */
  feenableexcept(FE_DIVBYZERO | FE_INVALID);

  /* Note: FPU exceptions appear not to work properly 
   * when the Intel C-Compiler for Linux is used
   */

  /* set core-dump size to infinity */
  getrlimit(RLIMIT_CORE, &rlim);
  rlim.rlim_cur = RLIM_INFINITY;
  setrlimit(RLIMIT_CORE, &rlim);

  /* MPICH catches the signales SIGSEGV, SIGBUS, and SIGFPE....
   * The following statements reset things to the default handlers,
   * which will generate a core file.  
   */
  signal(SIGSEGV, SIG_DFL);
  signal(SIGBUS, SIG_DFL);
  signal(SIGFPE, SIG_DFL);
  signal(SIGINT, SIG_DFL);

}
#endif

void write_pid_file(void)
{
  pid_t my_pid;
  char mode[8], buf[500];
  FILE *fd;
  int i;

  my_pid = getpid();

  sprintf(buf, "%s%s", All.OutputDir, "PIDs.txt");

  for(i = 0; i < NTask; i++)
    {
      if(ThisTask == i)
        {
          if(ThisTask == 0)
            sprintf(mode, "w");
          else
            sprintf(mode, "a");

          if((fd = fopen(buf, mode)))
            {
              fprintf(fd, "%s %d\n", getenv("HOST"), (int) my_pid);
              fclose(fd);
            }
        }

      MPI_Barrier(MPI_COMM_WORLD);
    }
}
