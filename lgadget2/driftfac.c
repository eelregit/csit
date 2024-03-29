#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"


/*! \file driftfac.c
 *  \brief tabulates the cosmological drift and kick factors
 */

static double logTimeBegin; 
static double logTimeMax;

void init_drift_table(void)
{
#define WORKSIZE 100000
  int i,j;
  double result, abserr;
  gsl_function F;
  gsl_integration_workspace  *workspace;

  logTimeBegin= log(All.TimeBegin);
  logTimeMax= log(All.TimeMax);

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  for(i = 0; i < DRIFT_TABLE_LENGTH; i++)
    {
      for(j = 0; j < 3; j++){
        F.function = &drift_integ;
        F.params = &j;
        gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,  /* note: absolute error just a dummy */
          1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
        DriftTable[j][i] = result;
      }

      F.function = &gravkick_integ;
      gsl_integration_qag(&F, exp(logTimeBegin), exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1)), All.Hubble,	/* note: absolute error just a dummy */
			  1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
      GravKickTable[i] = result;

    }
  
  gsl_integration_workspace_free(workspace);
}



/*
double drift_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * a * a * a);
}
*/
double drift_integ(double a, void *param)
{
  double h;
  double anifac[3];
  int *axes = (int *)param;

  h = hubble_function(a);
  eval_aniss(a, anifac);

  return 1 / (h * a * anifac[*axes] * anifac[*axes]);
}


double gravkick_integ(double a, void *param)
{
  double h;

  h = hubble_function(a);

  return 1 / (h * a * a);
}

double growthfactor_integ(double a, void *param)
{
  double s;
 
  s = hubble_function(a) / All.Hubble * sqrt(a * a * a);

  return pow(sqrt(a) / s, 3);
}





double get_drift_factor(int time0, int time1, int axes)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * DriftTable[axes][0];
  else
    df1 = DriftTable[axes][i1 - 1] + (DriftTable[axes][i1] - DriftTable[axes][i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * DriftTable[axes][0];
  else
    df2 = DriftTable[axes][i2 - 1] + (DriftTable[axes][i2] - DriftTable[axes][i2 - 1]) * (u2 - i2);

  return df2 - df1;
}


double get_gravkick_factor(int time0, int time1)
{
  double a1, a2, df1, df2, u1, u2;
  int i1, i2;

  /* note: will only be called for cosmological integration */

  a1 = logTimeBegin + time0 * All.Timebase_interval;
  a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * GravKickTable[0];
  else
    df1 = GravKickTable[i1 - 1] + (GravKickTable[i1] - GravKickTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * GravKickTable[0];
  else
    df2 = GravKickTable[i2 - 1] + (GravKickTable[i2] - GravKickTable[i2 - 1]) * (u2 - i2);

  return df2 - df1;
}
