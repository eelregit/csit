#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

/* This routine checks for when a particles redshift (from the origin)
   is less than the particles time and write out the particle to a
   lightcone.  I would like to modify this to allow for the lightcone
   center to be placed arbitrarily (i.e. center for a full sky lightcone) */

#ifdef LIGHTCONE

static double logTimeBegin; 
static double logTimeMax;
static double a_DriftTable[DRIFT_TABLE_LENGTH];
static double glX[LC_NGAULEG_STEPS];
static double glW[LC_NGAULEG_STEPS];

static void gauleg_coeffs(int N, double *x, double *w);
static double RofZGL(double Z);

void setup_lightcone_indexing(void)
{
  int i;
  
  /*set the number corresponding to which lightcones we're making*/
  for(i=0;i<NLCS;i++)
    LCNum[i] = 0;
  
#ifdef OCTANT1
  LCNum[0] = 1;
#endif
#ifdef OCTANT2
  LCNum[1] = 1;
#endif
#ifdef OCTANT3
  LCNum[2] = 1;
#endif
#ifdef OCTANT4
  LCNum[3] = 1;
#endif
#ifdef OCTANT5
  LCNum[4] = 1;
#endif
#ifdef OCTANT6
  LCNum[5] = 1;
#endif
#ifdef OCTANT7
  LCNum[6] = 1;
#endif
#ifdef OCTANT8
  LCNum[7] = 1;
#endif
#ifdef SPHERE
  LCNum[8] = 1;
#endif
  
  if(ThisTask == 0)
    {
      printf("\nmaking light cones %d %d %d %d %d %d %d %d %d.\n",
	     LCNum[0], LCNum[1], LCNum[2], LCNum[3], LCNum[4],
	     LCNum[5], LCNum[6], LCNum[7], LCNum[8]);
#ifdef LIGHTCONE_OUTPUT_AINT
      printf("\nfinalized light cone files will record particle scale factors.\n");
#endif
    }
}  

/*the initial setup to make sure all the files are present, etc...*/
void setup_lightcone(void)
{
  int i;
  char buf[1000];
  
  /*initialize the output buffer to zero*/
  for(i=0;i<NLCS;i++)
    NPart_in_buf[i] = 0;

  //make lightcone and lightcone_buff directories
  if(ThisTask == 0)
    {
      sprintf(buf, "%slightcone", All.OutputDir);
      mkdir(buf, 02755);
      
      for(i=0;i<NLCS;i++)
	{
	  if(LCNum[i] == 0)
	    continue;
	  
	  sprintf(buf, "%slightcone/lightcone%03d", All.OutputDir, i);
	  mkdir(buf, 02755);
	}
      
      sprintf(buf, "%slightcone_buff", All.OutputDir);
      mkdir(buf, 02755);
      
      for(i=0;i<NLCS;i++)
	{
	  if(LCNum[i] == 0)
	    continue;
	  
	  sprintf(buf, "%slightcone_buff/lightcone%03d", All.OutputDir, i);
	  mkdir(buf, 02755);
	  
	  sprintf(buf, "%slightcone_buff/lightcone%03d/p", All.OutputDir, i);
	  mkdir(buf, 02755);
	  
	  sprintf(buf, "%slightcone_buff/lightcone%03d/v", All.OutputDir, i);
	  mkdir(buf, 02755);
	  
	  sprintf(buf, "%slightcone_buff/lightcone%03d/a", All.OutputDir, i);
	  mkdir(buf, 02755);
	  
	  sprintf(buf, "%slightcone_buff/lightcone%03d/i", All.OutputDir, i);
	  mkdir(buf, 02755);
	}
    }
  
  /* put a barrier here just in case */
  ////////////////////////////////////
  MPI_Barrier(MPI_COMM_WORLD);
  ////////////////////////////////////
  
  /*Remove old output files if they exist and we're not restarting*/
  char fname_p[500], fname_v[500], fname_i[500], fname_a[500];
  FILE *fp;
  int nprocgroup, groupTask, masterTask;
  if(RestartFlag == 0) //this is a fresh run so init everything
    {
      nprocgroup = NTask / All.NumFilesWrittenInParallel;
      if((NTask % All.NumFilesWrittenInParallel))
	nprocgroup++;
      masterTask = (ThisTask / nprocgroup) * nprocgroup;

      for(groupTask = 0; groupTask < nprocgroup; groupTask++)
	{
	  if(ThisTask == (masterTask + groupTask))
	    {
	      for(i=0;i<NLCS;i++)
		{
		  if(LCNum[i] == 0)
		    continue;
		  
		  //remove old files
		  if(ThisTask == 0) 
		    printf("Removing any old Lightcone files...\n");
		  sprintf(fname_p, "%slightcone_buff/lightcone%03d/p/%s_Lightcone.p.%d.%d", All.OutputDir, i, All.SnapshotFileBase, ThisTask, i);
		  sprintf(fname_v, "%slightcone_buff/lightcone%03d/v/%s_Lightcone.v.%d.%d", All.OutputDir, i, All.SnapshotFileBase, ThisTask, i);
		  sprintf(fname_i, "%slightcone_buff/lightcone%03d/i/%s_Lightcone.i.%d.%d", All.OutputDir, i, All.SnapshotFileBase, ThisTask, i);
		  sprintf(fname_a, "%slightcone_buff/lightcone%03d/a/%s_Lightcone.a.%d.%d", All.OutputDir, i, All.SnapshotFileBase, ThisTask, i);
		  remove(fname_p);
		  remove(fname_v);
		  remove(fname_i);
		  remove(fname_a);
			      
		  //now touch the files as blank so that the finalize routine won't fail
		  fp = fopen(fname_p,"a");
		  fclose(fp);
		  fp = fopen(fname_v,"a");
		  fclose(fp);
		  fp = fopen(fname_i,"a");
		  fclose(fp);
		  fp = fopen(fname_a,"a");
		  fclose(fp);
		}
	    }
	  
	  /* wait inside the group */
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      /*we also want to clear the output data if we're not restarting*/
      //for(i=0;i<All.MaxPart;i++)
      //P[i].LC_out = 0;
      for(i=0;i<NLCS;i++)
	NPart_Written[i] = 0;
    }
  
  //drift table factors
  logTimeBegin = log(All.TimeBegin);
  logTimeMax = log(All.TimeMax);
  
  //init gaus leg quad weights
  gauleg_coeffs(LC_NGAULEG_STEPS,glX,glW);
  
  /*generate the lookup table for calculating redshift of a distance*/
  MakeZofRTable();

  /*The largest time needed of the box*/
  double rbox = sqrt(3*All.BoxSize*All.BoxSize);
  double zbox = ZofR(rbox);
  double abox = 1./(1+zbox);
  if (ThisTask == 0) 
    printf("On the opposite corner of the box:  r = %f, z = %f, a = %f\n", rbox, zbox, abox);
  
  //reset (or init really) profile info for all tasks here
  LCProf.NumLCTests = 0;
  LCProf.NumLCFound = 0;
  LCProf.CPU_CheckLC = 0.0;
  LCProf.CPU_IO = 0.0;
  
  //init interp for dt to alpha factor
  for(i=0;i<DRIFT_TABLE_LENGTH;++i)
    a_DriftTable[i] = exp(logTimeBegin + ((logTimeMax - logTimeBegin) / DRIFT_TABLE_LENGTH) * (i + 1));
}

/*
  Gauss Leg. integration routines
  from SDSSIDL package by E. Sheldon (C) 2005
  http://code.google.com/p/sdssidl/ 
  lifted under GNU GPL v2
*/
static double abszdiff(double zdiff) 
{
  if (zdiff < 0.0) 
    return(-1.*zdiff);
  else
    return(zdiff);
}

static void gauleg_coeffs(int N, double *x, double *w)
{
  int i, j, m;
  double xm, xl, z1, z, p1, p2, p3, pp=0, pi, EPS;
  double x1 = 0.0;
  double x2 = 1.0;
  
  EPS = 3.e-11;
  pi = M_PI;

  m = (N + 1)/2;

  xm = (x1 + x2)/2.0;
  xl = (x2 - x1)/2.0;
  z1 = 0.0;

  for (i=1; i<= m; ++i) 
    {

      z=cos( pi*(i-0.25)/(N+.5) );

      while (abszdiff(z-z1) > EPS) 
        {
          p1 = 1.0;
          p2 = 0.0;
          for (j=1; j <= N;++j)
            {
              p3 = p2;
              p2 = p1;
              p1 = ( (2.0*j - 1.0)*z*p2 - (j-1.0)*p3 )/j;
            }
          pp = N*(z*p1 - p2)/(z*z -1.);
          z1=z;
          z=z1 - p1/pp;

        }

      x[i-1] = xm - xl*z;
      x[N+1-i-1] = xm + xl*z;
      w[i-1] = 2.0*xl/( (1.-z*z)*pp*pp );
      w[N+1-i-1] = w[i-1];

    }
}

static double get_drift_factor_lc(double a1, double a2)
{
  double df1, df2, u1, u2;
  int i1, i2;
  
  /* note: will only be called for cosmological integration */
  
  //a1 = logTimeBegin + time0 * All.Timebase_interval;
  //a2 = logTimeBegin + time1 * All.Timebase_interval;

  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * DriftTable[0];
  else
    df1 = DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);


  u2 = (a2 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i2 = (int) u2;
  if(i2 >= DRIFT_TABLE_LENGTH)
    i2 = DRIFT_TABLE_LENGTH - 1;

  if(i2 <= 1)
    df2 = u2 * DriftTable[0];
  else
    df2 = DriftTable[i2 - 1] + (DriftTable[i2] - DriftTable[i2 - 1]) * (u2 - i2);
  
  return df2 - df1;
}

static double get_drift_factor_lc_absval(double a1)
{
  double df1, u1;
  int i1;
  
  /* note: will only be called for cosmological integration */
  
  //a1 = logTimeBegin + time0 * All.Timebase_interval;
  
  u1 = (a1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;

  if(i1 <= 1)
    df1 = u1 * DriftTable[0];
  else
    df1 = DriftTable[i1 - 1] + (DriftTable[i1] - DriftTable[i1 - 1]) * (u1 - i1);
  
  return df1;
}

static double AofdT(double dt, int i0, int i1)
{
  int i;
  double slope = 0.0;
  
  for(i=i0;i<=i1;++i)
    {
      if(dt < DriftTable[i])
	break;
    }
  slope = (a_DriftTable[i]-a_DriftTable[i-1])/(DriftTable[i]-DriftTable[i-1]);
  
  return a_DriftTable[i-1] + slope*(dt - DriftTable[i-1]);
}

static void putPartInLC(float LCPos[3], int part, float a, int LC)
{
  int k;
  double t0,t1;
  
  /*transfer information to output buffer*/
  ++(LCProf.NumLCFound);
  for(k=0;k<3;k++)
    {
      out_particles[LC][NPart_in_buf[LC]].Pos[k] = LCPos[k];
      out_particles[LC][NPart_in_buf[LC]].Vel[k] = P[part].Vel[k]/a;
    }
  out_particles[LC][NPart_in_buf[LC]].ID = P[part].ID;
  out_particles[LC][NPart_in_buf[LC]].a = a;
  NPart_in_buf[LC]++;
  if(NPart_in_buf[LC] == LC_OUTPUT_BUF)
    {
      //here we track the buffer I/O time separately 
      //we subtract it from the particle checking time since that profiling code contains this one
      t0 = second();
      
      dump_lc_buf(LC);
      
      t1 = second();
      LCProf.CPU_IO += timediff(t0,t1);
      LCProf.CPU_CheckLC -= timediff(t0,t1);
    }
  
  //NOT doing this anymore...
  /*update markers to indicate that particle has been written*/
  //P[part].LC_out += 1<<LC;
}

static void checkPutPartInLC(float LCPos[3], int part, float a)
{
  
  //to make sure we do not get duplicates on edges some how
  //always use >= 0.0 and < 0.0
  
#ifdef OCTANT1
  //if(P[part].LC_out % 2 == 0)
  if(LCPos[0] >= 0.0 && LCPos[1] >= 0.0 && LCPos[2] >= 0.0)
    putPartInLC(LCPos,part,a,0);
#endif
#ifdef OCTANT2
  //if ((P[part].LC_out >> 1) % 2 == 0)
  if(LCPos[0] < 0.0 && LCPos[1] >= 0.0 && LCPos[2] >= 0.0)
    putPartInLC(LCPos,part,a,1);
#endif
#ifdef OCTANT3
  //if ((P[part].LC_out >> 2) % 2 == 0)
  if(LCPos[0] >= 0.0 && LCPos[1] < 0.0 && LCPos[2] >= 0.0)
    putPartInLC(LCPos,part,a,2);
#endif
#ifdef OCTANT4
  //if ((P[part].LC_out >> 3) % 2 == 0)
  if(LCPos[0] >= 0.0 && LCPos[1] >= 0.0 && LCPos[2] < 0.0)
    putPartInLC(LCPos,part,a,3);
#endif
#ifdef OCTANT5
  //if ((P[part].LC_out >> 4) % 2 == 0)
  if(LCPos[0] < 0.0 && LCPos[1] < 0.0 && LCPos[2] >= 0.0)
    putPartInLC(LCPos,part,a,4);
#endif
#ifdef OCTANT6
  //if ((P[part].LC_out >> 5) % 2 == 0)
  if(LCPos[0] < 0.0 && LCPos[1] >= 0.0 && LCPos[2] < 0.0)
    putPartInLC(LCPos,part,a,5);
#endif
#ifdef OCTANT7
  //if ((P[part].LC_out >> 6) % 2 == 0)
  if(LCPos[0] >= 0.0 && LCPos[1] < 0.0 && LCPos[2] < 0.0)
    putPartInLC(LCPos,part,a,6);
#endif
#ifdef OCTANT8
  //if ((P[part].LC_out >> 7) % 2 == 0)
  if(LCPos[0] < 0.0 && LCPos[1] < 0.0 && LCPos[2] < 0.0)
    putPartInLC(LCPos,part,a,7);
#endif

#ifdef SPHERE
  //if ((P[part].LC_out >> 8) % 2 == 0)
  putPartInLC(LCPos,part,a,8);
#endif
}

void init_lc_interp(int time0, int time1)
{
  int i;
  double dloga,a,loga,R2[LC_NQUADFIT],alpha[LC_NQUADFIT];
  int nadd;
  double chisq,b,l1err;
  int i0,i1;
  double u0,u1;
  
  //by default will try quadratic interp
  LCInterp.doLinearInterp = 0;
  
  //init range of table
  LCInterp.a0 = All.TimeBegin * exp(time0*All.Timebase_interval);
  LCInterp.loga0 = log(LCInterp.a0);
  LCInterp.R2Prev = RofZGL(1.0/LCInterp.a0 - 1.0);
  LCInterp.RPrevRCurr = LCInterp.R2Prev;
  LCInterp.R2Prev = LCInterp.R2Prev*LCInterp.R2Prev;
  
  LCInterp.a1 = All.TimeBegin * exp(time1*All.Timebase_interval);
  LCInterp.loga1 = log(LCInterp.a1);
  LCInterp.R2Curr = RofZGL(1.0/LCInterp.a1 - 1.0);
  LCInterp.RPrevRCurr *= LCInterp.R2Curr;
  LCInterp.R2Curr = LCInterp.R2Curr*LCInterp.R2Curr;
  
  LCInterp.dt_drift = get_drift_factor_lc(LCInterp.loga0,LCInterp.loga1);
  LCInterp.dt_a0 = get_drift_factor_lc_absval(LCInterp.loga0);
  
  u0 = (LCInterp.loga0 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i0 = (int) u0;
  i0 = i0-3;
  if(i0 <= 0)
    i0 = 1;
  LCInterp.i0 = i0;
  
  u1 = (LCInterp.loga1 - logTimeBegin) / (logTimeMax - logTimeBegin) * DRIFT_TABLE_LENGTH;
  i1 = (int) u1;
  i1 = i1+1;
  if(i1 >= DRIFT_TABLE_LENGTH)
    i1 = DRIFT_TABLE_LENGTH - 1;
  LCInterp.i1 = i1;

  //get quadratic fit to R2(alpha)
  dloga = (LCInterp.loga1 - LCInterp.loga0)/(LC_NQUADFIT-3);
  nadd = 0;
  b = 0.0;
  for(i=0;i<LC_NQUADFIT;++i)
    {
      loga = LCInterp.loga0 + dloga*(i-1);
      a = exp(loga);
      R2[i] = RofZGL(1.0/a-1.0);
      R2[i] = R2[i]*R2[i];
      alpha[i] = get_drift_factor_lc(LCInterp.loga0,loga)/LCInterp.dt_drift;
      
      if(alpha[i] != 0 && alpha[i] != 1.0 && gsl_finite(alpha[i]))
	{
	  b += (R2[i] - LCInterp.R2Prev - (LCInterp.R2Curr - LCInterp.R2Prev)*alpha[i]*alpha[i])
	    /(alpha[i] - alpha[i]*alpha[i]);
	  ++nadd;
	}
    }
  
  if(nadd == 0)
    {
      if(ThisTask == 0)
	printf("doing linear LC interp: could not get quadratic fit for LC R2(alpha)\n!");
      LCInterp.doLinearInterp = 1;
      return;
    }
  
  b /= nadd;
  LCInterp.p[0] = LCInterp.R2Curr - LCInterp.R2Prev - b;
  LCInterp.p[1] = b;
  LCInterp.p[2] = LCInterp.R2Prev;
  
  chisq = 0.0;
  l1err = 0.0;
  for(i=0;i<LC_NQUADFIT;++i)
    {
      b = sqrt(LCInterp.p[0]*alpha[i]*alpha[i] + LCInterp.p[1]*alpha[i] + LCInterp.p[2]) - sqrt(R2[i]);
      chisq += b*b;
      l1err += fabs(b);
    }
  
  if(ThisTask == 0)
    printf("LC R2(alpha) fit is %g*alpha^2 + %g*alpha + %g, RMS,L1 err in R(alpha) = %g|%g (# of pts = %d).\n",
	   LCInterp.p[0],LCInterp.p[1],LCInterp.p[2],sqrt(chisq/LC_NQUADFIT),l1err/LC_NQUADFIT,LC_NQUADFIT);  
}

/*
  This routine checks if a particle has crossed the light cone surface.  The parameter imageCovergae sets 
  how many images outside of the 8 cubes around 0,0,0 are checked.
*/
void check_particle(int part, float LCPrevPos[3], int imageCoverage)
{
  int i,j,k;
  float fbox;
  double R2_LCPrev,R2_LCCurr;
  float posShift[3];
  float LCCurrPos[3];
  float LCCurrPosShifted[3],LCPrevPosShifted[3];
  int lowlim_boxloop,uplim_boxloop;
  double t0,t1;
  float LCPos[3];
  double alpha,vdt2,a,b,c,q,aexpn,xdotdx,dtabsval;
  
#ifdef LIGHTCONE_INTERNAL_TIMER
  t0 = second();
#endif
  
  fbox = (float) All.BoxSize;

  vdt2 = LCInterp.dt_drift*LCInterp.dt_drift*(P[part].Vel[0]*P[part].Vel[0] + P[part].Vel[1]*P[part].Vel[1] + P[part].Vel[2]*P[part].Vel[2]);
  
  /* Algorithm for light cone generator - Matthew R. Becker and Andrey V. Kravtsov, U of C 2011
    1) We first make sure the particle's previous position (LCPrevPos) is in the primary image domain (i.e. [0,L]).
       We also make sure to shift the new position (LCCurrPos) the same amount (might not be to the same image!).
    2) Then we move both positions into each of the periodic images in 4X4X4 cube around the primary image.
       The primary image is at (2,2,2) with each box in cube labeled from 0-3 (i.e. it's zero index).
    2a) In each of these cubes, we check if the particle has crossed the light cone surface with origin at (0,0,0).
    2b) Then we add it to the appropriate buffer if it has crossed the light cone surface.
    
    NOTE: We are no longer marking if the particle has been output already since a particle 
    could appear twice in the light cone and that is OK.
    
    NOTE: The parameter imageCovergae controls how many extra layers of cubes about the 8 octants around (0,0,0) 
    are used.  
  */
  
  //1) 
  LCCurrPos[0] = P[part].Pos[0];
  LCCurrPos[1] = P[part].Pos[1];
  LCCurrPos[2] = P[part].Pos[2];
  
  for(j=0;j<3;j++)
      {
        while(LCPrevPos[j] < 0)
          {
	    LCPrevPos[j] += fbox;
	    LCCurrPos[j] += fbox;
	  }
	
	while(LCPrevPos[j] >= fbox)
          {
	    LCPrevPos[j] -= fbox;
	    LCCurrPos[j] -= fbox;
	  }
      }

  //2)
  lowlim_boxloop = -1 - imageCoverage;
  uplim_boxloop = 0 + imageCoverage;
  for(i=lowlim_boxloop;i<=uplim_boxloop;++i)
    {
      posShift[0] = i*fbox;
      
      for(j=lowlim_boxloop;j<=uplim_boxloop;++j)
	{
	  posShift[1] = j*fbox;
	  
	  for(k=lowlim_boxloop;k<=uplim_boxloop;++k)
	    {
	      posShift[2] = k*fbox;
	      
	      LCPrevPosShifted[0] = LCPrevPos[0]+posShift[0];
	      LCPrevPosShifted[1] = LCPrevPos[1]+posShift[1];
	      LCPrevPosShifted[2] = LCPrevPos[2]+posShift[2];
	      
	      R2_LCPrev =
		LCPrevPosShifted[0]*LCPrevPosShifted[0] + 
		LCPrevPosShifted[1]*LCPrevPosShifted[1] + 
		LCPrevPosShifted[2]*LCPrevPosShifted[2];
	      
	      LCCurrPosShifted[0] = LCCurrPos[0]+posShift[0];
	      LCCurrPosShifted[1] = LCCurrPos[1]+posShift[1];
	      LCCurrPosShifted[2] = LCCurrPos[2]+posShift[2];
	      
	      R2_LCCurr =
		LCCurrPosShifted[0]*LCCurrPosShifted[0] + 
		LCCurrPosShifted[1]*LCCurrPosShifted[1] + 
		LCCurrPosShifted[2]*LCCurrPosShifted[2];
	      
	      ++(LCProf.NumLCTests);
	      
	      //2a
	      if(R2_LCPrev < LCInterp.R2Prev && R2_LCCurr >= LCInterp.R2Curr) //particle is inside of LC and then goes outside of LC
		{
		  /* WOOT! WOOT! We found a particle in the light cone!
		     The LC generator will will try to use an O(a^3) quadrtaic approx to r(a) 
		     to do the LC interp.  If this cannot be done, then it reverts to an order O(a^2) method.
		     
		     Note that we always force the LC generator to put the LC inteprolated position between the prev and curr
		     positions along the velocity vector for either the O(a^3) or O(a^2) method.
		  */
		  
		  if(LCInterp.doLinearInterp)
		    {
		      xdotdx = LCPrevPosShifted[0]*LCCurrPosShifted[0] + LCPrevPosShifted[1]*LCCurrPosShifted[1] + LCPrevPosShifted[2]*LCCurrPosShifted[2]
			- R2_LCPrev;
		      alpha = (LCInterp.R2Prev - R2_LCPrev)/2.0/(xdotdx + LCInterp.R2Prev - LCInterp.RPrevRCurr);
		    }
		  else
		    {
		      a = vdt2 - LCInterp.p[0];
		      b = 2.0*LCInterp.dt_drift*(P[part].Vel[0]*LCPrevPosShifted[0] + P[part].Vel[1]*LCPrevPosShifted[1] + P[part].Vel[2]*LCPrevPosShifted[2]) 
			- LCInterp.p[1];
		      c = R2_LCPrev - LCInterp.p[2];
		      
		      q = -0.5*(b + b/fabs(b)*sqrt(b*b - 4.0*a*c));
		      if(fabs(q/a) > fabs(c/q))
			alpha = c/q;
		      else
			alpha = q/a;
		    }
		  
		  if(alpha > 1.0)
		    alpha = 1.0;
		  if(alpha < 0.0)
		    alpha = 0.0;
		  
		  LCPos[0] = alpha*LCCurrPosShifted[0] + (1.0-alpha)*LCPrevPosShifted[0];
		  LCPos[1] = alpha*LCCurrPosShifted[1] + (1.0-alpha)*LCPrevPosShifted[1];
		  LCPos[2] = alpha*LCCurrPosShifted[2] + (1.0-alpha)*LCPrevPosShifted[2];
                  
		  dtabsval = alpha*LCInterp.dt_drift + LCInterp.dt_a0;
		  aexpn = AofdT(dtabsval,LCInterp.i0,LCInterp.i1);
		  
		  /* never get here due to if statements above - commented out
		  if(alpha > 1.0 || alpha < 0.0)
		    {
		      printf("alpha is not correct for light cone interpolation: alpha = %g\n",alpha);
		      printf("r_i = %g, x_i = %g, r_i+1 = %g, x_i+1 = %g\n",sqrt(LCInterp.R2Prev),sqrt(R2_LCPrev),sqrt(LCInterp.R2Curr),sqrt(R2_LCCurr));
		      printf("x_i = %g|%g|%g, x_i+1 = %g|%g|%g\n",LCPrevPosShifted[0],LCPrevPosShifted[1],LCPrevPosShifted[2],
			     LCCurrPosShifted[0],LCCurrPosShifted[1],LCCurrPosShifted[2]);
		      printf("vdt2 = %g, v = %g|%g|%g, dt = %g\n",vdt2,P[part].Vel[0],P[part].Vel[1],P[part].Vel[2],LCInterp.dt_drift);
		      printf("alpha = %g, a,b,c, = %g|%g|%g, p = %g|%g|%g, s1,s2=%g|%g\n",alpha,a,b,c,LCInterp.p[0],LCInterp.p[1],LCInterp.p[2],q/a,c/q);
		      printf("LCPos = %g|%g|%g\n",LCPos[0],LCPos[1],LCPos[2]);
		      endrun(666);
		    }
		  */
		  
		  //going to cut the particles down to a 2x2x2 cube
		  if(fabs(LCPos[0]) <= fbox && fabs(LCPos[1]) <= fbox && fabs(LCPos[2]) <= fbox)
		    {
		      //2b
		      checkPutPartInLC(LCPos,part,aexpn);
		    } //cut down to 2x2x2 box
		} //cross LC surface check
	    } //k loop
	} //j loop
    } //i loop
  
#ifdef LIGHTCONE_INTERNAL_TIMER
  t1 = second();
  LCProf.CPU_CheckLC += timediff(t0,t1);
#endif
}

/*When a task has enough particles, dump the information to disk*/
void dump_lc_buf(int LC)
{
  int i;
  FILE *fd_p, *fd_v, *fd_i;
  FILE *fd_a;
  char fname_p[500], fname_v[500], fname_i[500];
  char fname_a[500];
      
  /*actually write the particle*/
  sprintf(fname_p, "%slightcone_buff/lightcone%03d/p/%s_Lightcone.p.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
  sprintf(fname_v, "%slightcone_buff/lightcone%03d/v/%s_Lightcone.v.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
  sprintf(fname_i, "%slightcone_buff/lightcone%03d/i/%s_Lightcone.i.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
  sprintf(fname_a, "%slightcone_buff/lightcone%03d/a/%s_Lightcone.a.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
    
  if(!(fd_p = fopen(fname_p, "a")))
    {
      printf("Task %d can't open file `%s' for writing lightcone.\n",ThisTask, fname_p);
      endrun(123);
    }
  if(!(fd_v = fopen(fname_v, "a")))
    {
      printf("Task %d can't open file `%s' for writing lightcone.\n",ThisTask, fname_v);
      endrun(123);
    }
  if(!(fd_i = fopen(fname_i, "a")))
    {
      printf("Task %d can't open file `%s' for writing lightcone.\n",ThisTask, fname_i);
      endrun(123);
    }
  if(!(fd_a = fopen(fname_a, "a")))
    {
      printf("Task %d can't open file `%s' for writing lightcone.\n",ThisTask, fname_a);
      endrun(123);
    }
  for (i=0;i<NPart_in_buf[LC];i++)
    {
      my_fwrite(&out_particles[LC][i].Pos[0], 3*sizeof(float), 1, fd_p);
      my_fwrite(&out_particles[LC][i].Vel[0], 3*sizeof(float), 1, fd_v);
      my_fwrite(&out_particles[LC][i].ID, sizeof(long long), 1, fd_i);
      my_fwrite(&out_particles[LC][i].a, sizeof(float), 1, fd_a);
    }

  fclose(fd_p);
  fclose(fd_v);
  fclose(fd_i);
  fclose(fd_a);

  NPart_Written[LC] += NPart_in_buf[LC];
  NPart_in_buf[LC] = 0;
}

/*Calculates z given distance from a lookup table*/
double ZofR(double R)
{
  int i;
  double slope = 0.;

  for(i=1;i<NZofRSteps;i++)
    {
      if (ZofRTableR[i] <= R)
	break;
    }
  slope = (ZofRTableZ[i-1] - ZofRTableZ[i])/(ZofRTableR[i-1] - ZofRTableR[i]);
  return ZofRTableZ[i] + slope*(R-ZofRTableR[i]);
}			

double RofZ(double Z)
{
  int i;
  double slope = 0.;

  for(i=1;i<NZofRSteps;i++)
    {
      if (ZofRTableZ[i] <= Z)
	break;
    }
  slope = (ZofRTableR[i-1] - ZofRTableR[i])/(ZofRTableZ[i-1] - ZofRTableZ[i]);
  return ZofRTableR[i] + slope*(Z-ZofRTableZ[i]);
}

static double RofZGL(double Z)
{
  double r = 0.0;
  int i;
  
  for(i=0;i<LC_NGAULEG_STEPS;++i)
    r += glW[i]/hubble_function(1.0/(1.0 + Z*glX[i]));
  
  return Z*r*C/All.UnitVelocity_in_cm_per_s;
}

/*generates the lookup table for calculating z given r*/
void MakeZofRTable(void)
{
  int i;
  double Thisa,logTimeEnd = 0.0;
  
  for(i=0;i<NZofRSteps;++i)
    {
      Thisa = exp(logTimeBegin + ((logTimeEnd - logTimeBegin) / (NZofRSteps-1)) * i);
      ZofRTableA[i] = Thisa;
      ZofRTableZ[i] = 1./Thisa - 1.;
      ZofRTableR[i] = RofZGL(ZofRTableZ[i]);
    }
}

/* Task 0 reads in the individual task lightcone information and writes as a 
   Gadget output file */
void Finalize_LightCone(void)
{
  char fname_p[500], fname_v[500], fname_i[500];
  char sfx[50];
  int buf, LC, i;
  FILE *fp, *fv, *fi, *fd;
  float tPos[3];
  long long tId;
  long long NInLC[NLCS];
  long long NInLCThisTask[NLCS];
  struct io_header outheader;
  char fname[500];
  long lnum;
  
  int nprocgroup, groupTask, masterTask;
  
#ifdef LIGHTCONE_OUTPUT_AINT
  char fname_a[500];
#endif
  
  //vars to make sure only  All.NumFilesWrittenInParallel tasks write to disk at once
  nprocgroup = NTask / All.NumFilesWrittenInParallel;
  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;
  masterTask = (ThisTask / nprocgroup) * nprocgroup;
  
  if(ThisTask == 0)
    printf("\nFinalizing the LC output\n");
  
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))
	{
	  /*First we need to write out remainder of buffer*/
	  for(LC=0;LC<NLCS;LC++)
	    {
	      if(LCNum[LC] == 0)
		continue;
	      dump_lc_buf(LC);
	    }
	}
      
      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  //Count how many particles have been output 
  //use the vars in mem 
  //will check each file below to make sure it is consistent with this value - bug check for restart and or I/O errors
  MPI_Allreduce(NPart_Written, NInLC, NLCS, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(ThisTask == 0) 
    {
      for(LC=0;LC<NLCS;LC++)
	{
	  if(LCNum[LC] == 0)
	    continue;
	 
	  printf("Light cone %d will have %lld particles.\n", LC, NInLC[LC]);
	}
    }
  
  if(ThisTask == 0)
    printf("Starting to write out LC files...\n");
  
  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (masterTask + groupTask))
        {
	  
	  for(LC=0;LC<NLCS;LC++)
	    {
	      //Are we even interested in this lightcone geometry?
	      if(LCNum[LC] == 0)
		continue;
	      
	      //open the file
	      sprintf(sfx, "%03d", LC);
	      sprintf(fname, "%slightcone/lightcone%s/%s_Lightcone_%s.%d", All.OutputDir, sfx, All.SnapshotFileBase,sfx,ThisTask);
	      if(!(fd = fopen(fname, "w")))
		{
		  printf("can't open file `%s' for making final light cone.\n",fname);
		  endrun(123);
		}
	      
	      //write the header with buffer
	      outheader = header;
	      
	      for(i=0;i<6;++i)
		{
		  outheader.npart[i] = 0;
		  outheader.npartTotal[i] = 0;
		}
	      
	      outheader.npart[1] = (int) (NPart_Written[LC]);
	      if(NInLC[LC] >> 32)
		{
		  outheader.npartTotal[2] = NInLC[LC] >> 32;
		  lnum = outheader.npartTotal[2];
		  outheader.npartTotal[1] = NInLC[LC] - (lnum << 32);
		}
	      else
		{
		  outheader.npartTotal[1] = NInLC[LC];
		  outheader.npartTotal[2] = 0;
		}
	      
	      buf = 256;
	      my_fwrite(&buf, sizeof(int), 1, fd);
	      my_fwrite(&outheader, sizeof(struct io_header), 1, fd);
	      my_fwrite(&buf, sizeof(int), 1, fd);
	      
	      //the buffer for the position and velocity information
	      buf = outheader.npart[1]*3*sizeof(float);

	      //read in position information and write it out
	      my_fwrite(&buf, sizeof(int),1,fd);
	      sprintf(fname_p, "%slightcone_buff/lightcone%03d/p/%s_Lightcone.p.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
	      if(!(fp = fopen(fname_p, "r")))
		{
		  printf("can't open file `%s' for finishing light cone.\n",fname_p);
		  endrun(123);
		}
	      NInLCThisTask[LC] = 0;
	      while(1)
		{
		  fread(&tPos[0], sizeof(float), 3, fp);
		  if (feof(fp))
		    break;
		  my_fwrite(&tPos[0], sizeof(float), 3, fd);
		  NInLCThisTask[LC]++;
		}
	      fclose(fp);
	      my_fwrite(&buf, sizeof(int),1,fd);
	      if(NInLCThisTask[LC] != NPart_Written[LC])
		{
		  printf("LIGHT CONE CORRUPTION: # of particles in file '%s' (%lld particles) does not match count in memory (%lld particles)!\n",
			 fname_p,NInLCThisTask[LC],NPart_Written[LC]);
		  endrun(666);
		}
	      
	      //read in velocity information and write it out
	      my_fwrite(&buf, sizeof(int),1,fd);
	      sprintf(fname_v, "%slightcone_buff/lightcone%03d/v/%s_Lightcone.v.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
	      if(!(fv = fopen(fname_v, "r")))
		{
		  printf("can't open file `%s' for finishing light cone.\n",fname_v);
		  endrun(123);
		}
	      NInLCThisTask[LC] = 0;
              while(1)
                {
                  fread(&tPos[0], sizeof(float), 3, fv);
                  if (feof(fv))
                    break;
                  my_fwrite(&tPos[0], sizeof(float), 3, fd);
                  NInLCThisTask[LC]++;
                }
              fclose(fv);
              my_fwrite(&buf, sizeof(int),1,fd);
              if(NInLCThisTask[LC] != NPart_Written[LC])
                {
                  printf("LIGHT CONE CORRUPTION: # of particles in file '%s' (%lld particles) does not match count in memory (%lld particles)!\n",
                         fname_v,NInLCThisTask[LC],NPart_Written[LC]);
                  endrun(666);
                }
	      
	      //and finally the ids
	      buf = ((int) outheader.npart[1])*sizeof(tId);
	      my_fwrite(&buf, sizeof(int),1,fd);
	      sprintf(fname_i, "%slightcone_buff/lightcone%03d/i/%s_Lightcone.i.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
	      if(!(fi = fopen(fname_i, "r")))
		{
		  printf("can't open file `%s' for finishing light cone.\n",fname_i);
		  endrun(123);
		}
	      NInLCThisTask[LC] = 0;
              while(1)
                {
		  fread(&tId, sizeof(tId), 1, fi);
                  if (feof(fi))
                    break;
		  my_fwrite(&tId, sizeof(tId), 1, fd);
                  NInLCThisTask[LC]++;
                }
              fclose(fi);
              my_fwrite(&buf, sizeof(int),1,fd);
              if(NInLCThisTask[LC] != NPart_Written[LC])
                {
                  printf("LIGHT CONE CORRUPTION: # of particles in file '%s' (%lld particles) does not match count in memory (%lld particles)!\n",
                         fname_i,NInLCThisTask[LC],NPart_Written[LC]);
                  endrun(666);
                }
	            
#ifdef LIGHTCONE_OUTPUT_AINT
	      //read in a int information and write it out
	      buf = ((int) outheader.npart[1])*sizeof(float);
	      my_fwrite(&buf, sizeof(int),1,fd);
	      sprintf(fname_a, "%slightcone_buff/lightcone%03d/a/%s_Lightcone.a.%d.%d", All.OutputDir, LC, All.SnapshotFileBase, ThisTask, LC);
	      if(!(fv = fopen(fname_a, "r")))
		{
		  printf("can't open file `%s' for finishing light cone.\n",fname_a);
		  endrun(123);
		}
	      NInLCThisTask[LC] = 0;
              while(1)
                {
                  fread(&tPos[0], sizeof(float), 1, fv);
                  if (feof(fv))
                    break;
                  my_fwrite(&tPos[0], sizeof(float), 1, fd);
                  NInLCThisTask[LC]++;
                }
              fclose(fv);
              my_fwrite(&buf, sizeof(int),1,fd);

              if(NInLCThisTask[LC] != NPart_Written[LC])
                {
                  printf("LIGHT CONE CORRUPTION: # of particles in file '%s' (%lld particles) does not match count in memory (%lld particles)!\n",
                         fname_a,NInLCThisTask[LC],NPart_Written[LC]);
                  endrun(666);
                }
#endif
	      
	      //close the finalized light cone file
	      fclose(fd);
	      
	      //too verbose printf("Finished output for octant # %d on Task %d.\n", LC, ThisTask);
	    } //end of loop over light cones
	  
	} //end of if for tasks currently writing
      
      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }
}

#endif
