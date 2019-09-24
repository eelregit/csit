
#ifdef MEMWATCH
#include "memwatch.h"
#endif

/*! \file allvars.h
 *  \brief declares global variables.
 *
 *  This file declares all global variables. Further variables should
 *  be added here, and declared as 'extern'. The actual existence of
 *  these variables is provided by the file 'allvars.c'. To produce
 *  'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>


#define PERIODIC
#ifndef PMGRID
#define PMGRID 64
#endif


#define  GADGETVERSION   "Lean-Gadget 2.0"  /*!< code version string */

#define  TIMEBASE    (1<<28)  /*!< The simulated timespan is mapped onto the
                               *   integer interval [0,TIMESPAN], where
                               *   TIMESPAN needs to be a power of 2. Note
                               *   that (1<<28) corresponds to 2^29.
                               */	


#define  MAX_REAL_NUMBER  1e37
#define  MIN_REAL_NUMBER  1e-37


/* ... often used physical constants (cgs units) */
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.99792458e10      //2.9979e10 - Matt Becker changed this to the exact SI value
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18	/* in h/sec */
#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifndef ASMTH
/*! ASMTH gives the scale of the short-range/long-range force split in
  units of FFT-mesh cells */
#define ASMTH 1.25	
#endif
#ifndef RCUT
/*! RCUT gives the maximum distance (in units of the scale used for
 * the force split) out to which short-range forces are evaluated in
 * the short-range tree walk.
 */
#define RCUT  4.5	
#endif

#ifndef LINKLENGTH
#define LINKLENGTH 0.2          /*!< FoF linking length in units of mean particle spaceing */
#endif

#ifndef GROUP_MIN_LEN
#define GROUP_MIN_LEN 32        /*!< FoF groups with fewer particles than this are not stored */
#endif

#define MAXLEN_OUTPUTLIST 350	/*!< maxmimum number of entries in output list */

#define RNDTABLE 256            /*!< length of table with precomputed random
                                     numbers. This is introduced to make sure
                                     that always the same random number is
                                     used for a given particle independent of
                                     the number of processors. */

#define DRIFT_TABLE_LENGTH  1000  /*!< length of look-up table for cosmological drift/kick factor */


#ifdef DOUBLEPRECISION
#define FLOAT double
#else
#define FLOAT float
#endif

#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
typedef unsigned short int uint4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
typedef unsigned int uint4byte;
#endif

typedef  long long  peanokey;
 

#ifndef DOMAINLEVELS
#define DOMAINLEVELS 4                  /*!< Number of levels of top-level tree (may not be larger than 10)*/
#endif

#ifndef HASHBITS
#define HASHBITS 6                      /*!< Number of bits used per dimension to construct grid of hash-cells,
					     must be smaller than peano bits per dimension */
#endif

#define PEANO_BITS_PER_DIMENSION 18	/*!< Number of bits used for construction of Peano-Hilbert keys.
                                             If the type 'peanokey' is defined only as int, 10 is the maximum 
					     that can be used. For a long int, one can go larger.  */

#define DOMAINGRID (1 << DOMAINLEVELS)              /*!< Dimension of effective top-level grid used to partition the global BH tree */
#define PEANOGRID  (1 << PEANO_BITS_PER_DIMENSION)  /*!< Dimension of fiducial Peano-Hilbert grid */

#define TAG_N          10   /*!< Various tags used for labelling MPI messages */ 
#define TAG_PDATA      11
#define TAG_WORK       12
#define TAG_MAXMIN     13
#define TAG_FLAG       14
#define TAG_FOF_IMPORT 15
#define TAG_FOF_LINK_A 16
#define TAG_FOF_LINK_B 17
#define TAG_FOF_GRAV_A 18
#define TAG_FOF_GRAV_B 19
#define TAG_FOF_PSEUDO 20
#define TAG_FOF_EXCH   21
#define TAG_TWOPOINT   22
#define TAG_GRAV_A     23
#define TAG_GRAV_B     24
#define TAG_PSEUDO     25
#define TAG_PM_A       26
#define TAG_PM_B       27
#define TAG_PM_FOLD    28



extern double DomainCenter[3],   /*!< Centre of global BH tree. For L-Gadget2 equal to the box centre. */
              DomainLen,         /*!< Side-length of root-cell of global BH tree. For L-Gadget2 equal to the box size. */
              DomainFac;         /*!< Factor to convert coordinates to cell numbers in top-level grid */
extern int DomainMyStart,        /*!< First cell in top-level grid assigned to this processor */
           DomainMyLast;         /*!< Last cell in top-level grid assigned to this processor */
extern int *DomainStartList,     /*!< List of first top-level cells assigned to each processor */ 
           *DomainEndList;       /*!< List of last top-level cells assigned to each processor */ 

extern double DomainWork[DOMAINGRID * DOMAINGRID * DOMAINGRID];   /*!< used to sum work in each top-level cell */
extern int DomainCount[DOMAINGRID * DOMAINGRID * DOMAINGRID];     /*!< used to count particles in each top-level cell */
extern int DomainTask[DOMAINGRID * DOMAINGRID * DOMAINGRID];      /*!< tells for each top-level cell to which processor it was assigned to */
extern int DomainNodeIndex[DOMAINGRID * DOMAINGRID * DOMAINGRID]; /*!< points to the node in the BH tree that holds this top-level cell's particle data */
extern int DomainNextnode[DOMAINGRID * DOMAINGRID * DOMAINGRID];  /*!< for pseudo-particles, this field is used to point to the next tree node in the tree walk */
extern int DomainPeanoMap[DOMAINGRID * DOMAINGRID * DOMAINGRID];  /*!< this table holds precomputed Peano-Hilbert keys for the top-level grid */



extern double RndTable[RNDTABLE];  /*!< This is a table with random numbers */



extern struct DomainNODE   /*!< Holds the multipole moment of one cell in the
                             top-level grid. Needed for data exchange, such
                             that a global top-level tree can be
                             constructed */
{
  FLOAT s[3];
  int4byte mass;
}
DomainMoment[DOMAINGRID * DOMAINGRID * DOMAINGRID];




extern int ThisTask;		/*!< the local processors  */
extern int NTask,               /*!< the number of processors  */
           PTask;	        /*!< the smallest integer such that NTask <= 2^PTask */

extern double CPUThisRun;	/*!< Sums CPU time of current process */

extern int RestartFlag;         /*!< Flags the kind of start the code is
                                  doing: Can be either (0) for start from
                                  intial conditions, (1) for resuming a run
                                  from restart files, or (2) for continuing
                                  from an old snapshot file */


extern gsl_rng *random_generator; /*!< the random number generator used */

extern int Flag_FullStep;      	/*!< Flags whether the current step is a full timestep */

extern int NumPart;		/*!< Number of particles on the local processor */



extern char ParameterFile[100];  /*!<  holds name of parameterfile */

extern FILE *FdInfo,             /*!<  file handle for "info.txt" file */
            *FdEnergy,           /*!<  file handle for "energy.txt" file */
            *FdTimings,          /*!<  file handle for "timings.txt" file which collects statistics on the gravitational tree computation */
            *FdCPU;              /*!<  file handle for log-file for CPU consumption in various code parts */
#ifdef DARKENERGY
extern FILE *FdDE;  /*!< file handle for darkenergy.txt log-file. */
#endif


extern double DriftTable[DRIFT_TABLE_LENGTH];       /*!<  table for cosmological drift factors */
extern double GravKickTable[DRIFT_TABLE_LENGTH];    /*!<  table for cosmological kick factors */


extern char *CommBuffer;	/*!< communication buffer, used at a number of places */

/*! this structure contains data which is the SAME for all tasks
 * (mostly code parameters read from the parameter file).  Holding
 * this data in a structure is convenient for writing/reading the
 * restart file, and it allows the introduction of new global
 * variables in a simple way. The only thing to do is to introduce
 * them into this structure.
 */
extern struct global_data_all_processes
{
  unsigned long long TotNumPart;		/*!<  total particle numbers (global value) */

  int MaxPart;			/*!< This gives the maxmimum number of
				   particles that can be stored on one
				   processor. */

  int ICFormat;			/*!< selects different versions of IC
				   file-format */

  int SnapFormat;		/*!< selects different versions of snapshot-file */


  int NumFilesWrittenInParallel;  /*!< puts an upper limit to the number of files the code may write/read simultaneously */

  int BufferSize;		/*!< size of communication buffer in MB */
  int MaxMemSize;               /*!< size of maximum memory consumption in MB */ 
  int BunchSizeForce;		/*!< number of particles fitting into the
				   buffer during the short-range gravity computation */
  int BunchSizeDomain;          /*!< number of particles fitting into a
				   buffer during the domain decomposition */
  int BunchSizeTwoPoint;	/*!< number of particles fitting into the
				   buffer during the computation of the two-point correlation function  */

  double PartAllocFactor;	/*!< in order to maintain work-load
				   balance, the particle load will
				   usually NOT be balanced.  Each
				   processor allocates memory for
				   PartAllocFactor times the average
				   number of particles to allow for
				   that */

  double TreeAllocFactor;	/*!< similarly for the tree: each processor
				   allocates a number of nodes which is
				   TreeAllocFactor times the maximum(!)
				   number of particles.  Note: A typical
				   local tree for N particles needs usually
				   ~0.65*N nodes */


  /* some force counters  */

  long long TotNumOfForces;	/*!< counts total number of force
				   computations  */


  /* system of units  */

  double UnitTime_in_s,         /*!< converts internal time unit to [sec/h] if desired */
    UnitMass_in_g,              /*!< converts internal mass unit to [gram/h] if desired */
    UnitVelocity_in_cm_per_s,   /*!< converts internal velocity unit to [cm/sec] if desired */
    UnitLength_in_cm,           /*!< converts internal length unit to [cm/h] if desired */
    UnitDensity_in_cgs,         /*!< converts internal density unit to [gram/cm^3 * h^2] if desired */
    UnitEnergy_in_cgs,
  GravityConstantInternal,      /*!< normally zero, otherwise used to override the physical value for G */
  G;                            /*!< gravity constant in internal code units */

  /* Cosmology */

  double Hubble;                 /*!< H_0 in internal code units. Note: this does NOT depend on h */
  double Omega0,                 /*!< total matter density relative to critical at z=0 */
         OmegaLambda,            /*!< dark energy density relative to critical at z=0 */
         OmegaBaryon,            /*!< baryon density relative to critical at z=0 */
         HubbleParam;	/*!< little `h', i.e. Hubble constant in units of 100
			  km/s/Mpc.  Note: Not actually used anywhere in
			  purely gravitational dynamics, i.e. only needed to
			  get physical values for cooling physics, for
			  example.
			 */

  double BoxSize;               /*!< Boxsize of periodic simulation volume */

  /* Code options */

  int ResubmitOn;                /*!< If set, the code will try to resubmit itself to a queuing system when it interrupts
				   a run because the CPU time limit is reached. */
  int TypeOfOpeningCriterion;    /*!< Gives the type of the timestep criterion that is used. Only type "O" is supported in L-Gadget2 */
  int OutputListOn;              /*!< If this is set, the desired output times are taken from a provided file */

  int TwoPointFlag;
  int PowerSpecFlag;

  /* parameters determining output frequency */

  int SnapshotFileCount;         /*!< This counts how many snapshot files have
                                  * already been produced for the run and is
                                  * used to number the snapshots in
                                  * consecutive order 
				  */

  double TimeBetSnapshot,
    TimeOfFirstSnapshot, CpuTimeBetRestartFile, TimeLastRestartFile, TimeBetStatistics, TimeLastStatistics;

  int NumCurrentTiStep;

  /* Current time of the simulation, global step, and end of simulation */

  double Time, TimeBegin, TimeStep, TimeMax;	/* marks end of the simulation */

  /* variables for organizing discrete timeline */

  double Timebase_interval;
  int Ti_Current;
  int Ti_nextoutput;


  int PM_Ti_endstep, PM_Ti_begstep;
  double Asmth[2], Rcut[2];
  double Corner[2][3], Xmintot[2][3], Xmaxtot[2][3];
  double TotalMeshSize[2];



  /* variables that keep track of cumulative CPU consumption */

  double TimeLimitCPU;
  double CPU_TreeConstruction;
  double CPU_TreeWalk;
  double CPU_Gravity;
  double CPU_Potential;
  double CPU_Domain;
  double CPU_Snapshot;
  double CPU_Total;
  double CPU_CommSum;
  double CPU_Imbalance;
  double CPU_Predict;
  double CPU_PM;
  double CPU_Peano;
  double CPU_FoF;
  double CPU_TwoPoint;
  double CPU_PowerSpec;


  /* tree code opening criterion */

  double ErrTolTheta;		/* BH-opening angle */
  double ErrTolForceAcc;	/* for new opening criterion */


  /* adjusts accuracy of time-integration */

  double ErrTolIntAccuracy;	/* for 1/a^{1/2} collisionless timestep
				   criterion */
  double MaxRMSDisplacementFac;

  double MinSizeTimestep, MaxSizeTimestep;

  double Dt_displacement;


  /* frequency of tree reconstruction/domain decomposition */


  /* gravitational and hydrodynamical softening lengths 
   * (given in terms of an `equivalent' Plummer softening length) 
   *
   * five groups of particles are supported 
   * 0=gas,1=halo,2=disk,3=bulge,4=stars 
   */

  double ComovSoftening;
  double Softening;
  double SofteningMaxPhys;

  /* If particle masses are all equal for one type, 
   *  the corresponding entry in MassTable is set to this value,
   * allowing the size of the snapshot files to be reduced
   */

  double PartMass;

  /* some filenames */
  char InitCondFile[200],
    OutputDir[200],
    SnapshotFileBase[100],
    EnergyFile[100],
    CpuFile[100], InfoFile[100], TimingsFile[100], RestartFile[100], ResubmitCommand[100], OutputListFilename[100];

  double OutputListTimes[MAXLEN_OUTPUTLIST];	/* was 200 in earlier version */
  int OutputListLength;

#ifdef DARKENERGY
  double DarkEnergyParam;       /*!< fixed w for equation of state */
#ifdef TIMEDEPDE
  char DarkEnergyFile[100];     /*!< tabelized w for equation of state */
#endif
#endif
#ifdef RESCALEVINI
  double VelIniScale;           /*!< Scale the initial velocities by this amount */
#endif


}
All;



/*! This structure holds all the information that is stored for each particle
 *  of the simulation.
 */
extern struct particle_data
{
  FLOAT Pos[3];			/*!< particle position at its current time */
  FLOAT Vel[3];			/*!< particle velocity at its current time */
  float OldAcc;			/*!< magnitude of old gravitational force. Used in relative opening criterion */
  float GravCost;		/*!< weight factor used to balance the work-load */
  int4byte Ti_endstep;          /*!< marks the beginning of this particles' timestep on the integer timeline */
  int4byte Ti_begstep;          /*!< marks the end of this particles' timestep on the integer timeline */
  long long ID;			/*!< unique particle identifier */
}
 *P,                            /*!< points to particles on this processor */
*DomainPartBuf;                 /*!< points to intermediate storage of particles during domain decomposition */


extern char *Exportflag;        /*!< in this auxiliary array, the tasks to which a given particles has to be exported are marked */


/*! This structure holds data that is exchanged in the communication during the gravity 
 * computation.
 */
extern struct gravdata_in
{
  union
  {
    FLOAT Pos[3];
    FLOAT Acc[3];
  }
  u;
 
  FLOAT OldAcc;
}
*GravDataIn, *GravDataGet, *GravDataResult, *GravDataOut;


/*! This structure holds auxiliary data that is used in the communication
 *algorithm used for the gravity computation.
 */
extern struct gravdata_index
{
  int Task;
  int Index;
  int SortIndex;
}
 *GravDataIndexTable;


/*! This structure holds the (partial) results of the gravitational force
 *  computation for particles that have been sent by other processors.
 */
extern struct gravdata_acctable
{
  FLOAT Acc[3];
  int Pindex;
}
 *GravDataAccTable;


/*! This structure holds data needed for the two-point correlation 
 *  function code.
 */
extern struct twopointdata_in
{
  FLOAT Pos[3];
  FLOAT Rs;
  int Task;
}
*TwoPointDataIn, *TwoPointDataGet;


/*! Header for the standard file format.
 */
extern struct io_header
{
  uint4byte npart[6];      /*!< npart[1] gives the number of particles in the present file, other particle types are ignored */
  double mass[6];          /*!< mass[1] gives the particle mass */
  double time;             /*!< time (=cosmological scale factor) of snapshot */
  double redshift;         /*!< redshift of snapshot */
  int4byte flag_sfr;       /*!< flags whether star formation is used (not available in L-Gadget2) */
  int4byte flag_feedback;  /*!< flags whether feedback from star formation is included */
  uint4byte npartTotal[6]; /*!< npart[1] gives the total number of particles in the run. If this number exceeds 2^32, the npartTotal[2] stores
                                the result of a division of the particle number by 2^32, while npartTotal[1] holds the remainder. */
  int4byte flag_cooling;   /*!< flags whether radiative cooling is included */
  int4byte num_files;      /*!< determines the number of files that are used for a snapshot */
  double BoxSize;          /*!< Simulation box size (in code units) */
  double Omega0;           /*!< matter density */
  double OmegaLambda;      /*!< vacuum energy density */
  double HubbleParam;      /*!< little 'h' */
  int4byte flag_stellarage;     /*!< flags whether the age of newly formed stars is recorded and saved */
  int4byte flag_metals;         /*!< flags whether metal enrichment is included */
  int4byte hashtabsize;         /*!< gives the size of the hashtable belonging to this snapshot file */
  //char fill[84];		/*!< fills to 256 Bytes */
  unsigned int npartTotalHighWord[6];  /*!< High word of the total number of par
ticles of each type */
  char fill[60];

}
header;  /*!< Global variable to hold the snapshot file header */



/*
 * Variables for Tree
 * ------------------
 */

/*! structure for tree node data */
extern struct NODE
{
  FLOAT len;			/*!< sidelength of treenode */
  FLOAT center[3];		/*!< geometrical center of node */
  union
  {
    int4byte suns[8];		/*!< temporary pointers to daughter nodes */
    struct
    {
      FLOAT s[3];               /*!< center of mass of node */
      int4byte mass;            /*!< mass of node */
      int4byte cost;            /*!< counts the number of interactions in which this node is used */
      int4byte sibling;         /*!< this gives the next node in the walk in case the current node can be used */
      int4byte nextnode;        /*!< this gives the next node in case the current node needs to be opened */
      int4byte father;          /*!< this gives the parent node of each node (or -1 if we have the root node) */
    }
    d;
  }
  u;
}
*Nodes_base,                    /*!< points to the actual memory allocted for the nodes */
*Nodes;                         /*!< this is a pointer used to access the nodes which is shifted such that Nodes[All.MaxPart] 
				  gives the first allocated node */



extern int MaxNodes;		/*!< maximum allowed number of internal tree-nodes on one processor */

extern int4byte *Nextnode;	/*!< gives next node in tree walk for each particle */
extern int4byte *Father;        /*!< gives parent node in the tree for each particle */


#ifdef LIGHTCONE

//#define LIGHTCONE_INTERNAL_TIMER /* define this to turn on profiling of just the light cone step, */
                                 /* as opposed to the prediction step only */

extern FILE *FdLC; /*!<  file handle for "lightcone.txt" file */
#define NLCS 9
#define NZofRSteps   100000
#define LC_OUTPUT_BUF 10000
#define LC_NQUADFIT 25
#define LC_NGAULEG_STEPS 500 //gives almost exact result out to z=999 - do not change

extern double ZofRTableR[NZofRSteps];
extern double ZofRTableZ[NZofRSteps];
extern double ZofRTableA[NZofRSteps];

extern struct LC_Data
{
  float Pos[3];
  float Vel[3];
  long long ID;
  float a;
} out_particles[NLCS][LC_OUTPUT_BUF];
extern long long NPart_in_buf[NLCS];
extern int LCNum[NLCS];
extern long long NPart_Written[NLCS];

extern struct LCProf_str
{
  long long NumLCTests;
  long long NumLCFound;
  double CPU_CheckLC;
  double CPU_IO;
} LCProf;

extern struct LCInterp_str
{
  int doLinearInterp;
  double RPrevRCurr;
  double loga0;
  double loga1;
  double a0;
  double a1;
  double R2Prev;
  double R2Curr;
  double dt_drift;
  double dt_a0;
  int i0;
  int i1;
  double p[3];
} LCInterp;
#endif

#endif

