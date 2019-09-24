/*! \file allvars.c
 *  \brief creates global variables.
 *
 *  This file creates all the global variables that are declared in allvars.h
 *  
 *  To produce 'allvars.c' from 'allvars.h', do the following:
 *
 *     - Erase all #define's
 *     - add #include "allvars.h" 
 *     - delete all keywords 'extern'
 *     - delete all struct definitions enclosed in {...}, e.g.
 *        "extern struct global_data_all_processes {....} All;"
 *        becomes "struct global_data_all_processes All;"
 */

#include "allvars.h"

double DomainCenter[3], DomainLen, DomainFac;
int DomainMyStart, DomainMyLast;
int *DomainStartList, *DomainEndList;

double DomainWork[DOMAINGRID * DOMAINGRID * DOMAINGRID];
int DomainCount[DOMAINGRID * DOMAINGRID * DOMAINGRID];
int DomainTask[DOMAINGRID * DOMAINGRID * DOMAINGRID];
int DomainNodeIndex[DOMAINGRID * DOMAINGRID * DOMAINGRID];
int DomainPeanoMap[DOMAINGRID * DOMAINGRID * DOMAINGRID];
int DomainNextnode[DOMAINGRID * DOMAINGRID * DOMAINGRID];

struct DomainNODE DomainMoment[DOMAINGRID * DOMAINGRID * DOMAINGRID];




int ThisTask;			/* the local processors  */
int NTask, PTask;		/* note: NTask = 2^PTask */

double CPUThisRun;		/* Sums CPU time of current process */


int RestartFlag;


gsl_rng *random_generator;      /* the random number generator used */

int Flag_FullStep;

int NumPart;			/* Note: this is the LOCAL processor value */




double RndTable[RNDTABLE];

/* variables for input/output ,  usually only used on process 0 
 */
char ParameterFile[100];
FILE *FdInfo, *FdEnergy, *FdTimings, *FdCPU;
#ifdef DARKENERGY
FILE *FdDE; // file handle for darkenergy.txt log-file
#endif


double DriftTable[DRIFT_TABLE_LENGTH];
double GravKickTable[DRIFT_TABLE_LENGTH];


char *CommBuffer;		/* communication buffer, used at a number of places */
char *Exportflag;

/* this structure contains data which is the SAME for all 
 * tasks (mostly code parameters read from the parameter file). 
 * Holding this data in a structure is convenient for writing/reading
 * the restart file, and it allows the introduction of new global
 * variables in a simple way. The only thing to do is to introduce them
 * into this structure.
 */
struct global_data_all_processes All;



/* The following structure holds all the information that is
 * stored for each particle of the simulation.
 */
struct particle_data *P, *DomainPartBuf;






/* Various structures for communication during the gravity 
 * computation.
 */



struct gravdata_in *GravDataIn, *GravDataGet, *GravDataResult, *GravDataOut;

struct gravdata_index *GravDataIndexTable;

struct gravdata_acctable *GravDataAccTable;



struct twopointdata_in *TwoPointDataIn, *TwoPointDataGet;









/* Header for the standard file format.
 */
struct io_header header;



/*
 * Variables for Tree
 * ------------------
 */


struct NODE *Nodes, *Nodes_base;



int MaxNodes;			/* maximum allowed number of internal nodes */



int4byte *Nextnode;		/* gives next node in tree walk  (nodes array) */
int4byte *Father;

#ifdef LIGHTCONE
FILE *FdLC;
double ZofRTableR[NZofRSteps];
double ZofRTableA[NZofRSteps];
double ZofRTableZ[NZofRSteps];
struct LC_Data out_particles[NLCS][LC_OUTPUT_BUF];
long long NPart_in_buf[NLCS];
int LCNum[NLCS];
long long NPart_Written[NLCS];
struct LCProf_str LCProf;
struct LCInterp_str LCInterp;
#endif
