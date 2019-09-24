#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef HPM
#include <libhpm.h>
#endif



/*! \file main.c
 *  \brief start of the program
 */

/*!
 *  This function initializes the MPI communication packages, and sets
 *  cpu-time counters to 0.  Then begrun() is called, which sets up
 *  the simulation either from IC's or from restart files.  Finally,
 *  run() is started, the main simulation loop, which iterates over
 *  the timesteps.
 */
int main(int argc, char **argv)
{
  double t0, t1;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

#ifdef HPM
  hpmInit(ThisTask, "L-Gadget2");
  hpmStart(1, "Main");
#endif


  if(NTask <= 1)
    {
      if(ThisTask == 0)
	printf("Number of processors MUST be a larger than 1.\n");
      endrun(0);
    }

  for(PTask = 0; NTask > (1 << PTask); PTask++);

  if(argc < 2)
    {
      if(ThisTask == 0)
	{
	  printf("Parameters are missing.\n");
	  printf("Call with <ParameterFile> [<RestartFlag>]\n");
	}
      endrun(0);
    }

  strcpy(ParameterFile, argv[1]);

  if(argc >= 3)
    RestartFlag = atoi(argv[2]);
  else
    RestartFlag = 0;


  All.CPU_TreeConstruction = All.CPU_TreeWalk = All.CPU_Gravity = All.CPU_Domain =
    All.CPU_Snapshot = All.CPU_Total = All.CPU_CommSum = All.CPU_Imbalance =
    All.CPU_Predict = All.CPU_FoF = All.CPU_PM = All.CPU_Peano = All.CPU_TwoPoint = All.CPU_PowerSpec = 0;

  CPUThisRun = 0;


#ifdef DEBUG
  enable_core_dumps_and_fpu_exceptions();
#endif

  t0 = second();

  begrun();			/* set-up run  */

  t1 = second();
  CPUThisRun += timediff(t0, t1);
  All.CPU_Total += timediff(t0, t1);

  run();			/* main simulation loop */

#ifdef HPM
  hpmStop(1);
  hpmTerminate(ThisTask);
#endif

  MPI_Finalize();		/* clean up & finalize MPI */

  return 0;
}




/* ----------------------------------------------------------------------
   The rest of this file contains documentation for compiling and
   running the code, in a format appropriate for 'doxygen'. To
   (re)generate the documentation, run 'doxygen' in this directory. 
   ----------------------------------------------------------------------
 */

/*! \mainpage Reference documentation for L-GADGET2

\author Volker Springel \n
        Max-Planck-Institute for Astrophysics \n
        Garching, Germany \n
        volker@mpa-garching.mpg.de \n
\n

<b>L-GADGET2</b> is a special "lean" version of GADGET2, optimized
for having low memory footprint (while still being efficient) when 
running pure N-body simulations of cosmological structure
formation. The key features of L-GADGET2 are as follows:

- L-GADGET2 is a TreePM gravity solver that uses an explicit
  long-range/short-range force split. Time integration is done by a
  generalized symplectic integration algorithm. The short-range
  timestep is individual and adaptive on a per-particle basis.

- The code allows particle numbers larger than 2^32=4.3e9 (64-bit
  arithmetic is used where needed for this).

- The code can be run for an arbitrary number of processors (provided one has
  at least 2). The results of the computation are in principle completely
  independent of the number of processors used. However, since the order in
  which arithmetic operations are carried out depends on the number of CPUs,
  floating point round-off can be different for different numbers of CPUs.

- The grid-size selected for the PM part of the algorithm is arbitrary (but
  needs to be at least 4 per dimension). Note however that the FFT will be
  fastest for mesh-dimensions that are products of small prime factors.

- The code contains a fast parallel FOF groupfinder, which is (optionally)
  run each time a snapshot file is produced.  The data of the snapshot
  files and the group catalogues is structured such that individual groups
  can be rapidly and selectively accessed, which is particularly important
  in simplifying the handling of very large simulation data sets. A small
  library with functions for accessing these data structures conveniently
  from within IDL is included with the source code. Also, examples for
  plotting the halo mass function or particles of individual groups in IDL
  are included.

- The code can measure the two-point correlation function upon request
  (each time an output is generated).

- Also, the code contains power spectrum estimation routines. For modes
  between the fundamental of the box and the Nyquist frequency of the PM
  mesh, the long-range mesh of the PM force computation is used, which
  essentially does not noy incur additional computational cost. For modes
  on smaller scales, the technique of "self-folding" is used, if desired in
  several statges and with a smaller Fourier mesh.

- The code can also be used to generate a 'glass-file'.

- Unlike GADGET2, the special version L-GADGET2 can only be used for 
  cosmological simulations of a collisionless component represented by
  equal-mass particles in periodic boxes.

Note that the documentation you are reading here is not exhaustive. It
serves as annotated and cross-referenced code-documentation, and the few
other pages included here mostly cover differences to the GADGET2
code. Familiarity with GADGET2 is hence assumed. Please refer to the
documentation of the latter for more detailed background information.

\n

\section install Compilation 

L-GADGET2 is written in standard ANSI-C, and should hence compile without
problem on every modern computer with a unix-like operating
system. However, L-GADGET2 needs the following non-standard libraries for
compilation:

- \b MPI - the Message Passing Interface (version 1.0 or higher). Many
  vendor supplied versions exist, in addition to excellent open source
  implementations, e.g.  MPICH
  (http://www-unix.mcs.anl.gov/mpi/mpich/) or LAM
  (http://www.lam-mpi.org/).

- \b GSL - the GNU scientific library. This open-source package can be
  obtained at http://www.gnu.org/software/gsl, for example. GADGET
  needs this library for a few simple cosmological
  integrations at start-up, and for random number generation.

- \b FFTW - the <em>Fastest Fourier Transform in the West</em>. This
  open-source package can be obtained at http://www.fftw.org, for
  example. Note that you need to use the MPI-enabled version 2.3.1 -- the
  new version 3.0 uses a different application interface and does not yet
  support MPI.

Note that FFTW needs to be compiled with parallel support enabled.  This
can be achieved by passing the option <b>--enable-mpi</b> to the configure
script. When at it, you might as well add <b>--enable-type-prefix</b> to
obtain the libraries in both a single and double precision version. If this
has not been done, you should set the option \b NOTYPEPREFIX_FFTW in
<b>L-Gadget2</b>'s \ref Gadget-Makefile "Makefile".

Note that if any of the above libraries is not installed in standard
locations on your system, the \ref Gadget-Makefile "Makefile" provided
with the code may need slight adjustments. Similarly, compiler
options, particularly with respect to optimisations, may need
adjustment to the C-compiler that is used. Finally, the \ref
Gadget-Makefile "Makefile" contains a number of compile-time options
that need to be set appropriately for the type of simulation that is
simulated.

The provided makefile is compatible with GNU-make, i.e. typing \b make
or \b gmake should then build the executable <b>L-Gadget2</b>.  
If your computer does not have GNU-make, you might have to modify 
the syntax of the makefile at a few places.

\section howtorun Running the code

In order to start the code, a \ref parameterfile "parameterfile" needs to
be specified as an argument to the executable. An optional second parameter
can be either 1 or 2, depending on the type of restart-option that is
desired. 

If "1" is selected, the code will be started from a set of restart
files, which are typically periodically generated by the running simulation
code, or when the run terminates, either because the specified CPU-time
limit is reached or because a file named "stop" was found in the output
directory. Generating the latter with a command like "echo > stop" is hence
a possibility to interrupt a running simulation gracefully.

If "2" is selected, the simulation can be restarted from a snapshot file,
which must be specified as initial conditions file in the parameterfile.
In this case, the starting redshift needs not be modified. Instead, it will
be taken from the header of the snapshot file. Also, the code will append
to all log files of the code, and carry on with the numbering of the
snapshot files beginning from the number of the file that was read. 

If no additional parameter is specified, the simulation is started from
scratch from the initial conditions file that is specified in the
parameterfile.

*/











/*! \page parameterfile  Parameterfile of L-GADGET2

- \b InitCondFile \n Filename of the initial conditions file. If it's
a multi-file snapshot, omit the trailing digits that distinguish
between the different parts of the snapshot file. The code will
automatically recognize the number of files, provided the file header has been
correctly set up.

- \b OutputDir                  \n  The output directory of the code.    
 

- \b EnergyFile                  \n        energy.txt

- \b InfoFile                    \n        info.txt

- \b TimingsFile                 \n        timings.txt

- \b CpuFile                    \n         cpu.txt

- \b RestartFile                 \n        restart

- \b SnapshotFileBase           \n         snapshot

- \b TimeLimitCPU                \n        430000

- \b ResubmitOn                  \n        0

- \b ResubmitCommand            \n         xyz

- \b NumFilesWrittenInParallel   \n        2

- \b TimeBegin                  \n         0.02

- \b TimeMax                    \n         1

- \b Omega0                     \n         0.3

- \b OmegaLambda               \n          0.7

- \b OmegaBaryon               \n          0

- \b HubbleParam               \n          0.7

- \b BoxSize                   \n          50000

- \b OutputListFilename        \n          parameterfiles/outputs.txt

- \b OutputListOn              \n          1

- \b TimeBetSnapshot           \n          0

- \b TimeOfFirstSnapshot       \n          0

- \b CpuTimeBetRestartFile     \n          7200

- \b TimeBetStatistics         \n          0.1

- \b ErrTolIntAccuracy         \n          0.05

- \b MaxSizeTimestep           \n          0.15

- \b MinSizeTimestep           \n          0

- \b MaxRMSDisplacementFac     \n          0.25

- \b ErrTolTheta               \n          0.5

- \b TypeOfOpeningCriterion    \n          1

- \b ErrTolForceAcc           \n           0.005

- \b PartAllocFactor         \n            2.5

- \b TreeAllocFactor         \n            0.7

- \b BufferSize              \n            150

- \b UnitLength_in_cm       \n             3.08568e+21

- \b UnitMass_in_g          \n             1.989e+43

- \b UnitVelocity_in_cm_per_s   \n         100000

- \b GravityConstantInternal    \n         0

- \b Softening            \n               7

- \b SofteningMaxPhys    \n                7

*/










/*! \page Gadget-Makefile  Makefile of L-GADGET2

Unlike GADGET2, the special version L-GADGET2 can only be used for
cosmological N-body simulations in periodic boxes. All particles are
assumed to be of equal mass, and the PM-Tree algorithm is always used.
The desired PM mesh-size must be passed at compile time via the
makefile.
 
Some key features of L-GADGET2 are controlled with compile-time
options in the makefile rather than by the parameterfile. This was
done in order to allow the generation of highly optimised binaries by
the compiler, even when the underlying source allows for many
different ways to run the code. Unfortunately, this technique has the
disadvantage that different simulations may require different binaries
of GADGET. If several simulations are run concurrently, there is hence
the danger that a simulation is started/resumed with the `wrong'
binary. Note that while GADGET checks the plausibility of some
important code options, this is not done for all of them. To minimise
the risk of using the wrong code for a simulation, my recommendation
is therefore to produce a separate executable for each simulation that
is run. For example, a good strategy is to make a copy of the whole
code together with its makefile in the output directory of each
simulation run, and then to use this copy to compile the code and to
run the simulation.

The makefile contains a dummy list of all available compile-time code
options, with most of them commented out by default. To activate a
certain feature, the corresponding parameter should be commented in,
and given the desired value, where appropriate. 

<b>Important Note:</b>Whenever you change one of the makefile options
described below, a full recompilation of the code is necessary. To
guarantee that this is done, you should give the command <b>make
clean</b>, which will erase all object files, followed by <b>make</b>.


- \b PMGRID=128 \n This sets the mesh-size used in the TreePM method,
   i.e. the long-range force is computed with a PM-algorithm, and the
   short range force with the tree. The parameter has to be set to the
   size of the mesh that should be used, e.g.~64, 96, 128, etc. The
   mesh dimensions need not necessarily be a power of two, but the FFT
   is fastest for such a choice. 

- \b FOF \n If this is set, the code will automaticall produce a
   friends-of-friends group catalogue for each snapshot file that is
   produced.

- \b WALLCLOCK \n If set, a wallclock timer is used by the code to
   measure internal time consumption (see cpu-log file).  Otherwise, a
   timer that measures consumed processor ticks is used. It is
   recommended to set this option unless you have to share CPUs with
   other processes.

- \b DOUBLEPRECISION \n This makes the code store and compute internal
   particle data in double precision. Note that output files are
   nevertheless written by converting the values that are saved to
   single precision.

- \b NOTREERND \n If this is not set, the tree construction will
   succeed even when there are a few particles at identical
   locations. This is done by `rerouting' particles once the node-size
   has fallen below \f$10^{-3}\f$ of the softening length. When this
   option is activated, this will be suppressed and the tree
   construction will always fail if there are particles at extremely
   close or identical coordinates.

- \b NOSTOP_WHEN_BELOW_MINTIMESTEP \n If this is
   activated, the code will not terminate when the timestep falls
   below the value of \b MinSizeTimestep specified in the
   parameterfile. This is useful for runs where one wants to enforce a
   constant timestep for all particles. This can be done by activating
   this option, and by setting \b MinSizeTimestep} and \b
   MaxSizeTimestep to an equal value.

- \b T3E \n The code assumes that \b sizeof(int)=4 holds. A few
   machines (like Cray T3E) have \b sizeof(int)=8. In this case, set
   the T3E flag.

- \b NOTYPEPREFIX_FFTW \n If this is set, the fftw-header/libraries
   are accessed without type prefix (adopting whatever was chosen as
   default at compile-time of fftw). Otherwise, the type prefix 'd'
   for double-precision is used.

*/
