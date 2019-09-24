
/*! \file proto.h
 *  \brief  contains function prototypes
 */


#ifndef ALLVARS_H
#include "allvars.h"
#endif
#include "forcetree.h"
#include "domain.h"
#include "fof.h"

void reallocate_particle_memory_MaxPart(void);
void reallocate_particle_memory_NumPart(void);
void pm_only(void);

double growth_int(double a, void *param);
double GrowthFactor(double astart, double aend);
double growth(double a);

#ifdef SYSMALLOC

#define mymalloc(n) malloc(n)
#define myfree(n) free(n)
#define myrealloc(p,n) realloc(p,n)

#else

void mymalloc_init(size_t n);
void *mymalloc(size_t n);
void myfree(void *p);
void *myrealloc(void *p, size_t n);

#endif

void rayleigh_save(void);
void prepare_rayleigh(void);

void write_pid_file(void);

void dump_potential(void);

void powerspec(int flag);
void powerspec_save(void);
double PowerSpec_Efstathiou(double k);

void foldonitself(void);
void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master, int *last);

void read_ic(char *fname);
void read_file(char *fname, int readTask, int lastTask);
int find_files(char *fname);

void enable_core_dumps_and_fpu_exceptions(void);
void  glass_step(void);

void seed_glass(void);


double compute_mean_rms_velocity(void);
int grav_tree_compare_key(const void *a, const void *b);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
int compare_key(const void *a, const void *b);
void free_commbuffers(void);

void move_particles(int time0, int time1);
void ngb_treefind_variable_check(FLOAT xyz[3], FLOAT hguess, int **ngblistback, int ngb_found);
void find_dt_displacement_constraint(void);
void save_snapshot(int num);

void find_next_sync_point_and_drift(void);
void kick_particle(int i, FLOAT *acc);

void  set_units_sfr(void);

void  gravity_forcetest(void);
void init_peano_map(void);

void   allocate_commbuffers(void);
void   allocate_memory(void);
void   begrun(void);
void   check_omega(void);
void   close_outputfiles(void);
void   compute_accelerations(int mode);
void   compute_global_quantities_of_system(void);
void   compute_potential(void);
void   construct_timetree(void);
void   cooling_and_starformation(void);
void   count_hot_phase(void);
void   delete_node(int i);
void   density(void);
void   density_decouple(void);
void   determine_interior(void);
int    dissolvegas(void);
double dmax(double,double);
double dmin(double,double);
void   do_box_wrapping(void);
double drand48();
double enclosed_mass(double R);
void   endrun(int);
void   energy_statistics(void);
void   ensure_neighbours(void);

void   every_timestep_stuff(void);
void    ewald_corr(double dx, double dy, double dz, double *fper);

void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_force_ni(int iii, int jjj, int kkk, double x[3], double force[3]);

void   ewald_init(void);
double ewald_psi(double x[3]);
double  ewald_pot_corr(double dx, double dy, double dz);
int    find_ancestor(int i);
int    find_next_outputtime(int time);
void   find_next_time(void);
int    find_next_time_walk(int node);
void   free_memory(void);
void   advance_and_find_timesteps(void);
int    get_timestep(int p, double *a, int flag);


double get_starformation_rate(int i);
void   gravity_tree(void);
void   hydro_force(void);
int    imax(int,int);
int    imin(int,int);
void   init(void);
void   init_clouds(void);
void   insert_node(int i);
void   integrate_sfr(void);
int    mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void   open_outputfiles(void);
void   peano_hilbert_order(void);
double pot_integrand(double xx);
void   predict(double time);
void   predict_collisionless_only(double time);
void   predict_sph_particles(double time);
void   prepare_decouple(void);
void   read_ic(char *fname);
void   read_ic_cluster(char *fname);
void   read_ic_cluster_gas(char *fname);
void   read_ic_cluster_wimp(char *fname);
int    read_outputlist(char *fname);
void   read_parameter_file(char *fname);
void   rearrange_particle_sequence(void);
void   reorder_gas(void);
void   reorder_particles(void);
void   restart(int mod);
void   run(void);
void   savepositions(int num);
void   savepositions_ioformat1(int num);
double second(void);
void   set_softenings(void);
void   set_sph_kernel(void);
void   set_units(void);
void   setup_smoothinglengths(void);

void   statistics(void);
double timediff(double t0,double t1);
void   veldisp(void);
void   veldisp_ensure_neighbours(int mode);

double get_random_number(int id);
void set_random_numbers(void);



double get_hydrokick_factor(int time0, int time1);
double get_gravkick_factor(int time0, int time1);
double drift_integ(double a, void *param);
double gravkick_integ(double a, void *param);
double hydrokick_integ(double a, void *param);
void   init_drift_table(void);
double get_drift_factor(int time0, int time1);

int ngb_clear_buf(FLOAT searchcenter[3], FLOAT hguess, int numngb);

void force_costevaluate(void);


/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */

double pow(double, double);




void long_range_init(void);
void long_range_force(void);
void pm_init_periodic(void);
void pmforce_periodic(void);
void pm_init_regionsize(void);
void pm_init_nonperiodic(void);
int  pmforce_nonperiodic(int grnr);

void pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

void force_treeevaluate_potential_shortrange(int target);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);



void twopoint(void);
void count_local(int target, int mode, double rs);
void twopoint_ngb_treefind_variable(FLOAT searchcenter[3], FLOAT rsearch, int mode);
int twopoint_compare_key(const void *a, const void *b);
void twopoint_save(void);


#ifdef LIGHTCONE
void setup_lightcone_indexing(void);
void setup_lightcone(void);
void MakeZofRTable(void);
double ZofR(double r);
double RofZ(double z);
void Finalize_LightCone(void);
void check_particle(int part, float LCPrevPos[3], int imageCoverage);
void dump_lc_buf(int LC);
void init_lc_interp(int tim0, int time1);
#endif

double INLINE_FUNC hubble_function(double a);
#ifdef DARKENERGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef TIMEDEPGRAV
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif /* TIMEDEPGRAV */
#ifdef EXTERNALHUBBLE
double INLINE_FUNC hubble_function_external(double a);
#endif /* EXTERNALHUBBLE */
#endif /* TIMEDEPDE */
#endif /* DARKENERGY */

