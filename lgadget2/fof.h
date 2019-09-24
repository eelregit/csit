
/*! \file fof.h
 *  \brief function prototypes for FoF code
 */

void fof_fof(int num);

void fof_check_cell(int p, int i, int j, int k);
double fof_periodic(double x);
double fof_periodic_wrap(double x);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);
void fof_find_groups(void);
void fof_course_binning(void);
void fof_find_minids(void);
void fof_import_ghosts(void);
int fof_link_accross(void);
void fof_compile_catalogue(void);
int fof_grid_compare(const void *a, const void *b);
void fof_exchange_id_lists(void);
