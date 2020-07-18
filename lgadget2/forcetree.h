
/*! \file forcetree.h
 *  \brief function prototypes for routines of tree code
 */


#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif
void force_update_hmax_of_parent_node(int no);

void force_treeupdate_toplevel(int no, int level, int *nextfree);
void force_exchange_pseudodata(void);
void force_treeinit(void);

double dIm_func(int, double, double);
double Im_func(int, double, double);


int force_treeevaluate_shortrange(int target, int mode, FLOAT *acc);
void force_insert_pseudo_particles(void);

int force_treeevaluate_ewald_correction(int target);
void force_create_empty_nodes(int no, int level, int nx, int ny, int nz, int *nodecount, int *nextfree);
void force_flag_localnodes(void);

void force_add_star_to_tree(int igas, int istar);

void   force_costevaluate(void);
int    force_getcost_single(void);
int    force_getcost_quadru(void);
void   force_resetcost(void);
void   force_setupnonrecursive(int no);
size_t force_treeallocate(int maxnodes, int maxpart);  
int    force_treebuild(void);
int    force_treebuild_single(int startnode, int *typelist);
void   force_treeevaluate(int target);
void   force_treeevaluate_direct(int target);
void   force_treeevaluate_single(int tree, int targetpart, double epsilon);
void   force_treeevaluate_single_BH(int tree, int targetpart, double epsilon);
void   force_treeevaluate_potential(int target);
void   force_treeevaluate_potential_single(int tree, int targetpart, double epsilon);
void   force_treeevaluate_potential_single_BH(int tree, int targetpart, double epsilon);
void   force_treefree(void);
void   force_update_node(int no, int flag);

void force_update_node_recursive(int no, int sib, int father);

void   force_update_size_of_parent_node(int no);

void   dump_particles(void);

FLOAT  INLINE_FUNC ngb_periodic(FLOAT x);
FLOAT  INLINE_FUNC ngb_periodic_longbox(FLOAT x);
FLOAT  ngb_select_closest(int k, int n, FLOAT *arr, int *ind);
void   ngb_treeallocate(int npart);
void   ngb_treebuild(void);
int    ngb_treefind_pairs(FLOAT xyz[3], FLOAT hsml);
int    ngb_treefind_variable(FLOAT xyz[3], FLOAT hguess);
void   ngb_treefree(void);
void   ngb_treesearch(int);
void   ngb_treesearch_pairs(int);
void   ngb_update_nodes(void);
void   ngb_treesearch_notsee(int no);















