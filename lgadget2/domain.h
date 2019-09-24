/*! \file domain.h
 *  \brief function prototypes for domain decomposition
 */

void   DomainDecomposition(int geom_flag); 

void domain_shiftSplit(void);

void domain_decomposeType(int geom_flag);
void domain_countToGo(void);
int  domain_findSplit(int cpustart, int ncpu, int first, int last);
int domain_findSplit_equalvolume(int cpustart, int ncpu, int first, int last);
void domain_findExchangeNumbers(int task, int partner, int *send, int *recv);
void domain_exchangeParticles(int partner, int send_count, int recv_count);
void domain_sumCost(void);
void domain_findExtent(void);



void domain_freebuffer(void);
void domain_allocatebuffer(void);








