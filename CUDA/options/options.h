#define MAXNFIBRES 4 
#define MAXNPARAMS 15
//NFIBRES*3 + S0 + d + d_std + f0 + rician = 15
#define MAXNDIR 640
//640 gradient directions
#define THREADS_BLOCK_MCMC 64
#define THREADS_BLOCK_FIT 64
#define MAXNDIRS_PER_THREAD 10 
//MAXNDIR/THREADS_BLOCK_MCMC


