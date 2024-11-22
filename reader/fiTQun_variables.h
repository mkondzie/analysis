//----------fiTQun variables----------//
int fqnmrfit;//                              Number of MR fit results that are available
int fqmrifit[50];//                    Fit ID of each MR fit result
int fqmrnring[50];//                   Number of rings for this fit [1-6]
int fqmrpcflg[50];//                   <0 if MINUIT did not converge during the fit
float fqmrnll[50];//                     Best-fit -lnL
float fqmrtotmu[50];//                   Best-fit total predicted charge
int fqmrpid[50][6];//                  Particle type index for each ring in this fit (Same convention as in 1R fit)
float fqmrmom[50][6];//                  Fit momentum of each ring
float fqmrdconv[50][6];//                Fit conversion length of each ring(always "0" in default mode)
float fqmreloss[50][6];//               Energy lost in the upstream track segment(for upstream tracks only)
float fqmrt0[50][6];//                   Fit creation time of each ring
float fqmrpos[50][6][3];//               Fit vertex position of each ring
float fqmrdir[50][6][3];//               Fit direction of each ring

int fqnse;

int fq1rpcflg[50][7];//                     Flag to indicate whether fiTQun believes the particle is exiting the ID(<0 if MINUIT did not converge)
float fq1rmom[50][7];//                     Fit momentum
float fq1rt0[50][7];//                      Fit particle creation time
float fq1rtotmu[50][7];//                   Best-fit total predicted charge
float fq1rnll[50][7];//                     Best-fit -lnL
float fq1rpos[50][7][3];//                  Fit vertex (0=X, 1=Y, 2=Z)
float fq1rdir[50][7][3];//                  Fit direction (0=X, 1=Y, 2=Z)
float fq1rdconv[50][7];//                   Fit conversion length (always 0 for 1R fits)
float fq1reloss[50][7];//                   Energy lost in the upstream track segment before the hadronic interaction(for upstream tracks only)

//----------fiTQun analysis variables----------//
int mer_id_fiTQun;
float fracmom_2_fiTQun;
float likelihood_fiTQun;
