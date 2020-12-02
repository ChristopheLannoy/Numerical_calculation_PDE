/* prototypes utilisés dans main.h */
double mytimer();
int prob(int m, int *n, int **ia, int **ja, double **a, double L,
		 double xratio, double yratio);
int primme(int primme_n, int* primme_ia, int* primme_ja, double* primme_a, 
           int nev, double *evals, double *evecs, int max);

int residu(int nev, int n, double *nrmres, double *evals, double* evecs);

int plot(int m, double *evecs, double L, double xratio, double yratio);

int euler_prog(int iter, double Dt, int n, double **usol);

int plot_temp (int iter, int m, int n, double *usol, double Dt, double L,
			   double xratio, double yratio);


//int dpjd_(int *nt,double *at, int *jat,int *iat,double *eigst,
//	double *rest,double *x,int *lxt,int *neigt,double *sigmat,int *isearcht,
//		int *ninitt,int *madspacet,int *itert,double *tolt,double *shiftt,
//		double *droptolt,double *memt,int *icntlt,int *iprintt,int *infot,
//		double *gapt);

