#include <stdlib.h>
#include <stdio.h>
#include <math.h>


int residu(int nev, int n, double *nrmres, double *evals, double *evecs)

{
	double *vy, *vz;
	int i;
	vy = malloc( nev * n * sizeof(double) );
	vz = malloc( nev * n * sizeof(double) );
	
    if (vy == NULL || vz == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice (residu)\n\n");
        return 1;
    }
	
	
	/*  vy =  A.evecs  */
	matvec_primme(evecs, vy, &nev);
	
	/* vz = vy - evals . evecs */
	
	for (i=0; i < n; i++) {
		vz[i] = vy[i] - evals[0] * evecs[i];	
	}

	/*   norme 2 de vz sur norme 2 de evecs*/
	double A = 0.0, B = 0.0;
	
	for (i=0; i < n; i++) {
		A+= vz[i] * vz[i];	
	}
	for (i=0; i < n; i++) {
		B+= evecs[i] * evecs[i];	
	}
	
	*nrmres = sqrt(A)/sqrt(B);
	
	free(vy);
	free(vz);  
	
	return 0;
}
