#include <math.h>
#include <stdio.h>
#include <stdlib.h>


int euler_prog(int iter, double Dt, int n, double **usol){
	
	int i, it;
	double D = 9.7e-5;
	
	
	// temp = A*u de l'itération précédante
	double *temp;      
	temp = malloc(n*sizeof(double));
	if (temp == NULL){
		printf("\n ERREUR : pas assez de mémoire pour temp (euler)\n\n");
		return 1;
	}
	
	int nev = 1;
	
	// condition initiale
	double T0 = 10;
	
	for (i = 0; i<n; i++){
		(*usol)[i] = T0;
	}
	
	
	// resolution par Euler progressif
	
	for (it = 1; it<=iter; it++){
		
		matvec_primme( &((*usol)[(it-1)*n]), temp, &nev);  //pointeur
		
		for(i = 0; i<n; i++){
			(*usol)[it*n + i] = (*usol)[(it-1)*n + i] - Dt*(D)*temp[i];
		}
		
		//printf("(*usol)[it*n + 9] %f \n", (*usol)[(it-1)*n + 9]);
		
	}
	free(temp);
	return 0;
	
	
	
	
}
