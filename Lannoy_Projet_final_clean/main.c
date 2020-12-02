#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#define D 9.7e-5

int main(int argc, char *argv[])
{
  /* déclarer les variables */
  int jj,m =150, nev = 1; 
  m=(m/8)*8+1; /*m<150 prob1  m<140 prob2 (limitation memoire)  m = x.8 + 1*/ 

  double L = 2, xratio = 5.0/8, yratio = 2.0/8;  //(prob1)
//double L = 1, xratio = 4.0/8, yratio = 4.0/8; //génère un autre problème (prob2)
     // intéressant de réduire tmax à 1000 car 
     //refroidissement plus rapide d'une membrane plus petite
//double L = 9, xratio = 1.0/8, yratio = 1.0/8; //(prob3) (m=200)
  
// si l'affichage de l'évolution de la température ne fonctionne pas
// diminuer m au valeurs recommendées (pour mon pc en particlier)  
  
  int n, *ia, *ja; 
  double *a;
  double *evals, *evecs;
  double t1, t2;
  double nrmres;

  /* générér le problème */
  if (prob(m, &n, &ia, &ja, &a, L, xratio, yratio))
     return 1;

  printf("\nPROBLEM: ");
  printf("m = %5d   n = %8d  nnz = %9d\n", m, n, ia[n] );
  
  /* affichage de la matrice A */
  
  /*
   
  printf("\n \n");
  for (jj=0; jj<=n; jj++){
            printf("%d \n ", ia[jj]);
  }
  printf("\n \n");
   for (jj=0; jj<ia[n]; jj++){
            printf("%d \n ", ja[jj]);
  }
  printf("\n \n");
   for (jj=0; jj<ia[n]; jj++){
            printf("%f \n ", a[jj]);
  }
   
  */
  
  
  /* allouer la memoire pour vecteurs & valeurs propres */
  evals = malloc(nev * sizeof(double));
  evecs = malloc(nev * n * sizeof(double));

  if (evals == NULL || evecs == NULL) {
      printf("\n ERREUR : pas assez de mémoire pour les vecteurs et valeurs propres\n\n");
      return 1;
  }

  /* primme - résolution */
  t1 = mytimer();
  if(primme(n, ia, ja, a, nev, evals, evecs, 0)) // 0= eigen pair for eigen value min , 1=>max
     return 1;
  t2 = mytimer();

  /* temps de solution */
  printf("\nTemps de solution (CPU) (PRIMME): %5.1f sec\n",t2-t1);
  printf("\nValeur propre minimale calculée (PRIMME): %5.1f\n",evals[0]);
  
  plot(m,evecs,L,xratio,yratio);
  
  
 //////////////////////////////////////////////////////////////////////
 ////////////      Q2  Norme résiduelle                  //////////////
 //////////////////////////////////////////////////////////////////////

  printf("\n\nQ2 Norme résiduelle \n\n");
  if(residu(nev, n, &nrmres, evals, evecs) ) 
	 return 1;
  printf("Norme résiduelle selon ma fct %e \n", nrmres);   
  
 
 
 
 
 //////////////////////////////////////////////////////////////////////
 ////////////    Q5 solveur alternatif:  JADAMILU        //////////////
 //////////////////////////////////////////////////////////////////////
printf("\n\n Q5 Solveur alternatif : JADAMILU");
// Conversion du format CSR en C (1) vers un format CSR en fortran de la 
// partie triangulaire supérieur (2).
// On suppose la matrice symétrique 
///////////////////////////////////////////////////////////////////////


int n1= n;                  //taille de la matrice A (1)
int nnz1 = 5*n - 4*(m-2);   //nombre d'elements non nul dans la matrice A (1)
double *a1 = a;
int *ja1 = ja;
int *ia1 = ia;    
 
 
// Exemple de matrice simple pour vérifier le fonctionnement de la conversion
//int n1= 3;     //taille de la matrice 1
//int nnz1 = 5;  //nombre d'elements non nul dans la matrice 1
//double a1[5] = {1,1,2,1,3}; //  1 0 1 
//int ja1[5] = {0,2,1,0,2};   //  0 2 0
//int ia1[4] = {0,2,3,5};     //  1 0 3



// définition de la taille de la matrice 2 //
int *ia2;
int nelem_Ja_A2= (nnz1-n1)/2 + n1; //taille des matrices ja2 et a2
int *ja2;
double *a2;
int n2 = n1;

ia2 = malloc( (n1+1)*sizeof(int));
ja2 = malloc( (nelem_Ja_A2)*sizeof(int));
a2 = malloc( (nelem_Ja_A2)*sizeof(double));

if (ia2 == NULL || ja2 == NULL || a2 == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer la matrice 2\n\n");
        return 1;
    }


//Remplissage de la matrice 2//
int i,j,nelem2=0,nelem1=0;
for(i=0;i<n1;i++){
	ia2[i]=nelem2+1;
	for(j=ia1[i];j<ia1[i+1];j++){
		if(ja1[j]>=i){      // Si on est dans la partie triangulaire sup
			a2[nelem2]= a1[nelem1];
			ja2[nelem2]=ja1[j] +1;
			nelem2 ++; nelem1++;
		}
		else{nelem1++;}
	}
}
ia2[n1]=nelem2+1; //dernière élément de ia2


//Affichage de la matrice 2 convertie
/*
int r;
for(r=0; r<=n1;r++)
	printf("ia2[%d]= %d \n",r,ia2[r]);

for(r=0; r<nelem_Ja_A2;r++){
	printf("ja2[%d]= %d \n",r,ja2[r]);
}
for(r=0; r<nelem_Ja_A2;r++){
	printf("a2[%d]= %e \n",r,a2[r]);
}
*/






////   Résolution avec jadamilu après conversion de la matrice   ///////
////////////////////////////////////////////////////////////////////////



////   Test du solveur jadamilu avec une matrice simple (débugage)   ////
//matrice test sous format CSR (triangulaire supérieure) (commenter la transformation de matrice)
//int ia2[6] = {1,3,5,7,9,10};
//int ja2[9] = {1,2,2,3,3,4,4,5,5};
//double a2[9] = {0.0,5.0,1.0,5.0,2.0,5.0,3.0,5.0,4.0};
//int n2 = 5;


int maxeigt=1,maxspt= 20;  
int lxt = n2*(3*maxspt+maxeigt+1) + 4*maxspt*maxeigt;

double eigst[maxeigt];
double rest[maxeigt];  
double xt[lxt]; 
int neigt, madspacet, isearcht, ninitt;
int icntlt[5];
int itert,iprintt, infot;
double sigmat, tolt, shiftt, gapt, memt, droptolt;

iprintt = 6;
isearcht = 0; //sm eigenvalues
memt = 20.0; 
droptolt = 1.0e-3; 
neigt = maxeigt;
ninitt = 0; // pas d'estimation initiale
madspacet = maxspt; 
itert = 1000; 
tolt = 1.0e-10;
icntlt[0] = 0; 
icntlt[1] = 0; 
icntlt[2] = 0; 
icntlt[3] = 0; 
icntlt[4] = 0; 

sigmat = 0;
shiftt = 0;


t1 = mytimer();
dpjd_(&n2,a2,ja2,ia2,eigst,rest,xt,&lxt,&neigt,
		&sigmat,&isearcht,&ninitt,&madspacet,&itert,&tolt,
		&shiftt,&droptolt,&memt,icntlt,&iprintt,&infot,
		&gapt); 
t2 = mytimer();
printf("\nTemps de solution (CPU) JADAMILU: %5.1f sec\n",t2-t1);
printf("Valeur propre minimale JADAMILU:: %5.1f\n\n\n",eigst[0]);
  
  
 //////////////////////////////////////////////////////////////////////
 ////////////        Question 4 Euler progressif         //////////////
 //////////////////////////////////////////////////////////////////////
 
 printf("\n\nQ4 Euler progressif: Détermination du pas de temps\n\n");
 //calcule la valeur propre maximale pour déterminer le pas de temps
 if(primme(n, ia, ja, a, nev, evals, evecs, 1))   // 0= eigen pair for eigen value min , 1=>max
     return 1;
 
 printf("Valeur propre maximale calculée: %5.1f \n",evals[0]);
 
 double Dt = 2.0/(evals[0]*D);    // pas de temps de résolution pour euler progressif

// Dt *= 1.05; //montre l'instabilité de eulerp si Dt trop grand
 
 printf("Pas de temps correspondant pour stabilité (Dt): %f\n",Dt);
 
 if(Dt>10)
	Dt = 10;   // pas de temps minimal pour pouvoir afficher 100 images
 
 double tmax = 2000; //arret de la simulation
 int iter = tmax/Dt;  //nombre d'iterations nécessaire à eulerp our atteindre tmax
 
 double *usol; //stockage de la solution de eulerp
 usol = malloc((iter+1)*n*sizeof(double));
 if (usol == NULL) {
      printf("\n ERREUR : pas assez de mémoire pour la matrice usol (via euler)\n\n");
      return 1;
 }
 
 if(euler_prog(iter, Dt, n, &usol)){
	 printf("Erreur lors de la résolution par Euler \n");
	 return 1;
 } 
 
 plot_temp(iter, m, n, usol, Dt, L, xratio, yratio);
  

  /* libérér la mémoire */
  free(usol);
  free(ia); free(ja); free(a); free(evals); free(evecs);
  return 0;
}

