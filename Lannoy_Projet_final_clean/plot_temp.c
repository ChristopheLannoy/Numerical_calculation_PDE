#include <stdio.h>
#include <math.h>

int plot_temp (int iter, int m, int n, double *usol,double Dt,double L,
			   double xratio, double yratio){
    
    double h = L/(m-1);
	int ix, iy, ind, it;
	
	int nx = m-2;
	
	int ixlim = ((m-1) * xratio) - 1;
    int iylim = ((m-1) * yratio) - 1;   /* coordonnées du point caractérisant la membrane */
    
    int npts; //nombre de point sautés (dans l'ordre lexicographique), c-a-d qui n'est
              //pas représenté dans notre probléme car on sait que sa température vaut 0
              //point dont l'équation n'a donc pas été rprésentée dans la matrice A
    
    
    
  //on veut afficher 100 graphiques peu importe le nombre d'iteration de euler prog
  //on afichera donc sur gnuplot que les vecteurs correspondant a t = i*saut*Dt pour i allant de 0 à 100
  //Dt correspondant au pas de temps de la résolution par euler progressif
    
    int saut = iter/100;         
    int nombregraph = iter/saut; //approximativement 100 (division entière)
    
    FILE *data = fopen("dataTemp.txt", "w"); //écriture de fichier


    for(it = 0; it <= iter; it+= saut){
		npts = 0;
		for (iy = -1; iy <= nx; iy++) {
			for (ix = -1; ix <= nx; ix++) {
				
				
	/* Points appartenant au bord du carré */
				if(ix == -1 || iy == -1 || ix == nx || iy == nx){
					fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 0.0);
					continue;
				}
	/* Reste les points appartenant à l'intérieur du carré */
            
			/* Points appartenant à l'intérieur du carré mais pas à l'intérieur de la membrane */
				if(iy >= iylim && ix >= ixlim){ 
					npts++;
					fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , 0.0);
					continue;
				}
			/*  Point intérieur à la membrane (donc seuls points potentiellement non nuls) */
				ind = ix + nx * iy -npts;   /* numéro de l'éauqtion */
			
				fprintf(data, "%f %f %f\n",(ix+1)*h, (iy+1)*h , usol[it*n + ind]);
            
	
			}
			//fprintf(data, "\n");
		}
		fprintf(data, "\n \n");
	}
	fclose(data);
	
	
	
	

	// fichier commande
	FILE *cmd = fopen("cmdTemp.txt", "w");
	fprintf(cmd,
		"unset key\n"
		"set cblabel 'Temperature (en Kelvin)'\n"
		
		"set cbrange [0:10] \n"
		
		
		"do for[i=0: %d :1]{\n"
		"	set multiplot \n"
		
		
		"		set size 0.1,1\n"
		"		set origin 0,0\n"
		
		
		"		set xrange [0:1]\n"
		"		set yrange [0:105]\n"
		"		set title '                            visualisation du temps (n ieme image)' \n"
		"		plot i \n"
		
		
		"		set size 0.9,0.9\n"
		"		set origin 0.1,0\n"
		
		
		"		set xrange [0: %d ]\n"
		"		set yrange [0: %d ]\n"
		
		"		set title 'Graphique de la distribution de temperature' \n" 
		"		plot 'dataTemp.txt' index i using 1:2:3 with image \n"
		
		"	unset multiplot \n"
		//"	pause 0.1 \n"     //nécessaire si petite résolution (m petit) (ou ordi très puissant :) )
		"}", nombregraph,(int)L,(int)L);
	fclose(cmd);

	//exécution
	system("gnuplot -persist cmdTemp.txt");
	
	return 0;
}
