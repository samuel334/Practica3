/*********************************************************************************************************
					C R E D I T S

This application uses DISLIN v11.1 plotting library. The DISLIN software is free for non-commercial use.
You can find more information about the project below.

Autor: Helmut Michels
DISLIN Home Page: http://www.dislin.de
/
**********************************************************************************************************/
#include "dislin.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

/*grafica datos obtenidos*/
void graph(char *gtit, char *xlabel, char *ylabel, int size, float max[size], float min[size], float avg[size], int lower, int upper){
	int i;
	float xstep, *x = malloc(sizeof(float)*size);
	for(i = 0;i<size;i++){
	  x[i] = (float)i;
	}
	/*empieza a configurar el histograma*/
	setpag("hp4l");	//page of 2718*1900 points
	metafl("cons");	//metafile format, cons = graphics output on the screen
	scrmod("revers");	//color background
	disini();	//initializes DISLIN
	complx();	//font
	pagera();	//plots a border around the page
	incmrk(1);	//selects line or symbol mode for curve()
	axslen(2300, 1000); //length of the graph
	axspos(300, 1500); //position on the page
	
	name(xlabel, "x");
	name(ylabel, "y");
	titlin(gtit, 2); //plots a title over an axis system

	setgrf("name", "name", "line", "line"); //removes a part of an axis or a complete axis from an axis system
	axsbgd(intrgb(0.95, 0.95, 0.95));	//defines a background colour for axis systems
	polcrv("LINEAR"); //defines an interpolation method used by curve() to connect points

	xstep = ((size<10) ? 1.0 : (float)size/10);	//step between labels in x-axis
	graf(0.0, (float)size - 1.0, 0.0, xstep, 0.0, (float)upper, 0.0, 0.1); //plots a two-dimensional polar axis system
	grid(1, 1); //numbers of grid lines between labels

	thkcrv(3); //defines the thickness of curves

	hsymbl(15);	//defines the size of symbols
	marker(21); //select the symbols used to plot points	
	color("blue"); //defines the colours used for plotting text and lines
	curve(x, min, size); //connects data points with lines or plots them with symbols
	
	hsymbl(15);	//defines the size of symbols
	marker(21); //select the symbols used to plot points
	color("green"); //defines the colours used for plotting text and lines
	curve(x, avg, size); //connects data points with lines or plots them with symbols
	
	hsymbl(15);	//defines the size of symbols
	marker(21); //select the symbols used to plot points
	color("red"); //defines the colours used for plotting text and lines
	curve(x, max, size); //connects data points with lines or plots them with symbols
	
	color("fore");
	height(50);	//defines the character heigh
	title();	//This routine plots a title over an axis system

	endgrf();	//terminates an axis system and sets the level to 1

	disfin();	//terminates DISLIN
}

/*el maximo de los elementos de un arreglo*/
float maxValue(int size, float array[size]){
	int i;
	float max = array[0];
	for(i = 1;i<size;i++){
		if(array[i]>max){
			max = array[i];
		}
	}
	return max;
}

/*el minimo de los elementos de un arreglo*/
float minValue(int size, float array[size]){
	int i;
	float min = array[0];
	for(i = 1;i<size;i++){
		if(array[i]<min){
			min = array[i];
		}
	}
	return min;
}

/*el promedio de los elementos de un arreglo*/
float avgValue(int size, float array[size]){
	int i;
	float avg = 0.0;
	for(i = 0;i<size;i++){
		avg += array[i];
	}
	avg /= (float)size;
	return avg;
}

/*cuadrado de x*/
int square(int x){
	return x*x;
}

/*fenotipo*/
int getPhenotype(int c_lgth, const int chromosome[c_lgth]){
	int i, phenotype = 0;
	for(i = 0;i<c_lgth;i++){
		if(chromosome[i]==1){
			phenotype += (int)pow(2, c_lgth - 1 - i);
		}
	}
	return phenotype;
}

/*genotipo*/
void getChromosome(int c_lgth, int chromosome[c_lgth], int phenotype){
	int i;
	for(i = 0;i<c_lgth;i++){
		chromosome[i] = (phenotype & (int)pow(2, c_lgth - 1 - i))>>(c_lgth - 1 - i);
	}
}

/*seleccion por ruleta*/
void selection(float *max, float *min, float *avg, int p_lgth, int c_lgth, int parent[p_lgth][c_lgth], const int population[p_lgth][c_lgth], int (*f_fitness)(int lgth, const int chromosome[lgth])){
	int i, j, k;
	int *fitness, f_sum = 0; //aptitud y suma de aptitudes
	float r, *ve, cve = 0.0; //valor esperado, suma de valores esperados y valor esperado acumulado
	/*obtener aptitud*/
	fitness = malloc(sizeof(int)*p_lgth);
	ve = malloc(sizeof(float)*p_lgth);
	for(i = 0;i<p_lgth;i++){
		fitness[i] = f_fitness(c_lgth, population[i]);
		f_sum += fitness[i];
	}
	cve = 0.0;
	/*obtener valor esperado(probabilidad)*/
	for(i = 0;i<p_lgth;i++){
		ve[i] = (float)fitness[i]/f_sum;
		cve += ve[i];
	}
	/*seleccionar padres*/
	srand(time(0));
	for(i = 0;i<p_lgth;i++){
		r = (float)((double)rand()/(double)RAND_MAX);
		cve = 0.0;
		for(j = 0;j<p_lgth;j++){
			cve += ve[j];
			if(cve>r){
				for(k = 0;k<c_lgth;k++){
					parent[i][k] = population[j][k];
				}
				break;
			}
		}
	}
	/*obtener valor esperado maximo de la generacion actual*/
	*max = maxValue(p_lgth, ve);
	/*obtener valor esperado minimo de la generacion actual*/
	*min = minValue(p_lgth, ve);
	/*obtener valor esperado promedio de la generacion actual*/
	*avg = avgValue(p_lgth, ve);
	free(fitness);
	free(ve);
}

/*cruza dos individuos y obtiene un descendiente*/
int crossover(int c_lgth, int child[c_lgth], const int parent1[c_lgth], const int parent2[c_lgth], int cp){
	int i;
	for(i = 0;i<c_lgth;i++){
		if(i<cp){
			child[i] = parent1[i];
		}else{
			child[i] = parent2[i];
		}
	}
}

/*cruzar poblacion*/
void breeding(int p_lgth, int c_lgth, int child[p_lgth][c_lgth], const int parent[p_lgth][c_lgth]){
	int i = 0, k, lower = 1, upper = c_lgth - 2;
	int cp; //punto de cruza
	srand(time(0));
	while(i<p_lgth){
		cp = (rand() % (upper - lower + 1)) + lower;
		crossover(c_lgth, child[i], parent[i], parent[i+1], cp);
		crossover(c_lgth, child[i+1], parent[i+1], parent[i], cp);
		i += 2;
	}
}

/*mutacion*/
void mutation(int p_lgth, int c_lgth, int population[p_lgth][c_lgth], float pp, float bp){
	int i, j;
	srand(time(0));
	for(i = 0;i<p_lgth;i++){
		if((float)((double)rand()/(double)RAND_MAX)<=pp){
			for(j = 0;j<c_lgth;j++){
				if((float)((double)rand()/(double)RAND_MAX)<=bp){
					population[i][j] = ((population[i][j]==1) ? 0 : 1);
				}
			}
		}
	}
}

/*funcion de aptitud*/
int getFitness(int c_lgth, const int chromosome[c_lgth]){
	return square(getPhenotype(c_lgth, chromosome));
}

int main(int argc, char *argv[]){
	if(argc>2 || argc<2){
		exit(1);
	}
	int gen = atoi(argv[1]); //numero de generaciones
	int i, j, k, lower = 0, upper = 31;
	int p_lgth = 31; //numero de individuos
	int c_lgth = 5;	//tamaño de la cadena binaria
	int phenotype[p_lgth];	//individuos
	int chromosome[p_lgth][c_lgth];	//cromosomas
	int parent[p_lgth][c_lgth], child[p_lgth][c_lgth];	//padres y descendencia
	float *max, *min, *avg;
	/*generar poblacion inicial*/
	srand(time(0));
	for(i = 0;i<p_lgth;i++){
		phenotype[i] = (rand() % (upper - lower + 1)) + lower;
		getChromosome(c_lgth, chromosome[i], phenotype[i]);
	}
	max = malloc(sizeof(float)*gen);
	min = malloc(sizeof(float)*gen);
	avg = malloc(sizeof(float)*gen);
	/*algoritmo genetico simple*/
	printf("%d generaciones\n", gen);
	for(i = 0;i<gen;i++){
		/*seleccion*/
		selection(&max[i], &min[i], &avg[i], p_lgth, c_lgth, parent, chromosome, getFitness);
		/*cruza*/
		breeding(p_lgth, c_lgth, child, parent);
		/*mutar*/
		mutation(p_lgth, c_lgth, child, 0.02, 0.2);
		/*insertar en nueva generacion*/
		for(j = 0;j<p_lgth;j++){
			for(k = 0;k<c_lgth;k++){
				chromosome[j][k] = child[j][k];
			}
			phenotype[j] = getPhenotype(c_lgth, chromosome[j]);
		}
	}
	/*imprimir maximos, minimos y promedios resultantes*/
	for(i = 0;i<gen;i++){
		printf("min[%d]: %f, max[%d]:%f, avg[%d]: %f\n", i, min[i], i, max[i], i, avg[i]);
	}
	/*graficar histograma de maximos, minimos y promedios*/
	graph("Histograma", "Generaciones", "Funcion Aptitud", gen, max, min, avg, 0, 1);
	free(max);
	free(min);
	free(avg);
	return 1;
}
