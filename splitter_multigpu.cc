/*  splitter_multigpu.cc

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <sys/time.h>
#include "CUDA/xfibres_gpu.h"
#include "CUDA/init_gpu.h"

using namespace std;

int main(int argc, char *argv[]){


	char *subjdir=argv[1];
	int gflag=atoi(argv[2]);
	int nfibres=atoi(argv[3]);
	int slice=atoi(argv[4]);

	printf("\n-------------------------------\n");
	printf("SLICE %i\n",slice);

	int nextras=argc-5;
	if(gflag==1) nextras++;

	char **extras = new char*[nextras];

	for(int i=5;i<(argc);i++){ //copy all arguments to extras
		extras[i-5] = (char*)malloc(strlen(argv[i])*sizeof(char));	
		strcpy(extras[i-5],argv[i]);
	}

	init_gpu(); //exclusive mode ?

	char slice_str[8];
	char aux[8];


	sprintf(slice_str,"%d",slice);
	while(strlen(slice_str)<4){
		strcpy(aux,"0");
		strcat(aux,slice_str);
		strcpy(slice_str,aux);
	}

	if(gflag==1){
		    extras[nextras] = (char*)malloc((35+strlen(subjdir))*sizeof(char));	
		    strcpy(extras[nextras],"--gradnonlin=");
		    strcat(extras[nextras],subjdir);
		    strcat(extras[nextras],"/grad_dev_slice_");	
		    strcat(extras[nextras],slice_str);
	}

	xfibres_gpu(subjdir,slice,nextras,extras);

	printf("--------------------------------\n");
	printf("SLICE %i processed\n",slice);
	printf("--------------------------------\n");	

	return 0;
}
