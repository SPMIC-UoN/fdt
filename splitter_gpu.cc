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
	int nslices=atoi(argv[4]);

	int nextras=argc-5;
	if(gflag==1) nextras++;

	char **extras = new char*[nextras];

	for(int i=5;i<(argc);i++){ //copy all arguments to extras
		extras[i-5] = (char*)malloc(strlen(argv[i])*sizeof(char));
		strcpy(extras[i-5],argv[i]);
	}

	init_gpu();

	char slice_str[8];
	char aux[8];

	for(int slice=0;slice<nslices;slice++){
		sprintf(slice_str,"%d",slice);
		while(strlen(slice_str)<4){
			strcpy(aux,"0");
			strcat(aux,slice_str);
			strcpy(slice_str,aux);
		}

		if(gflag==1){
		    extras[nextras-1] = (char*)malloc((35+strlen(subjdir))*sizeof(char));	
		    strcpy(extras[nextras-1],"--gradnonlin=");
		    strcat(extras[nextras-1],subjdir);
		    strcat(extras[nextras-1],"/grad_dev_slice_");	
		    strcat(extras[nextras-1],slice_str);
		}

		xfibres_gpu(subjdir,slice,nextras,extras);

		printf("SLICE %i processed\n",slice);
	}
	return 0;
}
