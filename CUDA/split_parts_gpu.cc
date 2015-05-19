/*  split_parts_gpu.cc

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "newimage/newimageall.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>

void save_part(Matrix data, string name, int idpart){

	int nvox = data.Ncols();
	int ndir = data.Nrows();

	string file_name;
	file_name = name+"_"+num2str(idpart);

	ofstream out;
	out.open(file_name.data(), ios::out | ios::binary);
	out.write((char*)&data(1,1),nvox*ndir*sizeof(Real));
	out.close();
}


// parameters:
// 1. data.nii.gz
// 2. mask.nii.gz
// 3. grad_dev.nii.gz
// 4. gflag
// 5. nparts
// 6. output directory


int main(int argc, char *argv[]){

	std::string data_str = argv[1];
	std::string mask_str = argv[2];
	std::string grad_str = argv[3];
	int gflag = atoi(argv[4]);
	int nparts = atoi(argv[5]);
	std::string out_dir = argv[6];

	//printf("%s\n%s\n%s\n%i\n%i\n%s\n",data_str.data(),mask_str.data(),grad_str.data(),gflag,nparts,out_dir.data());
	NEWIMAGE::volume4D<float> data;
	NEWIMAGE::volume<float> mask;
    	read_volume4D(data,data_str);
    	read_volume(mask,mask_str);
	Matrix datam;
    	datam=data.matrix(mask); 

	int nvoxels=datam.Ncols();
	int ndirections=datam.Nrows();
	
	NEWIMAGE::volume4D<float> grad; 
	Matrix gradm;
	int dirs_grad=0;
	if(gflag){
		read_volume4D(grad,grad_str);
      		gradm=grad.matrix(mask);
		dirs_grad = gradm.Nrows();
	}
	

	int size_part=nvoxels/nparts;

	Matrix data_part;
	Matrix grad_part;
	string out_data;
	string out_grad;
	out_data.append(out_dir);
	out_grad.append(out_dir);
	out_data.append("/data");
	out_grad.append("/grad_dev");
	for(int i=0;i<(nparts-1);i++){
		data_part = datam.SubMatrix(1,ndirections,i*size_part+1,(i+1)*size_part);
		save_part(data_part,out_data,i);
		if(gflag){
			grad_part = gradm.SubMatrix(1,dirs_grad,i*size_part+1,(i+1)*size_part);
			save_part(grad_part,out_grad,i);
		}
	}
	// last part
	data_part = datam.SubMatrix(1,ndirections,(nparts-1)*size_part+1,nvoxels);
	save_part(data_part,out_data,(nparts-1));
	if(gflag){
		grad_part = gradm.SubMatrix(1,dirs_grad,(nparts-1)*size_part+1,nvoxels);
		save_part(grad_part,out_grad,(nparts-1));
	}
}
