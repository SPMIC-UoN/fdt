/*  xfibres_gpu.cc

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "CUDA/xfibres_gpu.h"
#include "xfibresoptions.h"
#include "newmat.h"
#include "newimage/newimageall.h"

#include <time.h>
#include <sys/time.h>
#include "CUDA/init_gpu.h"

#define SIZE_SUB_PART 9000 

using namespace Xfibres;

//////////////////////////////////////////////////////////
//       XFIBRES CPU PART. IT CALLS TO GPU PART
//////////////////////////////////////////////////////////
// last 2 parameters are idPart and nParts

int main(int argc, char *argv[]){

	struct timeval t1,t2;
	double time;
    	gettimeofday(&t1,NULL); 

	// Setup logging:
    	Log& logger = LogSingleton::getInstance();
    	xfibresOptions& opts = xfibresOptions::getInstance();
	opts.parse_command_line(argc-2,argv,logger);

	Matrix datam, bvals,bvecs;
    	NEWIMAGE::volume<float> mask;
    	NEWIMAGE::volume<int> vol2matrixkey;
    	bvals=read_ascii_matrix(opts.bvalsfile.value());
    	bvecs=read_ascii_matrix(opts.bvecsfile.value());
    	if(bvecs.Nrows()>3) bvecs=bvecs.t();
    	if(bvals.Nrows()>1) bvals=bvals.t();
    	for(int i=1;i<=bvecs.Ncols();i++){
      		float tmpsum=sqrt(bvecs(1,i)*bvecs(1,i)+bvecs(2,i)*bvecs(2,i)+bvecs(3,i)*bvecs(3,i));
      		if(tmpsum!=0){
			bvecs(1,i)=bvecs(1,i)/tmpsum;
			bvecs(2,i)=bvecs(2,i)/tmpsum;
			bvecs(3,i)=bvecs(3,i)/tmpsum;
      		}  
    	}

    	NEWIMAGE::volume4D<float> data;
    	read_volume4D(data,opts.datafile.value());
    	read_volume(mask,opts.maskfile.value());
    	datam=data.matrix(mask); 
	
	///////////////////////////////////////////
	////////// Read my part of data ///////////
	///////////////////////////////////////////
	int idPart = atoi(argv[argc-2]);
	int nParts = atoi(argv[argc-1]);
	int ndirections = bvals.Ncols();
	int dirs_grad = 0;
	int totalNvox = datam.Ncols();

	int size_part = totalNvox / nParts;
	Matrix mydatam;
	Matrix mygradm;
	if(idPart!=(nParts-1)){
		mydatam = datam.SubMatrix(1,ndirections,idPart*size_part+1,(idPart+1)*size_part);
	}else{
		mydatam = datam.SubMatrix(1,ndirections,idPart*size_part+1,totalNvox);
	}
	
	//Read Gradient Non_linearity Maps if provided
    	NEWIMAGE::volume4D<float> grad; Matrix gradm;
    	if (opts.grad_file.set()){
      		read_volume4D(grad,opts.grad_file.value());
      		gradm=grad.matrix(mask);
		dirs_grad = gradm.Nrows();
		if(idPart!=(nParts-1)){
			mygradm = gradm.SubMatrix(1,dirs_grad,idPart*size_part+1,(idPart+1)*size_part);
		}else{
			mygradm = gradm.SubMatrix(1,dirs_grad,idPart*size_part+1,totalNvox);
		}
    	}

	if(idPart==(nParts-1)) size_part = totalNvox - ((nParts-1)*size_part);
	 
	cout << "Number of Voxels to compute in this part: " << size_part << endl;  
	cout << "Number of Directions: " << ndirections << endl;  

    	if(opts.rician.value() && !opts.nonlin.value()) 
      	cout<<"Rician noise model requested. Non-linear parameter initialization will be performed, overriding other initialization options!"<<endl;

	//////////////////////////////////////////////////////////////
	////////// Divide the process in Subparts ////////////////////
	//////////////////////////////////////////////////////////////

	int nsubparts = (size_part/SIZE_SUB_PART); 			//SIZE_SUB_PART voxels per iteration
	int size_sub_part = SIZE_SUB_PART;
	if(size_part%SIZE_SUB_PART) nsubparts++;
	int last_sub_part = size_part - ((nsubparts-1)*SIZE_SUB_PART);

	if(last_sub_part<(SIZE_SUB_PART*0.5)){ 	//if is too small the last part we distribute it between the others
		size_sub_part = size_sub_part + last_sub_part/(nsubparts-1);
		nsubparts--;
		last_sub_part = size_part - ((nsubparts-1)*size_sub_part);		
	}

	Matrix mydatam_part;	
	Matrix mygradm_part;	
	
	for(int i=0;i<nsubparts-1;i++){
		cout << endl << "///////////////////////////////////////////////////////" << endl;
		cout <<         "/////////////////  SubPart " << i+1 << " of  " << nsubparts << "/////////////////////" << endl;
		cout <<  "///////////////////////////////////////////////////////" << endl;
		cout <<  "Processing " << size_sub_part << " voxels" << endl;
		mydatam_part = mydatam.SubMatrix(1,ndirections,i*size_sub_part+1,(i+1)*size_sub_part);
		if (opts.grad_file.set()) mygradm_part = mygradm.SubMatrix(1,dirs_grad,i*size_sub_part+1,(i+1)*size_sub_part);
		xfibres_gpu(mydatam_part,bvecs,bvals,mygradm_part,i);
		cout << endl;
	}

	cout << endl << "///////////////////////////////////////////////////////" << endl;
	cout <<         "/////////////////  SubPart " << nsubparts << " of  " << nsubparts << "/////////////////////" << endl;
	cout <<  "///////////////////////////////////////////////////////" << endl;
	cout <<  "Processing " << last_sub_part << " voxels" << endl;
	mydatam_part = mydatam.SubMatrix(1,ndirections,(nsubparts-1)*size_sub_part+1,size_part);
	if (opts.grad_file.set()) mygradm_part = mygradm.SubMatrix(1,dirs_grad,(nsubparts-1)*size_sub_part+1,size_part);
	xfibres_gpu(mydatam_part,bvecs,bvals,mygradm_part,nsubparts-1);
	cout << endl;

	//////////////////////////////////////////////////////////////
	////////// JOIN Results of the Subparts //////////////////////
	//////////////////////////////////////////////////////////////

	if(opts.modelnum.value()==1){
		join_subParts("mean_dsamples",size_part,nsubparts,size_sub_part,last_sub_part,true);
	}else if(opts.modelnum.value()==2){
		join_subParts("mean_dsamples",size_part,nsubparts,size_sub_part,last_sub_part,true);
		join_subParts("mean_d_stdsamples",size_part,nsubparts,size_sub_part,last_sub_part,true);
		//join_subParts("dsamples",size_part,nsubparts,size_sub_part,last_sub_part,false);
		//join_subParts("d_stdsamples",size_part,nsubparts,size_sub_part,last_sub_part,false);
	}	
	if (opts.f0.value()){
      		join_subParts("mean_f0samples",size_part,nsubparts,size_sub_part,last_sub_part,true);
		//join_subParts("f0samples",size_part,nsubparts,size_sub_part,last_sub_part,false);
    	}
    	if (opts.rician.value()){
		join_subParts("mean_tausamples",size_part,nsubparts,size_sub_part,last_sub_part,true);
    	}

	for(int f=0;f<opts.nfibres.value();f++){
		join_subParts("th"+num2str(f+1)+"samples",size_part,nsubparts,size_sub_part,last_sub_part,false);
		join_subParts("ph"+num2str(f+1)+"samples",size_part,nsubparts,size_sub_part,last_sub_part,false);
		join_subParts("f"+num2str(f+1)+"samples",size_part,nsubparts,size_sub_part,last_sub_part,false);

		//join_subParts("mean_f"+num2str(f+1)+"samples",size_part,nsubparts,size_sub_part,last_sub_part,true);
		//join_subParts("dyads"+num2str(f+1),size_part,nsubparts,size_sub_part,last_sub_part,false);
	}	

	join_subParts("mean_S0samples",size_part,nsubparts,size_sub_part,last_sub_part,true); //the last one to control with monitor
		
	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1); 
	cout << "Part processed in: " << time << " seconds" << endl;

  	return 0;
}

void join_subParts(string name, int size_part, int nsubparts, int size_sub_part, int last_sub_part, bool mean){
	Log& logger = LogSingleton::getInstance();
    	xfibresOptions& opts = xfibresOptions::getInstance();

	int nsamples = opts.njumps.value()/opts.sampleevery.value();
	if(mean) nsamples=1;

	Matrix tmp(nsamples,0);
	Matrix part; 
	ifstream in;
	ofstream out;

	string file_name;
		
	for(int i=0;i<nsubparts-1;i++){
		part.ReSize(nsamples,size_sub_part);
		file_name = logger.appendDir(name+"_"+num2str(i));
		in.open(file_name.data(), ios::in | ios::binary);
		in.read((char*)&part(1,1), size_sub_part*nsamples*sizeof(Real));
		in.close();
		remove (file_name.data());
		tmp = tmp | part;
	}
	part.ReSize(nsamples,last_sub_part);
	file_name = logger.appendDir(name+"_"+num2str(nsubparts-1));
	in.open (file_name.data(), ios::in | ios::binary);
	in.read((char*)&part(1,1), last_sub_part*nsamples*sizeof(Real));
	in.close();
	remove (file_name.data());
	tmp = tmp | part;

	file_name = logger.appendDir(name+"J");
	out.open(file_name.data(), ios::out | ios::binary);
	out.write((char*)&tmp(1,1),size_part*nsamples*sizeof(Real));
	out.close();

}
