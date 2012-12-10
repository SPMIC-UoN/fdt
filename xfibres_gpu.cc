/*  xfibres_gpu.cc

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#ifndef EXPOSE_TREACHEROUS
#define EXPOSE_TREACHEROUS
#endif

#include "CUDA/xfibres_gpu.h"
#include "xfibresoptions.h"
#include "newmat.h"
#include "newimage/newimageall.h"

#include <time.h>
#include <sys/time.h>
#include "CUDA/init_gpu.h"


using namespace Xfibres;


void Return_Qform(const Matrix& qform_mat, Matrix& QMat, const float xdim, const float ydim, const float zdim);

//////////////////////////////////////////////////////////
//       XFIBRES CPU PART. IT CALLS TO GPU PART
//////////////////////////////////////////////////////////

int xfibres_gpu(char *subjdir,int slice, int nextras, char** extras)
{
	char slice_str[8];
	char aux[8];
	sprintf(slice_str,"%d",slice);
	while(strlen(slice_str)<4){
		strcpy(aux,"0");
		strcat(aux,slice_str);
		strcpy(slice_str,aux);
	}
	string log_file(subjdir);		//logfile
	log_file += ".bedpostX/logs/log_";
	log_file += slice_str;

	ofstream myfile;
	myfile.open (log_file.data(),ios::out | ios::app );
	myfile << "----------------- SLICE " << slice << " -----------------" << endl;  

	struct timeval t1,t2;
	double time;
    	gettimeofday(&t1,NULL); 

	// Setup logging:
    	Log& logger = LogSingleton::getInstance();
    	xfibresOptions& opts = xfibresOptions::getInstance();
	char **argv=new char*[7+nextras];
    	argv[0] = (char*)malloc(7*sizeof(char));
    	argv[1] = (char*)malloc((40+strlen(subjdir))*sizeof(char));
    	argv[2] = (char*)malloc((40+strlen(subjdir))*sizeof(char));
    	argv[3] = (char*)malloc((40+strlen(subjdir))*sizeof(char));
    	argv[4] = (char*)malloc((40+strlen(subjdir))*sizeof(char));
    	argv[5] = (char*)malloc(10*sizeof(char));
    	argv[6] = (char*)malloc((60+strlen(subjdir))*sizeof(char));
	
    	strcpy (argv[0],"xfibres");    
    	strcpy (argv[1],"--data=");
	strcat (argv[1],subjdir);
	strcat (argv[1],"/data_slice_");
	strcat (argv[1],slice_str);
    	strcpy (argv[2],"--mask=");
	strcat (argv[2],subjdir);
	strcat (argv[2],"/nodif_brain_mask_slice_");
	strcat (argv[2],slice_str);
    	strcpy (argv[3],"--bvals=");
	strcat (argv[3],subjdir);
	strcat (argv[3],"/bvals");
    	strcpy (argv[4],"--bvecs=");
	strcat (argv[4],subjdir);
	strcat (argv[4],"/bvecs");
    	strcpy (argv[5],"--forcedir");
    	strcpy (argv[6],"--logdir=");
    	strcat (argv[6],subjdir);
    	strcat (argv[6],".bedpostX/diff_slices/data_slice_");
    	strcat (argv[6],slice_str);	
	for (int i=0;i<nextras;i++){
		argv[7+i] = (char*)malloc(strlen(extras[i])*sizeof(char));
		strcpy (argv[7+i],extras[i]); 
	}
    	opts.parse_command_line(7+nextras,argv,logger);

	Matrix datam, bvals,bvecs,matrix2volkey;
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
    	matrix2volkey=data.matrix2volkey(mask);
    	vol2matrixkey=data.vol2matrixkey(mask);
    	//Samples samples(vol2matrixkey,matrix2volkey,datam.Ncols(),datam.Nrows());	//in xfibres_gpu.cu

	//Read Gradient Non_linearity Maps if provided
    	NEWIMAGE::volume4D<float> grad; Matrix gradm, Qform, Qform_inv;
    	if (opts.grad_file.set()){
      		read_volume4D(grad,opts.grad_file.value());
      		gradm=grad.matrix(mask);
      		//Get the scale-free Qform matrix to rotate bvecs back to scanner's coordinate system 
      		Return_Qform(data.qform_mat(), Qform, data.xdim(), data.ydim(), data.zdim());
      		Qform_inv=Qform.i();
    	}
	
	myfile << "Number of Voxels: " << datam.Ncols() << endl;  
	myfile << "Number of Directions: " << datam.Nrows() << endl;  

    	if(opts.rician.value() && !opts.nonlin.value()) 
      		cout<<"Rician noise model requested. Non-linear parameter initialization will be performed, overriding other initialization options!"<<endl;

	xfibres_gpu(datam,bvecs,bvals,gradm,Qform,Qform_inv,vol2matrixkey,matrix2volkey,mask,slice,subjdir);

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1); 
	myfile << "Slice processed in: " << time << " seconds" << endl;
	myfile.close();  

  	return 0;
}



//Get the scale-free Qform matrix to rotate bvecs back to scanner's coordinate system 
//After applying this matrix, make sure to sign flip the x coordinate of the resulted bvecs!
//If you apply the inverse of the matrix, sign flip the x coordinate of the input bvecs
void Return_Qform(const Matrix& qform_mat, Matrix& QMat, const float xdim, const float ydim, const float zdim){ 
  	Matrix QMat_tmp; DiagonalMatrix Scale(3); 
  	QMat_tmp=qform_mat.SubMatrix(1,3,1,3);
  	Scale(1)=xdim; Scale(2)=ydim; Scale(3)=zdim;
  	QMat_tmp=Scale.i()*QMat_tmp;
  	QMat=QMat_tmp;
}
