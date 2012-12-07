#include "newmat.h"
#include "newimage/newimageall.h"
#include "xfibresoptions.h"

#include "xfibres_gpu.cuh"
#include "diffmodels.cuh"
#include "runmcmc.h"
#include "samples.h"
#include "options.h"

#include <host_vector.h>
#include <device_vector.h> 

#include <time.h>
#include <sys/time.h>
#include "init_gpu.h"
#include <fstream>

using namespace Xfibres;

void xfibres_gpu(	//INPUT
			const Matrix			datam,
			const Matrix			bvecs,
			const Matrix			bvals,
			const Matrix	 		gradm, 
			const Matrix 			Qform, 
			const Matrix 			Qform_inv,
			const NEWIMAGE::volume<int> 	vol2matrixkey,
			const NEWMAT::Matrix		matrix2volkey,
			const NEWIMAGE::volume<float>	mask,
			const int 			slice,
			const char*			subjdir)
{
	//write num of slice in a string for log file
	char slice_str[8];
	char aux[8];
	sprintf(slice_str,"%d",slice);
	while(strlen(slice_str)<4){
		strcpy(aux,"0");
		strcat(aux,slice_str);
		strcpy(slice_str,aux);
	}
	string gpu_log(subjdir);		//logfile
	gpu_log += ".bedpostX/logs/times_gpu_";
	gpu_log += slice_str;

	xfibresOptions& opts = xfibresOptions::getInstance();

	///// FIT /////
	thrust::host_vector<double> datam_host, bvecs_host, bvals_host, alpha_host, beta_host, params_host;
	thrust::host_vector<float> tau_host;
	vector<ColumnVector> datam_vec;
	vector<Matrix> bvecs_vec, bvals_vec;

	prepare_data_gpu_FIT(datam,bvecs,bvals,gradm,Qform,Qform_inv,datam_vec, bvecs_vec, bvals_vec, datam_host, bvecs_host,  bvals_host, alpha_host, beta_host, params_host, tau_host);	

	int nvox = datam.Ncols();
	int nfib= opts.nfibres.value(); 

	thrust::device_vector<double> datam_gpu=datam_host;
	thrust::device_vector<double> bvecs_gpu=bvecs_host;
	thrust::device_vector<double> bvals_gpu=bvals_host;	
	thrust::device_vector<double> params_gpu=params_host;
	thrust::host_vector<int> vox_repeat;	//contains the id's of voxels repeated
	vox_repeat.resize(nvox);
	int nrepeat=0;

	fit(datam_vec,bvecs_vec,bvals_vec,datam_host,bvecs_host,bvals_host,datam_gpu,bvecs_gpu,bvals_gpu,gpu_log,params_gpu,vox_repeat,nrepeat);

	if(opts.rician.value()){
		calculate_tau(datam_gpu,params_gpu,bvecs_gpu,bvals_gpu,vox_repeat,nrepeat,gpu_log, tau_host);
	}

	bvecs_gpu.clear();		//free bvecs_gpu
	bvecs_gpu.shrink_to_fit();
	
	//////   RUN MCMC  //////
	thrust::host_vector<double> signals_host,isosignals_host;
	thrust::host_vector<FibreGPU> fibres_host;
	thrust::host_vector<MultifibreGPU> multifibres_host;
	
	prepare_data_gpu_MCMC(nvox, nfib, signals_host, isosignals_host, fibres_host, multifibres_host);

	thrust::device_vector<double> signals_gpu=signals_host;
	thrust::device_vector<double> isosignals_gpu=isosignals_host;
	thrust::device_vector<FibreGPU> fibres_gpu=fibres_host;
	thrust::device_vector<MultifibreGPU> multifibres_gpu=multifibres_host;
	thrust::device_vector<float> tau_gpu = tau_host;
	thrust::device_vector<double> alpha_gpu=alpha_host;
	thrust::device_vector<double> beta_gpu=beta_host;

	init_Fibres_Multifibres(datam_gpu, params_gpu, tau_gpu, bvals_gpu, alpha_gpu, beta_gpu, gpu_log, fibres_gpu, multifibres_gpu, signals_gpu, isosignals_gpu);

	srand(opts.seed.value());
	//double seed1 = rand(); 
	//double seed2 = rand(); 

	runmcmc_burnin(datam_gpu, bvals_gpu, alpha_gpu, beta_gpu, rand(), gpu_log, fibres_gpu,multifibres_gpu, signals_gpu, isosignals_gpu);

	thrust::device_vector<int> multirecords_gpu;
	thrust::device_vector<float> rf0_gpu, rtau_gpu, rs0_gpu, rd_gpu, rdstd_gpu, rth_gpu, rph_gpu, rf_gpu, rlikelihood_energy_gpu;

	prepare_data_gpu_MCMC_record(nvox, multirecords_gpu, rf0_gpu, rtau_gpu, rs0_gpu, rd_gpu, rdstd_gpu, rth_gpu, rph_gpu, rf_gpu, rlikelihood_energy_gpu);

	runmcmc_record(datam_gpu, bvals_gpu, alpha_gpu,beta_gpu, fibres_gpu, multifibres_gpu, signals_gpu, isosignals_gpu, rand(), gpu_log, multirecords_gpu, rf0_gpu, rtau_gpu, rs0_gpu, rd_gpu, rdstd_gpu, rth_gpu, rph_gpu, rf_gpu, rlikelihood_energy_gpu);

	/////// FINISH ALL VOXELS  ///////
	record_finish_voxels(vol2matrixkey, matrix2volkey, mask, multirecords_gpu, rf0_gpu, rtau_gpu, rs0_gpu, rd_gpu, rdstd_gpu, rth_gpu, rph_gpu, rf_gpu, rlikelihood_energy_gpu, nvox);
}


// Correct bvals/bvecs accounting for Gradient Nonlinearities
// ColumnVector grad_nonlin has 9 entries, corresponding to the 3 components of each of the x,y and z gradient deviation
void correct_bvals_bvecs(const Matrix& bvals,const Matrix& bvecs, const ColumnVector& grad_nonlin, const Matrix& Qform, const Matrix& Qform_inv, Matrix& bvals_c, Matrix& bvecs_c){
  	bvals_c=bvals; bvecs_c=bvecs;
  	Matrix L(3,3);  //gradient coil tensor
  	float mag;
  	L(1,1)=grad_nonlin(1);  L(1,2)=grad_nonlin(4);  L(1,3)=grad_nonlin(7);
  	L(2,1)=grad_nonlin(2);  L(2,2)=grad_nonlin(5);  L(2,3)=grad_nonlin(8);
  	L(3,1)=grad_nonlin(3);  L(3,2)=grad_nonlin(6);  L(3,3)=grad_nonlin(9);

  	IdentityMatrix Id(3);
  
  	for (int l=1; l<=bvals.Ncols(); l++){
    		if (bvals(1,l)>0){ //do not correct b0s
      			//Rotate bvecs to scanner's coordinate system
      			ColumnVector bvec_tmp(3);
      			bvec_tmp=Qform*bvecs.Column(l);
      			bvec_tmp(1)=-bvec_tmp(1); //Sign-flip X coordinate

      			//Correct for grad-nonlin in scanner's coordinate system
     		 	bvecs_c.Column(l)=(Id+L)*bvec_tmp;//bvecs.Column(l);
      			mag=sqrt(bvecs_c(1,l)*bvecs_c(1,l)+bvecs_c(2,l)*bvecs_c(2,l)+bvecs_c(3,l)*bvecs_c(3,l));
      			if (mag!=0)
				bvecs_c.Column(l)=bvecs_c.Column(l)/mag;
      			bvals_c(1,l)=mag*mag*bvals(1,l);
      			bvec_tmp=bvecs_c.Column(l);

      			//Rotate corrected bvecs back to voxel coordinate system
      			bvec_tmp(1)=-bvec_tmp(1); //Sign-flip X coordinate
      			bvecs_c.Column(l)=Qform_inv*bvec_tmp;
    		}
  	}
}

//////   FIT  //////
void fit(	//INPUT
		const vector<ColumnVector> 	datam_vec, 
		const vector<Matrix> 		bvecs_vec,
		const vector<Matrix> 		bvals_vec,
		thrust::host_vector<double> 	datam_host,
		thrust::host_vector<double>	bvecs_host, 
		thrust::host_vector<double>	bvals_host,
		thrust::device_vector<double> 	datam_gpu, 
		thrust::device_vector<double>	bvecs_gpu, 
		thrust::device_vector<double>	bvals_gpu,
		string 				output_file,
		//OUTPUT
		thrust::device_vector<double>&	params_gpu,
		thrust::host_vector<int>&	vox_repeat,	//for get residuals with or withot f0
		int&				nrepeat)
{
	std::ofstream myfile;
	myfile.open (output_file.data(), ios::out | ios::app );
	myfile << "----------------------------------------------------- " << "\n"; 
   	myfile << "------------------- FIT IN GPU ---------------------- " << "\n"; 
   	myfile << "----------------------------------------------------- " << "\n"; 

	struct timeval t1,t2;
   	double time;
   	gettimeofday(&t1,NULL);

	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = datam_vec.size();
	int nfib= opts.nfibres.value();
	int nparams_fit = 2+3*opts.nfibres.value();
	if(opts.modelnum.value()==2) nparams_fit++;
	if(opts.f0.value()) nparams_fit++;

	if(opts.modelnum.value()==1){
		if(opts.nonlin.value()){ 
			fit_PVM_single(datam_vec,bvecs_vec,bvals_vec,datam_gpu,bvecs_gpu,bvals_gpu,opts.f0.value(),params_gpu);

			if (opts.f0.value()){
				float md,mf,f0;	

				for(int vox=0;vox<nvox;vox++){			
					md = params_gpu[vox*nparams_fit+(1)];
					mf = params_gpu[vox*nparams_fit+(2)];
					f0 = params_gpu[vox*nparams_fit+(nparams_fit-1)];
					if ((opts.nfibres.value()>0 && mf<0.05) || md>0.007 || f0>0.4){		//if true we need to repeat this voxel
						vox_repeat[nrepeat]=vox;
						nrepeat++;
					}
				}
				if(nrepeat>0){
					//prepare structures for the voxels that need to be reprocessed
					vector<ColumnVector> 	datam_repeat_vec; 
					vector<Matrix> 		bvecs_repeat_vec;
					vector<Matrix> 		bvals_repeat_vec;
					thrust::host_vector<double> 	datam_repeat_host;
					thrust::host_vector<double> 	bvecs_repeat_host;	
					thrust::host_vector<double> 	bvals_repeat_host;	
					thrust::host_vector<double> 	params_repeat_host;	
								
					prepare_data_gpu_FIT_repeat(datam_host, bvecs_host, bvals_host, vox_repeat, nrepeat, datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_host, bvecs_repeat_host, bvals_repeat_host, params_repeat_host);

					thrust::device_vector<double> datam_repeat_gpu=datam_repeat_host;
					thrust::device_vector<double> bvecs_repeat_gpu=bvecs_repeat_host;
					thrust::device_vector<double> bvals_repeat_gpu=bvals_repeat_host;	
					thrust::device_vector<double> params_repeat_gpu=params_repeat_host;
				
		 			fit_PVM_single(datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_gpu, bvecs_repeat_gpu, bvals_repeat_gpu, false, params_repeat_gpu);

					//mix all the parameteres: repeated and not repeated
					mix_params(params_repeat_gpu,vox_repeat, nrepeat, params_gpu);
				}
	  		}
		}else{
			fit_PVM_single_c(datam_vec,bvecs_vec,bvals_vec,datam_gpu,bvecs_gpu,bvals_gpu,opts.f0.value(),params_gpu);

			if (opts.f0.value()){
				float md,mf,f0;	
	
				for(int vox=0;vox<nvox;vox++){		
					md = params_gpu[vox*nparams_fit+(1)];
					mf = params_gpu[vox*nparams_fit+(2)];
					f0 = params_gpu[vox*nparams_fit+(nparams_fit-1)];
					if ((opts.nfibres.value()>0 && mf<0.05) || md>0.007 || f0>0.4){		//if true we need to repeat this voxel
						vox_repeat[nrepeat]=vox;
						nrepeat++;
					}
				}
				if(nrepeat>0){
					//prepare structures for the voxels that need to be reprocessed
					vector<ColumnVector> 	datam_repeat_vec; 
					vector<Matrix> 		bvecs_repeat_vec;
					vector<Matrix> 		bvals_repeat_vec;
					thrust::host_vector<double> 	datam_repeat_host;
					thrust::host_vector<double> 	bvecs_repeat_host;	
					thrust::host_vector<double> 	bvals_repeat_host;	
					thrust::host_vector<double> 	params_repeat_host;	
								
					prepare_data_gpu_FIT_repeat(datam_host, bvecs_host, bvals_host, vox_repeat, nrepeat, datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_host, bvecs_repeat_host, bvals_repeat_host, params_repeat_host);

					thrust::device_vector<double> datam_repeat_gpu=datam_repeat_host;
					thrust::device_vector<double> bvecs_repeat_gpu=bvecs_repeat_host;
					thrust::device_vector<double> bvals_repeat_gpu=bvals_repeat_host;	
					thrust::device_vector<double> params_repeat_gpu=params_repeat_host;
				
		 			fit_PVM_single_c(datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_gpu, bvecs_repeat_gpu, bvals_repeat_gpu, false, params_repeat_gpu);

					//mix all the parameteres: repeated and not repeated
					mix_params(params_repeat_gpu ,vox_repeat, nrepeat, params_gpu);
				}
	  		}
		}
	}else{
      		//model 2 : non-mono-exponential
		fit_PVM_single_c(datam_vec,bvecs_vec,bvals_vec,datam_gpu,bvecs_gpu,bvals_gpu,opts.f0.value(),params_gpu);
	
		fit_PVM_multi(datam_gpu,bvecs_gpu,bvals_gpu,nvox,opts.f0.value(),params_gpu);	

		if (opts.f0.value()){
				float md,mf,f0;	
	
				for(int vox=0;vox<nvox;vox++){			
					md = params_gpu[vox*nparams_fit+(1)];
					mf = params_gpu[vox*nparams_fit+(3)];
					f0 = params_gpu[vox*nparams_fit+(nparams_fit-1)];
					if ((opts.nfibres.value()>0 && mf<0.05) || md>0.007 || f0>0.4){		//if true we need to repeat this voxel
						vox_repeat[nrepeat]=vox;
						nrepeat++;
					}
				}
				if(nrepeat>0){
					//prepare structures for the voxels that need to be reprocessed
					vector<ColumnVector> 	datam_repeat_vec; 
					vector<Matrix> 		bvecs_repeat_vec;
					vector<Matrix> 		bvals_repeat_vec;
					thrust::host_vector<double> 	datam_repeat_host;
					thrust::host_vector<double> 	bvecs_repeat_host;	
					thrust::host_vector<double> 	bvals_repeat_host;	
					thrust::host_vector<double> 	params_repeat_host;		
								
					prepare_data_gpu_FIT_repeat(datam_host, bvecs_host, bvals_host, vox_repeat, nrepeat, datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_host, bvecs_repeat_host,  bvals_repeat_host, params_repeat_host);

					thrust::device_vector<double> datam_repeat_gpu=datam_repeat_host;
					thrust::device_vector<double> bvecs_repeat_gpu=bvecs_repeat_host;
					thrust::device_vector<double> bvals_repeat_gpu=bvals_repeat_host;	
					thrust::device_vector<double> params_repeat_gpu=params_repeat_host;
				
		 			fit_PVM_single_c(datam_repeat_vec, bvecs_repeat_vec, bvals_repeat_vec, datam_repeat_gpu, bvecs_repeat_gpu, bvals_repeat_gpu, false, params_repeat_gpu);

					fit_PVM_multi(datam_repeat_gpu, bvecs_repeat_gpu, bvals_repeat_gpu, nrepeat, false, params_repeat_gpu);
		
					//mix all the parameteres: repeated and not repeated
					mix_params(params_repeat_gpu ,vox_repeat, nrepeat,params_gpu);
				}
	  		}	
	}

	gettimeofday(&t2,NULL);
    	time=timeval_diff(&t2,&t1);
   	myfile << "TIME TOTAL: " << time << " seconds\n"; 
	myfile << "--------------------------------------------" << "\n\n" ; 
	myfile.close();

}

//prepare the structures for copy all neccesary data to FIT in GPU
void prepare_data_gpu_FIT(	//INPUT
				const Matrix				datam,
				const Matrix				bvecs,
				const Matrix				bvals,
				const Matrix	 			gradm, 
				const Matrix 				Qform, 
				const Matrix 				Qform_inv,
				//OUTPUT
				vector<ColumnVector>&			datam_vec,
				vector<Matrix>&				bvecs_vec,
				vector<Matrix>&				bvals_vec,
				thrust::host_vector<double>&   		datam_host,	//data prepared for copy to GPU
				thrust::host_vector<double>&		bvecs_host,				
				thrust::host_vector<double>&		bvals_host,
				thrust::host_vector<double>&		alpha_host,
				thrust::host_vector<double>&		beta_host,
				thrust::host_vector<double>&		params_host,
				thrust::host_vector<float>&		tau_host)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	int nvox = datam.Ncols(); 

	datam_vec.resize(nvox);
	datam_host.resize(nvox*NDIRECTIONS); 
	for(int vox=0;vox<nvox;vox++){
		datam_vec[vox]=datam.Column(vox+1);
		for(int j=0;j<NDIRECTIONS;j++){
			datam_host[vox*NDIRECTIONS+j]=datam(j+1,vox+1);
		}
	}

	bvecs_vec.resize(nvox);
	bvals_vec.resize(nvox);
	bvecs_host.resize(nvox*bvecs.Nrows()*bvecs.Ncols());
	bvals_host.resize(nvox*bvals.Ncols());

	alpha_host.resize(nvox*bvecs.Ncols());
	beta_host.resize(nvox*bvecs.Ncols());
	
	ColumnVector alpha,beta;

	if (opts.grad_file.set()){
		for(int vox=0;vox<nvox;vox++){
			correct_bvals_bvecs(bvals,bvecs, gradm.Column(vox+1),Qform,Qform_inv,bvals_vec[vox],bvecs_vec[vox]); //correct for gradient nonlinearities
 			MISCMATHS::cart2sph(bvecs_vec[vox],alpha,beta);
			
			for(int dir=0;dir<NDIRECTIONS;dir++){
				bvecs_host[vox*NDIRECTIONS*3+dir] = bvecs_vec[vox](1,dir+1);
				bvecs_host[vox*NDIRECTIONS*3+NDIRECTIONS+dir] = bvecs_vec[vox](2,dir+1);
				bvecs_host[vox*NDIRECTIONS*3+NDIRECTIONS*2+dir] = bvecs_vec[vox](3,dir+1);
				bvals_host[vox*NDIRECTIONS+dir] = bvals_vec[vox](1,dir+1);

				alpha_host[vox*NDIRECTIONS+dir] = alpha(dir+1);
        			beta_host[vox*NDIRECTIONS+dir] = beta(dir+1);
			}
		}
		
	}else{
 		MISCMATHS::cart2sph(bvecs,alpha,beta);

		for(int vox=0;vox<nvox;vox++){
			bvecs_vec[vox]=bvecs;
			bvals_vec[vox]=bvals;
			for(int dir=0;dir<NDIRECTIONS;dir++){
				bvecs_host[vox*NDIRECTIONS*3+dir] = bvecs(1,dir+1);
				bvecs_host[vox*NDIRECTIONS*3+NDIRECTIONS+dir] = bvecs(2,dir+1);
				bvecs_host[vox*NDIRECTIONS*3+NDIRECTIONS*2+dir] = bvecs(3,dir+1);
        			bvals_host[vox*NDIRECTIONS+dir] = bvals(1,dir+1);
			
				alpha_host[vox*NDIRECTIONS+dir] = alpha(dir+1);
        			beta_host[vox*NDIRECTIONS+dir] = beta(dir+1);
			}
		}
	}
	
	int nfib= opts.nfibres.value();
	int nparams;

	if(opts.f0.value()) nparams=3+nfib*3;
	else nparams=2+nfib*3;	
	if(opts.modelnum.value()==2) nparams++;

	params_host.resize(nvox*nparams);
	tau_host.resize(nvox);
}

//prepare the structures for copy all neccesary data to FIT in GPU when is repeated because f0. Only some voxels
void prepare_data_gpu_FIT_repeat(	//INPUT
					thrust::host_vector<double>   		datam_host,	
					thrust::host_vector<double>		bvecs_host,				
					thrust::host_vector<double>		bvals_host,
					thrust::host_vector<int>		vox_repeat,
					int					nrepeat,
					//OUTPUT
					vector<ColumnVector>&			datam_repeat_vec,
					vector<Matrix>&				bvecs_repeat_vec,
					vector<Matrix>&				bvals_repeat_vec,
					thrust::host_vector<double>&   		datam_repeat_host,	//data prepared for copy to GPU
					thrust::host_vector<double>&		bvecs_repeat_host,				
					thrust::host_vector<double>&		bvals_repeat_host,
					thrust::host_vector<double>&		params_repeat_host)
{
	xfibresOptions& opts = xfibresOptions::getInstance();

	ColumnVector datam(NDIRECTIONS);
	Matrix	bvecs(3,NDIRECTIONS);
	Matrix	bvals(1,NDIRECTIONS);

	datam_repeat_vec.resize(nrepeat);
	datam_repeat_host.resize(nrepeat*NDIRECTIONS); 
	bvecs_repeat_vec.resize(nrepeat);
	bvals_repeat_vec.resize(nrepeat);
	bvecs_repeat_host.resize(nrepeat*3*NDIRECTIONS);
	bvals_repeat_host.resize(nrepeat*NDIRECTIONS);

	for(int vox=0;vox<nrepeat;vox++){
		for(int dir=0;dir<NDIRECTIONS;dir++){
			datam(dir+1)= datam_host[vox_repeat[vox]*NDIRECTIONS+dir]; 
			datam_repeat_host[vox*NDIRECTIONS+dir]=datam_host[vox_repeat[vox]*NDIRECTIONS+dir];

			bvecs_repeat_host[vox*NDIRECTIONS*3+dir] = bvecs_host[vox_repeat[vox]*NDIRECTIONS*3+dir];
			bvecs_repeat_host[vox*NDIRECTIONS*3+NDIRECTIONS+dir] = bvecs_host[vox_repeat[vox]*NDIRECTIONS*3+NDIRECTIONS+dir];
			bvecs_repeat_host[vox*NDIRECTIONS*3+NDIRECTIONS*2+dir] = bvecs_host[vox_repeat[vox]*NDIRECTIONS*3+NDIRECTIONS*2+dir];
        		bvals_repeat_host[vox*NDIRECTIONS+dir] = bvals_host[vox_repeat[vox]*NDIRECTIONS+dir];
			
			bvecs(1,dir+1)= bvecs_host[vox_repeat[vox]*NDIRECTIONS*3+dir];
			bvecs(2,dir+1)= bvecs_host[vox_repeat[vox]*NDIRECTIONS*3+NDIRECTIONS+dir];
			bvecs(3,dir+1)= bvecs_host[vox_repeat[vox]*NDIRECTIONS*3+NDIRECTIONS*2+dir];
			bvals(1,dir+1)= bvals_host[vox_repeat[vox]*NDIRECTIONS+dir];
		}
		datam_repeat_vec[vox]=datam;	
		bvecs_repeat_vec[vox]=bvecs;
		bvals_repeat_vec[vox]=bvals;
	}
	
	int nfib= opts.nfibres.value();
	int nparams;

	nparams=2+nfib*3;	
	if(opts.modelnum.value()==2) nparams++;

	params_repeat_host.resize(nrepeat*nparams);
}


void mix_params(	//INPUT
			thrust::device_vector<double>   		params_repeat_gpu,
			thrust::host_vector<int>			vox_repeat,
			int						nrepeat,
			//INPUT-OUTPUT
			thrust::device_vector<double>&   		params_gpu)
{
	xfibresOptions& opts = xfibresOptions::getInstance();
	int nfib= opts.nfibres.value();
	int nparams;
	nparams=2+nfib*3;	
	if(opts.modelnum.value()==2) nparams++;

	for(int vox=0;vox<nrepeat;vox++){
		for(int par=0;par<nparams;par++){
			params_gpu[vox_repeat[vox]*(nparams+1)+par] = params_repeat_gpu[vox*nparams+par]; //(nparams+1) to count f0
		}
		params_gpu[vox_repeat[vox]*(nparams+1)+nparams] = 0.001;	//pvmf0=0.001
	}
}

void prepare_data_gpu_MCMC(	//INPUT
				int 					nvox,
				int 					nfib,
				//OUTPUT
				thrust::host_vector<double>&		signals_host,
				thrust::host_vector<double>&		isosignals_host,
				thrust::host_vector<FibreGPU>& 		fibres_host,
				thrust::host_vector<MultifibreGPU>& 	multifibres_host)
{ 	
	signals_host.resize(nvox*nfib*NDIRECTIONS);
	isosignals_host.resize(nvox*NDIRECTIONS);	
	fibres_host.resize(nvox*nfib);	
	multifibres_host.resize(nvox);
}

void prepare_data_gpu_MCMC_record(	//INPUT
					int 						nvox,
					//OUTPUT
					thrust::device_vector<int>&			multirecords_gpu,
					thrust::device_vector<float>&			rf0_gpu,
					thrust::device_vector<float>&			rtau_gpu,
					thrust::device_vector<float>&			rs0_gpu,
					thrust::device_vector<float>&			rd_gpu,
					thrust::device_vector<float>&			rdstd_gpu,
					thrust::device_vector<float>&			rth_gpu,
					thrust::device_vector<float>&			rph_gpu,
					thrust::device_vector<float>&			rf_gpu,
					thrust::device_vector<float>&			rlikelihood_energy_gpu)
{ 	
	xfibresOptions& opts = xfibresOptions::getInstance();

	int nfib = opts.nfibres.value();	
	int nrecords = (opts.njumps.value()/opts.sampleevery.value());   
	
	multirecords_gpu.resize(nvox*nrecords); 
	if(opts.f0.value()) rf0_gpu.resize(nvox*nrecords); 
	if(opts.rician.value()) rtau_gpu.resize(nvox*nrecords);  
	rs0_gpu.resize(nvox*nrecords);  
	rd_gpu.resize(nvox*nrecords);
	if(opts.modelnum.value()==2) rdstd_gpu.resize(nvox*nrecords);  
	rth_gpu.resize(nvox*nrecords*nfib);  
	rph_gpu.resize(nvox*nrecords*nfib);  
	rf_gpu.resize(nvox*nrecords*nfib);  
	rlikelihood_energy_gpu.resize(nvox*nrecords); 
}

void record_finish_voxels(	//INPUT
				const NEWIMAGE::volume<int> 			vol2matrixkey,
				const NEWMAT::Matrix				matrix2volkey,
				const NEWIMAGE::volume<float>			mask,
				thrust::device_vector<int>&			multirecords_gpu,
				thrust::device_vector<float>&			rf0_gpu,
				thrust::device_vector<float>&			rtau_gpu,
				thrust::device_vector<float>&			rs0_gpu,
				thrust::device_vector<float>&			rd_gpu,
				thrust::device_vector<float>&			rdstd_gpu,
				thrust::device_vector<float>&			rth_gpu,
				thrust::device_vector<float>&			rph_gpu,
				thrust::device_vector<float>&			rf_gpu,
				thrust::device_vector<float>&			rlikelihood_energy_gpu,
				int 						nvox)
{
	xfibresOptions& opts = xfibresOptions::getInstance();

	int nfib = opts.nfibres.value();	
	int nrecords = (opts.njumps.value()/opts.sampleevery.value());   

	thrust::host_vector<int> multirecords_host;
	thrust::host_vector<float> rf0_host,rtau_host,rs0_host,rd_host,rdstd_host,rth_host,rph_host,rf_host,rlikelihood_energy_host;

	multirecords_host.resize(nvox*nrecords);
	rf0_host.resize(nvox*nrecords);
	rtau_host.resize(nvox*nrecords);
	rs0_host.resize(nvox*nrecords);
	rd_host.resize(nvox*nrecords);
	rdstd_host.resize(nvox*nrecords);
	rth_host.resize(nvox*nfib*nrecords);
	rph_host.resize(nvox*nfib*nrecords);
	rf_host.resize(nvox*nfib*nrecords);
	rlikelihood_energy_host.resize(nvox*nrecords);

	thrust::copy(multirecords_gpu.begin(), multirecords_gpu.end(), multirecords_host.begin());
	if(opts.f0.value()) thrust::copy(rf0_gpu.begin(), rf0_gpu.end(), rf0_host.begin());
	if(opts.rician.value()) thrust::copy(rtau_gpu.begin(), rtau_gpu.end(), rtau_host.begin());
	thrust::copy(rs0_gpu.begin(), rs0_gpu.end(), rs0_host.begin());
	thrust::copy(rd_gpu.begin(), rd_gpu.end(), rd_host.begin());
	if(opts.modelnum.value()==2) thrust::copy(rdstd_gpu.begin(), rdstd_gpu.end(), rdstd_host.begin());
	thrust::copy(rth_gpu.begin(), rth_gpu.end(), rth_host.begin());
	thrust::copy(rph_gpu.begin(), rph_gpu.end(), rph_host.begin());
	thrust::copy(rf_gpu.begin(), rf_gpu.end(), rf_host.begin());	
	thrust::copy(rlikelihood_energy_gpu.begin(), rlikelihood_energy_gpu.end(), rlikelihood_energy_host.begin());	

	Samples samples(vol2matrixkey,matrix2volkey,nvox,NDIRECTIONS);

	float ard,arf0,artau,ardstd,ars0,arlikelihood_energy;	
	float *arth = new float[nfib];
    	float *arph = new float[nfib]; 
    	float *arf = new float[nfib];
	int samp;

	for(int vox=0;vox<nvox;vox++){
		for(int rec=0;rec<nrecords;rec++){	
			ard=rd_host[(vox*nrecords)+rec];
			if(opts.f0.value()){	
				arf0=rf0_host[(vox*nrecords)+rec];
			}

			if(opts.rician.value()){	
				artau=rtau_host[(vox*nrecords)+rec];
			}

			if(opts.modelnum.value()==2){	
				ardstd=rdstd_host[(vox*nrecords)+rec];
			}
		
			ars0=rs0_host[(vox*nrecords)+rec];
		
			arlikelihood_energy=rlikelihood_energy_host[(vox*nrecords)+rec];	

			for(int j=0;j<nfib;j++){
				arth[j]=rth_host[(vox*nfib*nrecords)+(j*nrecords)+rec];
				arph[j]=rph_host[(vox*nfib*nrecords)+(j*nrecords)+rec];
				arf[j]=rf_host[(vox*nfib*nrecords)+(j*nrecords)+rec];

			}

			samp=multirecords_host[(vox*nrecords)+rec];	
			samples.record(ard,arf0,artau,ardstd,ars0,arlikelihood_energy,arth,arph,arf,vox+1,samp);
		}	
		samples.finish_voxel(vox+1);
   	}

	samples.save(mask);
}

