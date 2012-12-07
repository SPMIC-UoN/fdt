#include "newmat.h"
#include "newimage/newimageall.h"
#include "xfibresoptions.h"


using namespace Xfibres;

////////////////////////////////////////////
//       MCMC SAMPLE STORAGE
////////////////////////////////////////////

class Samples{
  	xfibresOptions& opts;
  	Matrix m_dsamples;
  	Matrix m_d_stdsamples;
  	Matrix m_S0samples;
  	Matrix m_f0samples;
  	Matrix m_lik_energy;

	//   // storing signal
	//   Matrix m_mean_sig;
	//   Matrix m_std_sig;
	//   Matrix m_sig2;

  	vector<Matrix> m_thsamples;
  	vector<Matrix> m_phsamples;
  	vector<Matrix> m_fsamples;
  	vector<Matrix> m_lamsamples;

  	//for storing means
  	RowVector m_mean_dsamples;
  	RowVector m_mean_d_stdsamples;
  	RowVector m_mean_S0samples;
  	RowVector m_mean_f0samples;
  	RowVector m_mean_tausamples;
  	vector<Matrix> m_dyadic_vectors;
  	vector<RowVector> m_mean_fsamples;
  	vector<RowVector> m_mean_lamsamples;

  	//float m_sum_d;  changed GPU version
  	//float m_sum_d_std;  changed GPU version
  	//float m_sum_S0;  changed GPU version
  	//float m_sum_f0;  changed GPU version
  	//float m_sum_tau;  changed GPU version
  	//vector<SymmetricMatrix> m_dyad;  changed GPU version
  	//vector<float> m_sum_f;  changed GPU version
  	//vector<float> m_sum_lam;  changed GPU version
  	//ColumnVector m_vec;  changed GPU version

  	/////////////// GPU version /////////////////////
  	float *m_sum_d;
  	float *m_sum_S0;
  	float *m_sum_d_std;
  	float *m_sum_f0;
  	float *m_sum_tau;

  	vector<SymmetricMatrix> *m_dyad;
  	vector<float>  *m_sum_f;
  	vector<float> *m_sum_lam;
  	ColumnVector *m_vec;
  	////////////////////////////////////////////////
  
  	int m_nsamps;
  
  	NEWIMAGE::volume<int> m_vol2matrixkey;
  	Matrix m_matrix2volkey;
  	NEWIMAGE::volume<int> m_beenhere;

	public:

  	Samples(NEWIMAGE::volume<int> vol2matrixkey,Matrix matrix2volkey,int nvoxels,int nmeasures);
    
	void record(float rd,float rf0,float rtau,float rdstd,float rs0,float rlikelihood_energy, float *rth,float *rph, float *rf, int vox, int samp);
  
  	void finish_voxel(int vox);
  
  	void save(const NEWIMAGE::volume<float>& mask);
};
