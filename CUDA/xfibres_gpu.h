#include "newmat.h"
#include "newimage/newimageall.h"
#include <string>

int xfibres_gpu(char *subjdir,int slice, int nextras, char** extras);

void xfibres_gpu(	//INPUT
			const NEWMAT::Matrix			datam,
			const NEWMAT::Matrix			bvals,
			const NEWMAT::Matrix			bvecs,
			const NEWMAT::Matrix	 		gradm, 
			const NEWMAT::Matrix 			Qform, 
			const NEWMAT::Matrix 			Qform_inv,
			const NEWIMAGE::volume<int> 		vol2matrixkey,
			const NEWMAT::Matrix			matrix2volkey,
			const NEWIMAGE::volume<float>		mask,
			const int 				slice,
			const char*				subjdir);

