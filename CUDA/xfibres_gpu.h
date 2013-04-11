/*  xfibres_gpu.h

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

/*  CCOPYRIGHT  */

#include "newmat.h"
#include "newimage/newimageall.h"
#include <string>

//int xfibres_gpu(char *subjdir,int slice, int nextras, char** extras);

void join_subParts(string name, int size_part, int nsubparts, int size_sub_part, int last_sub_part, bool mean);

void xfibres_gpu(	//INPUT
			const NEWMAT::Matrix			datam,
			const NEWMAT::Matrix			bvals,
			const NEWMAT::Matrix			bvecs,
			const NEWMAT::Matrix	 		gradm, 
			int					idpart);

