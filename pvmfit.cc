/*  Copyright (C) 2009 University of Oxford  */



/*  CCOPYRIGHT  */

#include <iostream>
#include <cmath>
#include "miscmaths/miscmaths.h"
#include "miscmaths/nonlin.h"
#include "newmat.h"
#include "pvmfitOptions.h"
#include "newimage/newimageall.h"

const float two_pi=0.636619772;

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace PVMFIT;
using namespace NEWIMAGE;

inline SymmetricMatrix vec2tens(ColumnVector& Vec){
  SymmetricMatrix tens(3);
  tens(1,1)=Vec(1);
  tens(2,1)=Vec(2);
  tens(3,1)=Vec(3);
  tens(2,2)=Vec(4);
  tens(3,2)=Vec(5);
  tens(3,3)=Vec(6);
  return tens;
}
Matrix form_Amat(const Matrix& r,const Matrix& b)
{
  Matrix A(r.Ncols(),7);
  Matrix tmpvec(3,1), tmpmat;
  
  for( int i = 1; i <= r.Ncols(); i++){
    tmpvec << r(1,i) << r(2,i) << r(3,i);
    tmpmat = tmpvec*tmpvec.t()*b(1,i);
    A(i,1) = tmpmat(1,1);
    A(i,2) = 2*tmpmat(1,2);
    A(i,3) = 2*tmpmat(1,3);
    A(i,4) = tmpmat(2,2);
    A(i,5) = 2*tmpmat(2,3);
    A(i,6) = tmpmat(3,3);
    A(i,7) = 1;
  }
  return A;
}
// dtifit
class DTI {
public: 

  DTI(const ColumnVector& iY,
      const Matrix& iAmat):Y(iY),pinvAmat(iAmat){
  
    npts = Y.Nrows();
    v1.ReSize(3);
    v2.ReSize(3);
    v3.ReSize(3);

  }
  ~DTI(){}
  void fit();
  float get_fa()const{return fa;}
  float get_md()const{return md;}
  float get_s0()const{return s0;}
  float get_mo()const{return mo;}
  ColumnVector get_v1()const{return v1;}
  ColumnVector get_v2()const{return v2;}
  ColumnVector get_v3()const{return v3;}
  ColumnVector get_v(const int& i)const{if(i==1)return v1;else if(i==2)return v2;else return v3;}
  void print();

private:
  const ColumnVector& Y;
  const Matrix& pinvAmat;
  int npts;
  ColumnVector v1,v2,v3;
  float l1,l2,l3;
  float fa,s0,md,mo;

};
void DTI::fit(){
  ColumnVector logS(npts);
  ColumnVector Dvec(7);
  SymmetricMatrix tens;
  Matrix Vd;
  DiagonalMatrix Dd(3);

  for (int i=1;i<=npts; i++){
    if(Y(i)>0)
      logS(i)=log(Y(i));
    else
      logS(i)=0;
  }

  Dvec = -pinvAmat*logS;

  if(Dvec(7)>-23)
    s0=exp(-Dvec(7));
  else
    s0=Y.MaximumAbsoluteValue();
  
  for (int i=1;i<=Y.Nrows();i++){
    if(s0<Y.Sum()/Y.Nrows()){ s0=Y.MaximumAbsoluteValue();  }
    logS(i)=(Y(i)/s0)>0.01 ? log(Y(i)):log(0.01*s0);
  }
  Dvec = -pinvAmat*logS;
  s0=exp(-Dvec(7));
  
  if(s0<Y.Sum()/Y.Nrows()){ s0=Y.Sum()/Y.Nrows();  }
  tens = vec2tens(Dvec);
  
  EigenValues(tens,Dd,Vd);
  
  md = Dd.Sum()/Dd.Nrows();

  l1 = Dd(3,3);
  l2 = Dd(2,2);
  l3 = Dd(1,1);
  v1 = Vd.Column(3);
  v2 = Vd.Column(2);
  v3 = Vd.Column(1);

  float e1=l1-md, e2=l2-md, e3=l3-md;
  float n = (e1 + e2 - 2*e3)*(2*e1 - e2 - e3)*(e1 - 2*e2 + e3);
  float d = (e1*e1 + e2*e2 + e3*e3 - e1*e2 - e2*e3 - e1*e3);
  d = sqrt(MAX(0, d));
  d = 2*d*d*d;
  mo = MIN(MAX(d ? n/d : 0.0, -1),1);

  float numer=1.5*((l1-md)*(l1-md)+(l2-md)*(l2-md)+(l3-md)*(l3-md));
  float denom=(l1*l1+l2*l2+l3*l3);
  if(denom>0) fa=numer/denom;
  else fa=0;
  if(fa>0){fa=sqrt(fa);}
  else{fa=0;}


}
void DTI::print(){
  cout << "DTI FIT RESULTS " << endl;
  cout << "S0   :" << s0 << endl;
  cout << "MD   :" << md << endl;
  cout << "FA   :" << fa << endl;; 
  cout << "MO   :" << mo << endl;; 
}

// nonlinear optimisation

// nonlinear class for despot1_hifi
class PVMNonlinCF : public NonlinCF {
public:
  PVMNonlinCF(const ColumnVector& iY,
	      const ColumnVector& ibvals,
	      const ColumnVector& isina,const ColumnVector& icosa,const ColumnVector& ibeta,
	      const int& nfibres):Y(iY),bvals(ibvals),sinalpha(isina),cosalpha(icosa),beta(ibeta){

    npts    = Y.Nrows();
    nparams = nfibres*3 + 2;
    nfib    = nfibres;

 //    OUT(Y.t());
//     OUT(npts);
//     OUT(nparams);
//     OUT(bvals.t());
//     OUT(sinalpha.t());
//     OUT(beta.t());

  }
  ~PVMNonlinCF(){}
  
  NEWMAT::ReturnMatrix grad(const NEWMAT::ColumnVector& p)const;
  boost::shared_ptr<BFMatrix> hess(const NEWMAT::ColumnVector&p,boost::shared_ptr<BFMatrix> iptr)const;
  double cf(const NEWMAT::ColumnVector& p)const;
  NEWMAT::ReturnMatrix forwardModel(const NEWMAT::ColumnVector& p)const;

private:
  const ColumnVector& Y;
  const ColumnVector& bvals;
  const ColumnVector& sinalpha;
  const ColumnVector& cosalpha;
  const ColumnVector& beta;
  
  float npts;
  int   nparams;
  int   nfib;
};
NEWMAT::ReturnMatrix PVMNonlinCF::forwardModel(const NEWMAT::ColumnVector& p)const{
  ColumnVector ret(npts);
  
  double val;
  float angtmp;
  ColumnVector f(nfib);

  float sumf=0;
  for(int i=3,j=1;i<=p.Nrows();i+=3,j++){
    f(j) = abs(two_pi*atan(p(i)));  
    sumf  += f(j); 
  }      
  for(int i=1;i<=Y.Nrows();i++){
    val = 0.0;
    for(int k=3;k<=p.Nrows();k+=3){
      angtmp = cos(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)) + cosalpha(i)*cos(p(k+1));
      
      val += f(k/3)*exp(-bvals(i)*p(2)*angtmp*angtmp);
    }
    ret(i) = (p(1)*((1-sumf)*exp(-bvals(i)*p(2))+val)); 

  }

  ret.Release();
  return(ret);
}

double PVMNonlinCF::cf(const NEWMAT::ColumnVector& p)const{
  // p(1) = S0
  // p(2) = d
  // p(3) = f1
  // p(4) = th1
  // p(5) = ph1
  // etc.
  
  //cout << "CF" << endl;
  //OUT(p.t());

  double cfv = 0.0;
  double err;
  float angtmp;
  ColumnVector f(nfib);

  float sumf=0;
  for(int i=3,j=1;i<=p.Nrows();i+=3,j++){
    f(j)  = abs(two_pi*atan(p(i)));
    sumf += f(j);
  }

  for(int i=1;i<=Y.Nrows();i++){
    err = 0.0;
    for(int k=3;k<=p.Nrows();k+=3){
      angtmp = cos(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)) + cosalpha(i)*cos(p(k+1));
      //err += p(k)*exp(-bvals(i)*p(2)*angtmp*angtmp);
      err += f(k/3)*exp(-bvals(i)*p(2)*angtmp*angtmp);
    }
    err = (p(1)*((1-sumf)*exp(-bvals(i)*p(2))+err) - Y(i)); 

    cfv += err*err; 
  }
  //OUT(cfv);
  return(cfv);
}
NEWMAT::ReturnMatrix PVMNonlinCF::grad(const NEWMAT::ColumnVector& p)const{
  //cout << "gradv" << endl;
  //OUT(p.t());

  NEWMAT::ColumnVector gradv(p.Nrows());
  gradv = 0.0;
  double dval1;
  double dval2;
  float angtmp,tmpval;
  ColumnVector f(nfib);

  float sumf=0;
  for(int i=3,j=1;i<=p.Nrows();i+=3,j++){
    f(j) = abs(two_pi*atan(p(i)));
    sumf += f(j);
  }
      
  for(int i=1;i<=Y.Nrows();i++){
    // calculate difference between signal and data
    dval1 = 0.0;
    dval2 = 0.0;
    for(int k=3;k<=p.Nrows();k+=3){
      angtmp = cos(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)) + cosalpha(i)*cos(p(k+1));
      angtmp *= angtmp;
      //      tmpval = p(k)*exp(-bvals(i)*p(2)*angtmp);
      tmpval = f(k/3)*exp(-bvals(i)*p(2)*angtmp);
      dval1 += tmpval;
      dval2 += -bvals(i)*angtmp*tmpval;
    }
    dval1 = (p(1)*((1-sumf)*exp(-bvals(i)*p(2))+dval1) - Y(i)); 
    dval2 = (p(1)*(-bvals(i)*(1-sumf)*exp(-bvals(i)*p(2))+dval2)); 
    
    gradv(1) += 2 * dval1 * (dval1+Y(i))/p(1);
    gradv(2) += 2 * dval1 * dval2;

    for(int k=3;k<=p.Nrows();k+=3){
      angtmp = cos(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)) + cosalpha(i)*cos(p(k+1));
      tmpval = f(k/3)*exp(-bvals(i)*p(2)*angtmp*angtmp);

      gradv(k)   += 2 * dval1 * p(1) * (-exp(-bvals(i)*p(2)) + tmpval) / f(k/3) * sign(p(k)) / (1+p(k)*p(k));
      gradv(k+1) += 2 * dval1 * p(1) * tmpval * (-bvals(i)*p(2)*2*angtmp*(cos(p(k+2)-beta(i))*sinalpha(i)*cos(p(k+1)) - cosalpha(i)*sin(p(k+1))));
      gradv(k+2) += 2 * dval1 * p(1) * tmpval * (bvals(i)*p(2)*2*angtmp*(sin(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1))));
      
    }
  }
  
  //OUT(gradv.t());

  gradv.Release();
  return(gradv);
}
// this uses Gauss-Newton approximation
boost::shared_ptr<BFMatrix> PVMNonlinCF::hess(const NEWMAT::ColumnVector& p,boost::shared_ptr<BFMatrix> iptr)const{
  //cout << "hessian" << endl;
  //OUT(p.t());

  boost::shared_ptr<BFMatrix>   hessm;

  if (iptr && iptr->Nrows()==(unsigned int)p.Nrows() && iptr->Ncols()==(unsigned int)p.Nrows()) hessm = iptr;
  else hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(),p.Nrows()));

  double dval1;
  double dval2;
  float angtmp,tmpval;
  ColumnVector f(nfib);
  float sumf=0;
  for(int i=3,j=1;i<=p.Nrows();i+=3,j++){
    f(j) = abs(two_pi*atan(p(i)));
    sumf += f(j);
  }

  Matrix J(Y.Nrows(),nparams);
  for(int i=1;i<=Y.Nrows();i++){

    dval1 = 0.0;
    dval2 = 0.0;
    for(int k=3;k<=p.Nrows();k+=3){
      angtmp = cos(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)) + cosalpha(i)*cos(p(k+1));
      angtmp *= angtmp;
      //      tmpval = p(k)*exp(-bvals(i)*p(2)*angtmp);
      tmpval = f(k/3)*exp(-bvals(i)*p(2)*angtmp);
      dval1 += tmpval;
      dval2 += -bvals(i)*angtmp*tmpval;
    }
    dval1 = (p(1)*((1-sumf)*exp(-bvals(i)*p(2))+dval1) - Y(i)); 
    dval2 = (p(1)*(-bvals(i)*(1-sumf)*exp(-bvals(i)*p(2))+dval2)); 

    J(i,1) = (dval1+Y(i))/p(1);
    J(i,2) = dval2;

    for(int k=3;k<=p.Nrows();k+=3){
      angtmp = cos(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)) + cosalpha(i)*cos(p(k+1));
      tmpval = f(k/3)*exp(-bvals(i)*p(2)*angtmp*angtmp);

      J(i,k)   = p(1) * (-exp(-bvals(i)*p(2)) + tmpval)/ f(k/3) * sign(p(k)) / (1+p(k)*p(k));;
      J(i,k+1) = p(1) * tmpval * (-bvals(i)*p(2)*2*angtmp*(cos(p(k+2)-beta(i))*sinalpha(i)*cos(p(k+1)) - cosalpha(i)*sin(p(k+1))));
      J(i,k+2) = p(1) * tmpval * (bvals(i)*p(2)*2*angtmp*(sin(p(k+2)-beta(i))*sinalpha(i)*sin(p(k+1)))); 
    }
  }
  
  for (int i=1; i<=p.Nrows(); i++){
    for (int j=i; j<=p.Nrows(); j++){
      dval1 = 0.0;
      for(int k=1;k<=J.Nrows();k++)
	dval1 += J(k,i)*J(k,j);
      hessm->Set(i,j,dval1);
    }
  }

  for (int j=1; j<=p.Nrows(); j++) {
    for (int i=j+1; i<=p.Nrows(); i++) {
      hessm->Set(i,j,hessm->Peek(j,i));
    }
  }

  //hessm->Print();

  return(hessm);
}



int main(int argc, char** argv)
{
  //parse command line
  pvmfitOptions& opts = pvmfitOptions::getInstance();
  int success=opts.parse_command_line(argc,argv);
  if(!success) return 1;
   if(opts.verbose.value()){
    cout<<"data file "<<opts.datafile.value()<<endl;
    cout<<"mask file "<<opts.maskfile.value()<<endl;
    cout<<"bvecs     "<<opts.bvecsfile.value()<<endl;
    cout<<"bvals     "<<opts.bvalsfile.value()<<endl;
  }
  
  // Set random seed:
  Matrix r = read_ascii_matrix(opts.bvecsfile.value());
  if(r.Nrows()>3) r=r.t();
  for(int i=1;i<=r.Ncols();i++){
    float tmpsum=sqrt(r(1,i)*r(1,i)+r(2,i)*r(2,i)+r(3,i)*r(3,i));
    if(tmpsum!=0){
      r(1,i)=r(1,i)/tmpsum;
      r(2,i)=r(2,i)/tmpsum;
      r(3,i)=r(3,i)/tmpsum;
    }  
  }
  Matrix b = read_ascii_matrix(opts.bvalsfile.value());
  if(b.Nrows()>1) b=b.t();

  ColumnVector alpha,beta,sinalpha,cosalpha;
  ColumnVector bvals;

  bvals = b.Row(1).t();
  cart2sph(r,alpha,beta);
  sinalpha.ReSize(alpha.Nrows());
  cosalpha.ReSize(alpha.Nrows());
  for(int i=1;i<=alpha.Nrows();i++){
    sinalpha(i) = sin(alpha(i));
    cosalpha(i) = cos(alpha(i));
  }

  // for dti
  Matrix Amat;
  Amat = form_Amat(r,b);
  Amat = pinv(Amat);

  volume4D<float> data;
  volume<int> mask;

  if(opts.verbose.value()) cout<<"reading data"<<endl;
  read_volume4D(data,opts.datafile.value());

  if(opts.verbose.value()) cout<<"reading mask"<<endl;
  read_volume(mask,opts.maskfile.value());

  if(opts.verbose.value()) cout<<"ok"<<endl;
  int minx=0;
  int maxx=mask.xsize();
  int miny=0;
  int maxy=mask.ysize();
  int minz=0;
  int maxz=mask.zsize();
  cout<<minx<<" "<<maxx<<" "<<miny<<" "<<maxy<<" "<<minz<<" "<<maxz<<endl;

  if(opts.verbose.value()) cout<<"setting up vols"<<endl;
  volume<float> S0(maxx-minx,maxy-miny,maxz-minz);
  volume<float> dvol(maxx-minx,maxy-miny,maxz-minz);
  volume<float> tmpvol(maxx-minx,maxy-miny,maxz-minz);
  volume4D<float> tmpvol4D(maxx-minx,maxy-miny,maxz-minz,3);

  vector< volume<float> > fvol,thvol,phvol;
  vector< volume4D<float> > dyads;

  if(opts.verbose.value()) cout<<"copying input properties to output volumes"<<endl;
  copybasicproperties(data[0],S0);
  copybasicproperties(data[0],dvol);

  tmpvol = 0;
  tmpvol4D = 0;
  for(int i=0;i<opts.nfibres.value();i++){
    fvol.push_back(tmpvol);
    thvol.push_back(tmpvol);
    phvol.push_back(tmpvol);
    dyads.push_back(tmpvol4D);
  }

  if(opts.verbose.value()) cout<<"zeroing output volumes"<<endl;
  S0=0;dvol=0;

  if(opts.verbose.value()) cout<<"ok"<<endl;

  //int counter=0;
  ColumnVector S(bvals.Nrows());
  if(opts.verbose.value()) cout<<"starting the fits"<<endl;
  for(int k = minz; k < maxz; k++){
    cout<<k<<" slices processed"<<endl;
    for(int j=miny; j < maxy; j++){
      for(int i =minx; i< maxx; i++){
	if(mask(i,j,k)==0)continue;

	//counter++;
	//OUT(counter);
	//if(counter<1443)continue;

	for(int t=0;t < data.tsize();t++)
	  S(t+1)=data(i,j,k,t);
	
	// initialisation ///////////////////////////////////
	DTI dti(S,Amat);
	dti.fit();

	//dti.print();

	float th,ph;
	cart2sph(dti.get_v1(),th,ph);
	ColumnVector start(2+3*opts.nfibres.value());
	start(1) = dti.get_s0();
	start(2) = (dti.get_md()>0?dti.get_md():0.001);
	start(3) = tan(dti.get_fa())/two_pi;
	start(4) = th;
	start(5) = ph;
	float sumf=abs(two_pi*atan(start(3)));
	float tmpsumf = sumf;
	for(int ff=2;ff<=opts.nfibres.value();ff++){
	  float denom=2;
	  do{
	    start(ff*3) = tan(abs(two_pi*atan(start(3))/(denom)))/two_pi;
	    denom *= 2;
	    tmpsumf = sumf + abs(two_pi*atan(start(ff*3)));
	  }while(tmpsumf>1);
	  sumf += abs(two_pi*atan(start(ff*3)));
	  cart2sph(dti.get_v(ff),th,ph);
	  start(ff*3+1) = th;
	  start(ff*3+2) = ph;
	}
	//////////////////////////////////////////////////////

	PVMNonlinCF pvm_cf(S,bvals,sinalpha,cosalpha,beta,opts.nfibres.value());
	ColumnVector final_par;

	//cout << i << " " << j << " " << k << endl;
	//OUT(start.t());

	if(dti.get_fa()<1 & dti.get_md()>0){
	  NonlinParam  lmpar(start.Nrows(),NL_LM,start); // Levenberg-Marquardt
	  lmpar.SetGaussNewtonType(LM_L); // Levenberg
	  lmpar.SetStartingEstimate(start);
	
	  NonlinOut status = nonlin(lmpar,pvm_cf);
	  final_par = lmpar.Par();
	}
	else{
	  final_par = start;
	}

 	//OUT(final_par.t());

	
//  	OUT(pvm_cf.forwardModel(final_par).t());

	S0(i-minx,j-miny,k-minz)=final_par(1);
	dvol(i-minx,j-miny,k-minz)=final_par(2);
	
	for(int f=0;f<opts.nfibres.value();f++){
	  fvol[f](i-minx,j-miny,k-minz) = abs(2.0/M_PI*atan(final_par(3*(1+f))));
	  thvol[f](i-minx,j-miny,k-minz) = final_par(3*(1+f)+1);
	  phvol[f](i-minx,j-miny,k-minz) = final_par(3*(1+f)+2);
	}
	
	// where we reorder the PVFs
	if(opts.nfibres.value()>1){
	  vector< pair<float,int> > fpairs;
	  vector< pair<float,float> > apairs;
	  pair<float,int> mypair;
	  pair<float,float> mypair2;
	  for(int f=0;f<opts.nfibres.value();f++){
	    mypair.first = fvol[f](i-minx,j-miny,k-minz);
	    mypair.second = f;
	    fpairs.push_back(mypair);
	    mypair2.first = thvol[f](i-minx,j-miny,k-minz);
	    mypair2.second = phvol[f](i-minx,j-miny,k-minz);
	    apairs.push_back(mypair2);
	  } 
	  sort(fpairs.begin(),fpairs.end());
	  for(int f=0,ff=opts.nfibres.value()-1;f<opts.nfibres.value();f++,ff--){
	    fvol[ff](i-minx,j-miny,k-minz) = fpairs[f].first;
	    thvol[ff](i-minx,j-miny,k-minz) = apairs[fpairs[f].second].first;
	    phvol[ff](i-minx,j-miny,k-minz) = apairs[fpairs[f].second].second;
	  }
	}

	// create dyads?
	for(int f=0;f<opts.nfibres.value();f++){
	  dyads[f](i-minx,j-miny,k-minz,0) = sin(thvol[f](i-minx,j-miny,k-minz)) * cos(phvol[f](i-minx,j-miny,k-minz));
	  dyads[f](i-minx,j-miny,k-minz,1) = sin(thvol[f](i-minx,j-miny,k-minz)) * sin(phvol[f](i-minx,j-miny,k-minz));
	  dyads[f](i-minx,j-miny,k-minz,2) = cos(thvol[f](i-minx,j-miny,k-minz));
	}

      }
    }
  }
  
  
  if(opts.verbose.value())
    cout << "saving results" << endl;

  S0.setDisplayMaximumMinimum(S0.max(),0);
  save_volume(S0,opts.ofile.value()+"_S0");

  dvol.setDisplayMaximumMinimum(dvol.max(),0);
  save_volume(dvol,opts.ofile.value()+"_D");

  for(int f=1;f<=opts.nfibres.value();f++){
    fvol[f-1].setDisplayMaximumMinimum(1,0);
    save_volume(fvol[f-1],opts.ofile.value()+"_f"+num2str(f));
    thvol[f-1].setDisplayMaximumMinimum(thvol[f-1].max(),thvol[f-1].min());
    save_volume(thvol[f-1],opts.ofile.value()+"_th"+num2str(f));
    phvol[f-1].setDisplayMaximumMinimum(phvol[f-1].max(),phvol[f-1].min());
    save_volume(phvol[f-1],opts.ofile.value()+"_ph"+num2str(f));
    dyads[f-1].setDisplayMaximumMinimum(-1,1);
    save_volume4D(dyads[f-1],opts.ofile.value()+"_dyads"+num2str(f));

  }

  return 0;
}













