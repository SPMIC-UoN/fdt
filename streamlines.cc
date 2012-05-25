#include "streamlines.h"
#include "warpfns/fnirt_file_reader.h"
#include "warpfns/warpfns.h"

namespace TRACT{

  void read_ascii_files(const string& filename,vector<string>& content){
    ifstream fs(filename.c_str());
    string tmp;
    if(fs){
      fs>>tmp;
      do{
	content.push_back(tmp);
	fs>>tmp;
      }while(!fs.eof());
    }
    else{
      cerr<<"TRACT::read_ascii_files: "<<filename<<" does not exist"<<endl;
      exit(0);
    }
  }

  ColumnVector mean_sph_pol(ColumnVector& A, ColumnVector& B){
    // A and B contain th, ph f. 
    float th,ph;
    ColumnVector rA(3), rB(3);
    
    rA << (sin(A(1))*cos(A(2))) << (sin(A(1))*sin(A(2))) << (cos(A(1)));
    rB << (sin(B(1))*cos(B(2))) << (sin(B(1))*sin(B(2))) << (cos(B(1)));
    
    if(sum(SP(rA,rB)).AsScalar()>0)
      cart2sph((rA+rB)/2,th,ph);
    else
      cart2sph((rA-rB)/2,th,ph);
    
    A(1)=th; A(2)=ph;
    return A;
  }
  
  void imgradient(const volume<float>& im,volume4D<float>& grad){
    
    grad.reinitialize(im.xsize(),im.ysize(),im.zsize(),3);
    copybasicproperties(im,grad[0]);
    
    int fx,fy,fz,bx,by,bz;
    float dx,dy,dz; 
    for (int z=0; z<grad.zsize(); z++){
      fz = z ==(grad.zsize()-1) ? 0 :  1;
      bz = z == 0              ? 0 : -1;
      dz = (fz==0 || bz==0)    ? 1.0 :  2.0;
      for (int y=0; y<grad.ysize(); y++){
	fy = y ==(grad.ysize()-1) ? 0 :  1;
	by = y == 0              ? 0 : -1;
	dy = (fy==0 || by==0)    ? 1.0 :  2.0;
	for (int x=0; x<grad.xsize(); x++){
	  fx = x ==(grad.xsize()-1) ? 0 :  1;
	  bx = x == 0              ? 0 : -1;
	  dx = (fx==0 || bx==0)    ? 1.0 :  2.0;
	  grad[0](x,y,z) = (im(x+fx,y,z) - im(x+bx,y,z))/dx;
	  grad[1](x,y,z) = (im(x,y+fy,z) - im(x,y+by,z))/dy;
	  grad[2](x,y,z) = (im(x,y,z+fz) - im(x,y,z+bz))/dz;
	}
      }
    }
    
  }
  
  void make_unique(vector<int>& x,vector<float>& y){
    // sort
    vector< pair<int,int> > p;
    for(unsigned int i=0;i<x.size();i++){
      pair<int,int> pp;
      pp.first=x[i];
      pp.second=i;
      p.push_back(pp);
    }
    sort(p.begin(),p.end());
    vector<float> ycp(y.begin(),y.end());
    for(unsigned int i=0;i<x.size();i++){
      x[i]=p[i].first;
      y[i]=ycp[ p[i].second ];
    }
    // unique
    vector<int> xx;vector<float> yy;
    for(unsigned int i=0;i<x.size();i++){
      if(i>0){
	if( x[i]==x[i-1] )continue;
      }
      xx.push_back(x[i]);
      yy.push_back(y[i]);
    }
    x=xx;
    y=yy;
  }

  void make_unique(vector< pair<int,float> >&x){
    sort(x.begin(),x.end());
    vector< pair<int,float> > xx;
    for(unsigned int i=0;i<x.size();i++){
      if(i>0){
	if( x[i].first==x[i-1].first )continue;
      }
      xx.push_back(x[i]);
    }
    x=xx;
  }

  Streamliner::Streamliner(const CSV& seeds):opts(probtrackxOptions::getInstance()),
					     logger(LogSingleton::getInstance()),					     					     
					     vols(opts.usef.value()),
					     m_seeds(seeds)
  {    

    // the tracking mask - no tracking outside of this
    read_volume(m_mask,opts.maskfile.value());
    
    // particle
    m_part.initialise(0,0,0,0,0,0,opts.steplength.value(),
		      m_mask.xdim(),m_mask.ydim(),m_mask.zdim(),
		      false);
    
    // no-integrity
    if(opts.skipmask.value()!="") 
      read_volume(m_skipmask,opts.skipmask.value());
    
    // loop checking box
    m_lcrat=5; // hard-coded! box is five times smaller than brain mask
    if(opts.loopcheck.value()){
      m_loopcheck.reinitialize(int(ceil(m_mask.xsize()/m_lcrat)+1),
			       int(ceil(m_mask.ysize()/m_lcrat)+1),
			       int(ceil(m_mask.zsize()/m_lcrat)+1),3);
      m_loopcheck=0;
    }
    
    // exclusion mask
    // now in CSV format
    if(opts.rubbishfile.value()!=""){
      load_rubbish(opts.rubbishfile.value());   
    }
    
    // stopping mask
    // now in CSV format
    if(opts.stopfile.value()!=""){
      load_stop(opts.stopfile.value());
    }
     
    // waymasks in CSV format
    if(opts.waypoints.set()){
      load_waymasks(opts.waypoints.value());      
      m_waycond=opts.waycond.value();
    }
    
    // choose preferred fibre within this mask (prior knowledge)
    if(opts.prefdirfile.value()!=""){
      read_volume4D(m_prefdir,opts.prefdirfile.value());
    }    
    
    // Allow for either matrix transform (12dof affine) or nonlinear (warpfield)
    m_Seeds_to_DTI = IdentityMatrix(4);
    m_DTI_to_Seeds = IdentityMatrix(4);
    m_rotdir       = IdentityMatrix(3);

    m_IsNonlinXfm = false;
    if(opts.seeds_to_dti.value()!=""){
      if(!fsl_imageexists(opts.seeds_to_dti.value())){//presumably ascii file provided
	m_Seeds_to_DTI = read_ascii_matrix(opts.seeds_to_dti.value());
	m_DTI_to_Seeds = m_Seeds_to_DTI.i();
	m_rotdir       = m_Seeds_to_DTI.SubMatrix(1,3,1,3);
      }
      else{
	m_IsNonlinXfm = true;
	FnirtFileReader ffr(opts.seeds_to_dti.value());
	m_Seeds_to_DTI_warp = ffr.FieldAsNewimageVolume4D(true);
	if(opts.dti_to_seeds.value()==""){
	  cerr << "TRACT::Streamliner:: DTI -> Seeds transform needed" << endl;
	  exit(1);
	}
	FnirtFileReader iffr(opts.dti_to_seeds.value());
	m_DTI_to_Seeds_warp = iffr.FieldAsNewimageVolume4D(true);

	// now calculate the jacobian of this transformation (useful for rotating vectors)
	imgradient(m_Seeds_to_DTI_warp[0],m_jacx);
	imgradient(m_Seeds_to_DTI_warp[1],m_jacy);
	imgradient(m_Seeds_to_DTI_warp[2],m_jacz);
      }
    }
    
    vols.initialise(opts.basename.value(),m_mask);
    m_path.reserve(opts.nsteps.value());
    m_x_s_init=0;
    m_y_s_init=0;
    m_z_s_init=0;

    m_inmask3.reserve(opts.nsteps.value());
    m_inlrmask3.reserve(opts.nsteps.value());
    
  }
  
  void Streamliner::rotdir(const ColumnVector& dir,ColumnVector& rotdir,
			   const float& x,const float& y,const float& z){
    if(!m_IsNonlinXfm){
      rotdir = m_rotdir*dir;
      if(rotdir.MaximumAbsoluteValue()>0)
	rotdir = rotdir/std::sqrt(rotdir.SumSquare());
    }
    else{
      ColumnVector xyz_dti(3),xyz_seeds(3);
      xyz_seeds << x << y << z;
      xyz_dti = NewimageCoord2NewimageCoord(m_DTI_to_Seeds_warp,
					    false,m_seeds.get_refvol(),m_mask,xyz_seeds);

      Matrix F(3,3),Jw(3,3);
      int x=(int)round((float)xyz_dti(1));
      int y=(int)round((float)xyz_dti(2));
      int z=(int)round((float)xyz_dti(3));
      Jw << m_jacx(x,y,z,0) << m_jacx(x,y,z,1) << m_jacx(x,y,z,2)
	 << m_jacy(x,y,z,0) << m_jacy(x,y,z,1) << m_jacy(x,y,z,2)
	 << m_jacz(x,y,z,0) << m_jacz(x,y,z,1) << m_jacz(x,y,z,2);
      F = (IdentityMatrix(3) + Jw).i();
      rotdir = F*dir;
      if(rotdir.MaximumAbsoluteValue()>0)
	rotdir = rotdir/std::sqrt(rotdir.SumSquare());

    }
  }
  
  int Streamliner::streamline(const float& x_init,const float& y_init,const float& z_init, 
			      const ColumnVector& dim_seeds,const int& fibst,const int& seed_ind){ 
    //Tracer_Plus tr("Streamliner::streamline");
    //fibst tells tractvolsx which fibre to start with if there are more than one..
    //x_init etc. are in seed space...
    vols.reset(fibst);
    m_x_s_init=x_init; // seed x position in voxels
    m_y_s_init=y_init; // and y
    m_z_s_init=z_init; // and z
    ColumnVector xyz_seeds(3);
    xyz_seeds<<x_init<<y_init<<z_init;
    ColumnVector xyz_dti;
    ColumnVector th_ph_f;
    float xst,yst,zst,x,y,z,tmp2;
    float pref_x=0,pref_y=0,pref_z=0;
    int x_s,y_s,z_s;

    // find xyz in dti space
    if(!m_IsNonlinXfm)
      xyz_dti = vox_to_vox(xyz_seeds,dim_seeds,vols.dimensions(),m_Seeds_to_DTI);    
    else{
      xyz_dti = NewimageCoord2NewimageCoord(m_DTI_to_Seeds_warp,false,m_seeds.get_refvol(),m_mask,xyz_seeds);
    }

    xst=xyz_dti(1);yst=xyz_dti(2);zst=xyz_dti(3);
    m_path.clear();
    x=xst;y=yst;z=zst;
    m_part.change_xyz(x,y,z);


    float          pathlength=0;
    bool           rubbish_passed=false;
    bool           wayorder=true;
    vector<int>    waycrossed;
    vector<int>    crossedlocs3;
    int            cnt=-1;
    // loopcheck stuff
    float oldrx,oldry,oldrz;
    int   lcx,lcy,lcz;
    bool forcedir=false;
    //NB - this only goes in one direction!!

    if(opts.onewaycondition.value()){
      for(unsigned int i=0; i<m_way_passed_flags.size();i++)
	m_way_passed_flags[i]=0;
      if(opts.network.value()){
	m_net_passed_flags=0;
      }
    }


    for(int it=1;it<=opts.nsteps.value()/2;it++){
      
      if((m_mask((int)round(m_part.x()),(int)round(m_part.y()),(int)round(m_part.z()))!=0)){

	///////////////////////////////////
	//loopchecking
	///////////////////////////////////
	if(opts.loopcheck.value()){
	  lcx=(int)round(m_part.x()/m_lcrat);
	  lcy=(int)round(m_part.y()/m_lcrat);
	  lcz=(int)round(m_part.z()/m_lcrat);
	  oldrx=m_loopcheck(lcx,lcy,lcz,0);
	  oldry=m_loopcheck(lcx,lcy,lcz,1);
	  oldrz=m_loopcheck(lcx,lcy,lcz,2);
	  if(m_part.rx()*oldrx+m_part.ry()*oldry+m_part.rz()*oldrz<0){break;}
	  m_loopcheck(lcx,lcy,lcz,0)=m_part.rx();
	  m_loopcheck(lcx,lcy,lcz,1)=m_part.ry();
	  m_loopcheck(lcx,lcy,lcz,2)=m_part.rz();
	}
	
	if(opts.verbose.value()>1)
	  logger<<m_part;
	
	x=m_part.x();y=m_part.y();z=m_part.z();
	xyz_dti <<x<<y<<z;

	// now find xyz in seeds space
	if(cnt>=0){
	  if(!m_IsNonlinXfm)
	    xyz_seeds = vox_to_vox(xyz_dti,vols.dimensions(),dim_seeds,m_DTI_to_Seeds);    
	  else{
	    xyz_seeds = NewimageCoord2NewimageCoord(m_Seeds_to_DTI_warp,false,m_mask,m_seeds.get_refvol(),xyz_dti);
	  }
	}
	
	x_s =(int)round((float)xyz_seeds(1));
	y_s =(int)round((float)xyz_seeds(2));
	z_s =(int)round((float)xyz_seeds(3));
	

	// how prefdir works:
	// if a vector field is given (tested using 4th dim==3) then set prefdir to that orientation
	// eventually will be flipped to match the current direction of travel
	// if dim4<3, then :
	//    - use dim1 scalar as fibre sampling method (i.e. randfib 1 2 3 )
	//    - use dim2 scalar (when given) as which fibre to follow blindly, i.e. without flipping
	pref_x=0;pref_y=0;pref_z=0;
	int sample_fib = 0;
	if(opts.prefdirfile.value()!=""){
	  if(m_prefdir.tsize()==3){
	    pref_x = m_prefdir((int)x,(int)y,(int)z,0);
	    pref_y = m_prefdir((int)x,(int)y,(int)z,1);
	    pref_z = m_prefdir((int)x,(int)y,(int)z,2);
	  }
	  else{
	    if(m_prefdir(x_s,y_s,z_s,0)!=0)
	      sample_fib = (int)m_prefdir((int)x,(int)y,(int)z,0);
	    if(sample_fib>3)sample_fib=3;	    
	  }
	}

	// // attractor mask 
	// if(opts.pullmask.value()!=""){
	//   ColumnVector pulldir(3);float mindist;
	//   if(m_pullmask.is_near_surface(xyz_seeds,opts.pulldist.value(),pulldir,mindist)){
	//     m_prefx = pulldir(1);
	//     m_prefy = pulldir(2);
	//     m_prefz = pulldir(3);
	//   }	  
	// }


	// augment the path	
	m_path.push_back(xyz_seeds);cnt++;
	if(cnt>0)
	  pathlength += opts.steplength.value();


	// // // if first step and onewayonly
// // 	if(opts.onewayonly.value() && m_path.size()==2){
// // 	  Vec step(m_path[1](1)-m_path[0](1),
// // 		   m_path[1](2)-m_path[0](2),
// // 		   m_path[1](3)-m_path[0](3));

// // 	  if(m_seeds.coord_sign(seed_ind,m_path[0]+0.001*(m_path[1]-m_path[0]))<0){
// // 	    rubbish_passed=1;
// // 	    break;
// // 	  }	  
// // 	}

	
	
	// only test exclusion after at least one step
	if(opts.rubbishfile.value()!="" && cnt>0){
	  if(m_rubbish.has_crossed(m_path[cnt-1],m_path[cnt])){
	    rubbish_passed=1;	
	    break;
	  }
	}

	// update every passed_flag
	if(m_way_passed_flags.size()>0){
	  if(cnt>0){
	    waycrossed.clear();
	    m_waymasks.has_crossed_roi(m_path[cnt-1],m_path[cnt],waycrossed);

	    for(unsigned int wm=0;wm<waycrossed.size();wm++){
	      m_way_passed_flags[waycrossed[wm]]=1;
	    }
	    // check if order is respected
	    if(opts.wayorder.value() && wayorder){
	      for(unsigned int wm=0;wm<waycrossed.size();wm++){
		for(int wm2=0;wm2<waycrossed[wm];wm2++){
		  if(m_way_passed_flags[wm2]==0){ wayorder=false; break;}
		}
		if(!wayorder)break;
	      }
	    }	    
	  }
	}
	if(opts.network.value()){
	  if(cnt>0){
	    waycrossed.clear();
	    m_netmasks.has_crossed_roi(m_path[cnt-1],m_path[cnt],waycrossed);

	    for(unsigned int wm=0;wm<waycrossed.size();wm++){
	      m_net_passed_flags(waycrossed[wm]+1)=1;
	    }
	  }	  
	}

	// //////////////////////////////
	// update locations for matrix3
	if(opts.matrix3out.value() && cnt>0){
	  waycrossed.clear();crossedlocs3.clear();
	  if(m_mask3.has_crossed_roi(m_path[cnt-1],m_path[cnt],waycrossed,crossedlocs3)){	    
	    fill_inmask3(crossedlocs3,pathlength);
	  }
	  if(opts.lrmask3.value()!=""){
	    waycrossed.clear();crossedlocs3.clear();
	    if(m_lrmask3.has_crossed_roi(m_path[cnt-1],m_path[cnt],waycrossed,crossedlocs3)){	    
	      fill_inlrmask3(crossedlocs3,pathlength);
	    }	  
	  }
	}




	// only test stopping after at least one step
	if(opts.stopfile.value()!="" && m_path.size()>1){
	  if(m_path.size()==2 && opts.forcefirststep.value()){
	    // do nothing
	  }
	  else if(m_stop.has_crossed(m_path[cnt-1],m_path[cnt])){
	    break;	    
	  }	  
	}


	// //////////////////////////////

	// sample a new fibre orientation
	int sampled_fib,newx,newy,newz;	
	if(opts.skipmask.value() == ""){
	  //Tracer_Plus tr("sample");
	  th_ph_f = vols.sample(m_part.x(),m_part.y(),m_part.z(),    // sample at this location
				m_part.rx(),m_part.ry(),m_part.rz(), // choose closest sample to this
				pref_x,pref_y,pref_z,                // unless we have this prefered direction 
				sample_fib,                          // controls probabilities of sampling a given fibre
				sampled_fib,                         // which fibre has been sampled
				newx,newy,newz);                     // which voxel has been considered
	}
	else{
	  if(m_skipmask(x_s,y_s,z_s)==0)//only sample outside skipmask
	    th_ph_f=vols.sample(m_part.x(),m_part.y(),m_part.z(),
				m_part.rx(),m_part.ry(),m_part.rz(),
				pref_x,pref_y,pref_z,
				sample_fib,
				sampled_fib,
				newx,newy,newz);
	}
	

// 	// depending on which fibre has been sampled, decide whether to forcedir
// 	forcedir=false;
// 	if(opts.prefdirfile.value()!=""){
// 	  if(m_prefdir.tsize()==2){
// 	    if(m_prefdir(newx,newy,newz,1)==sampled_fib)
// 	      forcedir=true;
// 	  }
// 	}
 
	// jump
	tmp2=(float)rand()/(float)RAND_MAX;	

	if(th_ph_f(3)>tmp2){ //volume fraction criterion	  
	  if(!m_part.check_dir(th_ph_f(1),th_ph_f(2),opts.c_thr.value(),forcedir)){// curvature threshold
	    //cout<<"curv break"<<endl;
	    break;
	  }

	  if((th_ph_f(1)!=0 && th_ph_f(2)!=0)){
	    if( (m_mask((int)round(m_part.x()),(int)round(m_part.y()),(int)round(m_part.z())) != 0) ){

	      if(!opts.modeuler.value())
		m_part.jump(th_ph_f(1),th_ph_f(2),forcedir);

	      else{
		  ColumnVector test_th_ph_f;		  
		  m_part.testjump(th_ph_f(1),th_ph_f(2),forcedir);

		  test_th_ph_f=vols.sample(m_part.testx(),m_part.testy(),m_part.testz(),
					   m_part.rx(),m_part.ry(),m_part.rz(),
					   pref_x,pref_y,pref_z,sample_fib,sampled_fib,newx,newy,newz);
		  test_th_ph_f=mean_sph_pol(th_ph_f,test_th_ph_f);
		  m_part.jump(test_th_ph_f(1),test_th_ph_f(2),forcedir);
		}
	    }
	    else{
	      //cout<<"mask break 1"<<endl;
	      break; // outside mask	    
	    }
	  }
	}
	
	  }
      else{
	//cout<<"mask break 2"<<endl;
	break;
      }// outside mask

      
    } // Close Step Number Loop (done tracking sample)

    // reset loopcheck box
    if(opts.loopcheck.value()){
      m_loopcheck=0;
    }
    
    // rejflag = 0 (accept), 1 (reject) or 2 (wait for second direction)
    int rejflag=0;

    if(rubbish_passed) return(1);
    if(pathlength<opts.distthresh.value()) return(1);

    if(opts.network.value()){
      unsigned int numpassed=0;
      for(int i=1; i<m_net_passed_flags.Nrows();i++){
	if(m_net_passed_flags(i))numpassed++;
      }
      if(numpassed==0)rejflag=1;
      else if((int)numpassed<m_net_passed_flags.Nrows()){
	if(m_waycond=="AND"){
	  rejflag=opts.onewaycondition.value()?1:2;
	}
	else{
	  rejflag=0;
	}
      }
      else rejflag=0;
    }
    if(rejflag==1) return(rejflag);
    
    // Waypoints
    if(m_way_passed_flags.size()!=0){
      unsigned int numpassed=0;
      for(unsigned int i=0; i<m_way_passed_flags.size();i++){
	if(m_way_passed_flags[i]>0)numpassed++;
      }
      
      if(numpassed==0)rejflag=1;
      else if(numpassed<m_way_passed_flags.size()){
	if(m_waycond=="AND"){
	  if(opts.onewaycondition.value()){rejflag=1;}
	  else{
	    if(opts.wayorder.value()){rejflag=wayorder?2:1;}
	    else{rejflag=2;}
	  }
	}
	else{
	  rejflag=0;
	}
      }
      else{
	if(m_waycond=="OR" || !opts.wayorder.value())
	  rejflag=0;
	else
	  rejflag=wayorder?0:1;
      }
    }

    return rejflag;    
  }
  
  void Counter::initialise(){
   
    if(opts.simpleout.value()){
      initialise_path_dist();
    }
    if(opts.s2tout.value()){
      initialise_seedcounts();
    }
    if(opts.matrix1out.value()){
      initialise_matrix1();
    }
    if(opts.matrix2out.value()){
      initialise_matrix2();
    }
    if(opts.matrix3out.value()){
      initialise_matrix3();
    }
  }

  
  void Counter::initialise_seedcounts(){
    // now the CSV class does all the work
    m_targetmasks.reinitialize(m_stline.get_seeds().get_refvol());
    m_targetmasks.set_convention(opts.meshspace.value());
    m_targetmasks.load_rois(opts.targetfile.value());
    m_targetmasks.reset_values();

    // seeds are CSV-format
    if(!opts.simple.value()){
      m_s2t_count=m_stline.get_seeds();
      m_s2t_count.reset_values();
      for(int m=0;m<m_targetmasks.nRois();m++){
	m_s2t_count.add_map();
      }
    }
    // seeds are text
    if(opts.simple.value() || opts.s2tastext.value()){
      m_s2tastext.ReSize(m_numseeds,m_targetmasks.nRois());
      m_s2tastext=0;
    }

    m_s2trow=1;
    m_targflags.resize(m_targetmasks.nRois());
  }
  

  // matrix1 is nseeds X nseeds
  void Counter::initialise_matrix1(){
    m_ConMat1 = new SpMat<float> (m_numseeds,m_numseeds);
    m_Conrow1 = 1;

    vector<ColumnVector> coords = m_stline.get_seeds().get_locs_coords();
    vector<int> roicind         = m_stline.get_seeds().get_locs_coord_index();
    vector<int> roiind          = m_stline.get_seeds().get_locs_roi_index();
    
    Matrix CoordMat1(m_numseeds,5); 
    for (unsigned int i=0;i<coords.size();i++)
      CoordMat1.Row(i+1) << (float)coords[i](1) 
			 << (float)coords[i](2)
			 << (float)coords[i](3)
			 << roiind[i]
			 << roicind[i];
    
    applycoordchange(CoordMat1, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
    write_ascii_matrix(CoordMat1,logger.appendDir("coords_for_fdt_matrix1"));
     
  }
  
  // matrix2 is nseeds X nlrmask
  void Counter::initialise_matrix2(){
    // init lrmask-related 
    read_volume(m_lrmask,opts.lrmask.value());
    m_beenhere2.reinitialize(m_lrmask.xsize(),m_lrmask.ysize(),m_lrmask.zsize());
    m_lookup2.reinitialize(m_lrmask.xsize(),m_lrmask.ysize(),m_lrmask.zsize());
    m_lookup2=0;
    m_lrdim.ReSize(3);
    m_lrdim<<m_lrmask.xdim()<<m_lrmask.ydim()<<m_lrmask.zdim();

    int numnz=0;    
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask.value(Wx,Wy,Wz)!=0){
	    numnz++;
	    m_lookup2(Wx,Wy,Wz)=numnz;
	  }
    
    Matrix CoordMat_tract2(numnz,3);
    int mytrow=1;
    for(int Wz=m_lrmask.minz();Wz<=m_lrmask.maxz();Wz++)
      for(int Wy=m_lrmask.miny();Wy<=m_lrmask.maxy();Wy++)
	for(int Wx=m_lrmask.minx();Wx<=m_lrmask.maxx();Wx++)
	  if(m_lrmask(Wx,Wy,Wz)!=0){
	    CoordMat_tract2(mytrow,1)=Wx;
	    CoordMat_tract2(mytrow,2)=Wy;
	    CoordMat_tract2(mytrow,3)=Wz;
	    mytrow++;
	  }

    applycoordchange(CoordMat_tract2, m_lrmask.niftivox2newimagevox_mat().i());
    write_ascii_matrix(CoordMat_tract2,logger.appendDir("tract_space_coords_for_fdt_matrix2"));
    save_volume(m_lookup2,logger.appendDir("lookup_tractspace_fdt_matrix2"));

    
    // init matrix2-related
    m_ConMat2 = new SpMat<float> (m_numseeds,numnz);

    if( !opts.simple.value()){
      vector<ColumnVector> coords = m_stline.get_seeds().get_locs_coords();
      vector<int> roicind         = m_stline.get_seeds().get_locs_coord_index();
      vector<int> roiind          = m_stline.get_seeds().get_locs_roi_index();

      Matrix CoordMat2(m_numseeds,5);
      for (unsigned int i=0;i<coords.size();i++)
	CoordMat2.Row(i+1) << (float)coords[i](1) 
			   << (float)coords[i](2)
			   << (float)coords[i](3)
			   << roiind[i]
			   << roicind[i];
      
      applycoordchange(CoordMat2, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
      write_ascii_matrix(CoordMat2,logger.appendDir("coords_for_fdt_matrix2"));
    }


  }
  
  // the following will use the voxels of a given mask to initialise 
  // an NxN matrix. This matrix will store the number of samples from 
  // each seed voxel that have made it to each pair of voxels in the mask
  // now mask3 is in CSV format
  void Counter::initialise_matrix3(){
    m_stline.init_mask3();

    int nmask3 = m_stline.get_mask3().nLocs();
    int nlrmask3;
    if(opts.lrmask3.value()!="")
      nlrmask3 = m_stline.get_lrmask3().nLocs();
    else
      nlrmask3 = nmask3;

    // recalculate nmask3 if lowres surface provided
    m_ConMat3  = new SpMat<float> (nmask3,nlrmask3); 
    //OUT(m_ConMat3->Nrows());
    //OUT(m_ConMat3->Ncols());
    //exit(1);

    // save lookup tables...
    
    CSV  mask3(m_stline.get_mask3());
    vector<ColumnVector> coords = mask3.get_locs_coords();
    vector<int> roicind = mask3.get_locs_coord_index();
    vector<int> roiind = mask3.get_locs_roi_index();
    
    Matrix mat(mask3.nLocs(),5);
    for (unsigned int i=0;i<coords.size();i++)
      mat.Row(i+1) << coords[i](1) 
		   << coords[i](2)
		   << coords[i](3)
		   << roiind[i]
		   << roicind[i];
    
    applycoordchange(mat, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
    write_ascii_matrix(mat,logger.appendDir("coords_for_fdt_matrix3"));

    if(opts.lrmask3.value()!=""){
      CSV lrmask3(m_stline.get_lrmask3());
      vector<ColumnVector> lrcoords = lrmask3.get_locs_coords();
      vector<int> lrroicind = lrmask3.get_locs_coord_index();
      vector<int> lrroiind = lrmask3.get_locs_roi_index();
    
      Matrix lrmat(lrmask3.nLocs(),5);
      for (unsigned int i=0;i<lrcoords.size();i++)
	lrmat.Row(i+1) << lrcoords[i](1) 
		     << lrcoords[i](2)
		     << lrcoords[i](3)
		     << lrroiind[i]
		     << lrroicind[i];
    
      applycoordchange(lrmat, m_stline.get_seeds().get_refvol().niftivox2newimagevox_mat().i());
      write_ascii_matrix(lrmat,logger.appendDir("tract_space_coords_for_fdt_matrix3"));
    }

  }

  void Counter::count_streamline(){
    if(opts.save_paths.value()){
      add_path();
    }
    if(opts.simpleout.value()||opts.matrix1out.value()){
      update_pathdist();
    }
    if(opts.s2tout.value()){
      update_seedcounts();
    }
    if(opts.matrix2out.value()){
      update_matrix2_row();
    }
  }

  void Counter::count_seed(){
    if(opts.s2tout.value()){
      m_s2trow++;
    }
    if(opts.matrix1out.value()){
      m_Conrow1++;;
    }
  }

  void Counter::clear_streamline(){
    if(opts.simpleout.value()||opts.matrix1out.value()){
      reset_beenhere();
    }
    if(opts.s2tout.value()){
      reset_targetflags();
    }
    if(opts.matrix2out.value()){
      reset_beenhere2();
    }
    if(opts.matrix3out.value()){
      reset_beenhere3();
    }
  }

  void Counter::update_pathdist(){
    int x_s,y_s,z_s;
    float pathlength=0;
    vector<int> crossedrois,crossedlocs;
    int nlocs=0;
    //OUT(m_path.size());
    for(unsigned int i=0;i<m_path.size();i++){
      //OUT(i);
      //OUT(pathlength);
      //OUT(m_path[i].t());
      x_s=(int)round((float)m_path[i](1));
      y_s=(int)round((float)m_path[i](2));
      z_s=(int)round((float)m_path[i](3));
      // check here if back to seed
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	//cout<<"ste pathlength to zero"<<endl;
	pathlength=0;
      }

      if(m_beenhere(x_s,y_s,z_s)==0){
	if(!opts.pathdist.value())
	  m_prob(x_s,y_s,z_s)+=1; 
	else
	  m_prob(x_s,y_s,z_s)+=pathlength;
	m_beenhere(x_s,y_s,z_s)=1;

	if(opts.opathdir.value() && i>0){
	  ColumnVector v(3);
	  v=m_path[i]-m_path[i-1];
	  v/=std::sqrt(v.SumSquare());
	  m_localdir(x_s,y_s,z_s,0)+=v(1);
	  m_localdir(x_s,y_s,z_s,1)+=v(2);
	  m_localdir(x_s,y_s,z_s,2)+=v(3);
	}
      }

      // Fill alternative mask
      // This mask's values are:
      //  0: location not to be considered
      //  1: location has not been visited 
      //  2: location has been visited

      if(opts.pathfile.set()){
	//cout<<"here"<<endl;
	if(pathlength>0){
	  //OUT(m_path[i-1].t());OUT(m_path[i].t());
	  crossedrois.clear();crossedlocs.clear();
	  if(m_beenhere_alt.has_crossed_roi(m_path[i-1],m_path[i],crossedrois,crossedlocs)){
	    //cout<<"update"<<endl;
	    nlocs+=crossedlocs.size();
	    for(unsigned int i=0;i<crossedlocs.size();i++){
	      //OUT(m_beenhere_alt.get_value(crossedlocs[i]));
	      if(m_beenhere_alt.get_value(crossedlocs[i])==1){	      
		if(!opts.pathdist.value())
		  m_prob_alt.add_value(crossedlocs[i],1); 
		else
		  m_prob_alt.add_value(crossedlocs[i],pathlength);
		m_beenhere_alt.set_value(crossedlocs[i],2);
	      }
	    }
	  }
	}
      }
      pathlength+=opts.steplength.value();
      
    } 
    //cout<<"update"<<endl;
    //OUT(nlocs);
  }

  void Counter::reset_beenhere(){
    int x_s,y_s,z_s;
    vector<int> crossedlocs,crossedrois;
    bool hascrossed=false;
    for(unsigned int i=0;i<m_path.size();i++){
      x_s=(int)round((float)m_path[i](1));
      y_s=(int)round((float)m_path[i](2));
      z_s=(int)round((float)m_path[i](3));
      m_beenhere(x_s,y_s,z_s)=0;

      // back to first point? keep going
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0)
	continue;

      // alternative roi
      if(opts.pathfile.set()){
	if(i>0){
	  if(m_beenhere_alt.has_crossed_roi(m_path[i-1],m_path[i],crossedrois,crossedlocs))
	    hascrossed=true;
	}
      }
    }
    if(opts.pathfile.set() && hascrossed){
      //cout<<"reset"<<endl;
      //OUT(crossedlocs.size());
      for(unsigned int i=0;i<crossedlocs.size();i++){
	if(m_beenhere_alt.get_value(crossedlocs[i])>1)
	  m_beenhere_alt.set_value(crossedlocs[i],1);
      }      
    }

  }

  void Counter::add_path(){
    m_save_paths.push_back(m_path);
  }
  void Counter::save_paths(){
    string filename=logger.appendDir("saved_paths.txt");
    ofstream of(filename.c_str());
    if (of.is_open()){
      for(unsigned int i=0;i<m_save_paths.size();i++){
	stringstream flot;
	flot << "# " << m_save_paths[i].size()<<endl;
	for(unsigned int j=0;j<m_save_paths[i].size();j++){
	  flot << m_save_paths[i][j](1) << " " << m_save_paths[i][j](2) << " " << m_save_paths[i][j](3) << endl;
	}
	of<<flot.str();
      }
      of.close();
    }
    else{
      cerr<<"Counter::save_paths:error opening file for writing: "<<filename<<endl;
    }
  }
  

  void Counter::update_seedcounts(){
    vector<int> crossed;
    float       pathlength=0;
    int         cnt=0;
    for(unsigned int i=1;i<m_path.size();i++){
      pathlength+=opts.steplength.value();
      // check here if back to seed
      if((m_path[i]-m_path[0]).MaximumAbsoluteValue()==0)
	pathlength=0;

      // this would be more efficient if we only checked 
      // masks that haven't been crossed yet...
      crossed.clear();
      m_targetmasks.has_crossed_roi(m_path[i-1],m_path[i],crossed);
      for(unsigned int t=0;t<crossed.size();t++){
	if(m_targflags[crossed[t]])continue;

	if(!opts.simple.value()){
	  if(!opts.pathdist.value()){
	    m_s2t_count.add_map_value(m_curloc,1,crossed[t]);
	  }
	  else
	    m_s2t_count.add_map_value(m_curloc,pathlength,crossed[t]);
	}
	if(opts.simple.value() || opts.s2tastext.value()){
	  if(!opts.pathdist.value())
	    m_s2tastext(m_s2trow,crossed[t]+1)+=1;
	  else
	    m_s2tastext(m_s2trow,crossed[t]+1)+=pathlength;
	}

	m_targflags[crossed[t]]=true;cnt++;

      }
      // stop if they all have been crossed
      if(cnt==m_targetmasks.nRois())break;
    }
  }

  void Counter::update_matrix1(){
    // use path and has_crossed
    float pathlength=opts.steplength.value(),val=1;
    vector<int> crossedseeds,crossedlocs;
    vector<int> allcrossed;vector<float> allvals;
    for(unsigned int i=1;i<m_path.size();i++){
      // check here if back to seed (not a robust way to do it...)
      if((m_path[i]-m_path[0]).MaximumAbsoluteValue()==0){
	pathlength=opts.steplength.value();
	val = opts.pathdist.value()?pathlength:1;
	continue;
      }
      crossedseeds.clear();crossedlocs.clear();

      if(m_stline.get_seeds().has_crossed_roi(m_path[i-1],m_path[i],crossedseeds,crossedlocs)){
	allcrossed.insert(allcrossed.end(),crossedlocs.begin(),crossedlocs.end());
	vector<float> vals(crossedlocs.size(),val);
	allvals.insert(allvals.end(),vals.begin(),vals.end());
      }
      
      val = opts.pathdist.value()?pathlength:1;      
      pathlength+=opts.steplength.value();
    }
    // fill matrix1
    make_unique(allcrossed,allvals);
    for(unsigned int i=0;i<allcrossed.size();i++){
      m_ConMat1->AddTo(m_Conrow1,allcrossed[i]+1,allvals[i]);
    }
  }
  
  void Counter::update_matrix2_row(){
    //run this one every streamline - not every voxel..
    float d=opts.steplength.value();
    int x_lr,y_lr,z_lr,Concol2;
    ColumnVector xyz_lr(3);
    for(unsigned int i=0;i<m_path.size();i++){
      // check here if back to seed
      if(i>0 && (m_path[i]-m_path[0]).MaximumAbsoluteValue()==0)
	d=opts.steplength.value();

      xyz_lr=vox_to_vox(m_path[i],m_seedsdim,m_lrdim,m_I);	
      x_lr=(int)round((float)xyz_lr(1));
      y_lr=(int)round((float)xyz_lr(2));
      z_lr=(int)round((float)xyz_lr(3));
      Concol2=m_lookup2(x_lr,y_lr,z_lr);

      if(Concol2>0){
	if(m_beenhere2(x_lr,y_lr,z_lr)==0){	  
	  if(!opts.pathdist.value())
	    m_ConMat2->AddTo(m_curloc+1,Concol2,1);
	  else
	    m_ConMat2->AddTo(m_curloc+1,Concol2,d);
	  m_beenhere2(x_lr,y_lr,z_lr)=1;
	  d+=opts.steplength.value();
	}
      }
    } 
  }

  void Counter::update_matrix3(){
    vector< pair<int,float> >& inmask3 = m_stline.get_inmask3();
    vector< pair<int,float> >& inlrmask3 = m_stline.get_inlrmask3();
    bool uselr = (opts.lrmask3.value()!="");
    if(!uselr){
      if(inmask3.size()<2)return;
    }
    else{
      if(inmask3.size()<1 || inlrmask3.size()<1)return;
    }

    // remove duplicates
    make_unique(inmask3);

    if(!uselr){// where we update NxN matrix
      for(unsigned int i=0;i<inmask3.size();i++){
	for(unsigned int j=i+1;j<inmask3.size();j++){
	  if(!opts.pathdist.value())
	    m_ConMat3->AddTo(inmask3[i].first+1,inmask3[j].first+1,1);
	  else{
	    float val = fabs(inmask3[i].second-inmask3[j].second);
	    m_ConMat3->AddTo(inmask3[i].first+1,inmask3[j].first+1,val);
	  }
	}
      }
    }
    else{ // where we update Nxn matrix
      make_unique(inlrmask3);
      for(unsigned int i=0;i<inmask3.size();i++){
	for(unsigned int j=0;j<inlrmask3.size();j++){
	  if(!opts.pathdist.value())
	    m_ConMat3->AddTo(inmask3[i].first+1,inlrmask3[j].first+1,1);
	  else{
	    float val = fabs(inmask3[i].second-inlrmask3[j].second);
	    m_ConMat3->AddTo(inmask3[i].first+1,inlrmask3[j].first+1,val);
	  }
	}
      }
    }
  }  

  void Counter::reset_beenhere3(){
    m_stline.clear_inmask3();
    m_stline.clear_inlrmask3();
  }  
  
  void Counter::reset_beenhere2(){
    ColumnVector xyz_lr(3);
    for(unsigned int i=0;i<m_path.size();i++){
      xyz_lr=vox_to_vox(m_path[i],m_seedsdim,m_lrdim,m_I);
      m_beenhere2((int)round((float)xyz_lr(1)),
		  (int)round((float)xyz_lr(2)),
		  (int)round((float)xyz_lr(3)))=0;
    }    
  }

  void Counter::save_total(const int& keeptotal){
    // save total number of particles that made it through the streamlining
    ColumnVector keeptotvec(1);
    keeptotvec(1)=keeptotal;
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));

  }
  void Counter::save_total(const vector<int>& keeptotal){
    // save total number of particles that made it through the streamlining
    ColumnVector keeptotvec(keeptotal.size());
    for (int i=1;i<=(int)keeptotal.size();i++)
      keeptotvec(i)=keeptotal[i-1];
    write_ascii_matrix(keeptotvec,logger.appendDir("waytotal"));

  }

  void Counter::save(){
    if(opts.simpleout.value() && !opts.simple.value()){
      save_pathdist();
    }
    if(opts.save_paths.value())
      save_paths();
    if(opts.s2tout.value()){
      save_seedcounts();
    }
    if(opts.matrix1out.value()){
      save_matrix1();
    }
    if(opts.matrix2out.value()){
      save_matrix2();
    }
    if(opts.matrix3out.value()){
      save_matrix3();
    }
  }
  
  void Counter::save_pathdist(){  
    m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
    save_volume(m_prob,logger.appendDir(opts.outfile.value()));
    if(opts.pathfile.set()){
      m_prob_alt.save_rois(logger.appendDir(opts.outfile.value())+"_alt");
      //m_beenhere_alt.save_rois(logger.appendDir(opts.outfile.value())+"_beenhere");
    }
    if(opts.opathdir.value()){
      for(int z=0;z<m_prob.zsize();z++){
	for(int y=0;y<m_prob.ysize();y++){
	  for(int x=0;x<m_prob.xsize();x++){
	    if(m_prob(x,y,z)==0)continue;
	    double norm=0;
	    for(int t=0;t<3;t++)
	      norm+=m_localdir(x,y,z,t)*m_localdir(x,y,z,t);
	    norm=sqrt(norm);
	    if(norm==0)continue;
	    for(int t=0;t<3;t++)
	      m_localdir(x,y,z,t) /= norm;
	  }
	}
      }
      save_volume(m_prob,logger.appendDir(opts.outfile.value()));
      m_localdir.setDisplayMaximumMinimum(1,-1);
      save_volume4D(m_localdir,logger.appendDir(opts.outfile.value()+"_localdir"));
    }
  }
  
  void Counter::save_pathdist(string add){  //for simple mode
    string thisout=opts.outfile.value();
    make_basename(thisout);
    thisout+=add;
    m_prob.setDisplayMaximumMinimum(m_prob.max(),m_prob.min());
    save_volume(m_prob,logger.appendDir(thisout));
    if(opts.pathfile.set()){
      m_prob_alt.save_rois(logger.appendDir(thisout)+"_alt");
    }
  }

  void Counter::save_seedcounts(){
    for(int m=0;m<m_targetmasks.nRois();m++){
      string tmpname=m_targetmasks.get_name(m);
      int pos=tmpname.find("/",0);
      int lastpos=pos;
      while(pos>=0){
	lastpos=pos;
	pos=tmpname.find("/",pos);
	// replace / with _
	tmpname[pos]='_';
      }      
      //only take things after the last pos
      tmpname=tmpname.substr(lastpos+1,tmpname.length()-lastpos-1);

      if(m_s2t_count.nRois()>1){
	for(int i=0;i<m_s2t_count.nRois();i++)
	  m_s2t_count.save_map(i,m,logger.appendDir("seeds_"+num2str(i)+"_to_"+tmpname));
      }
      else{// keep this nomenclature for backward compatibility
	m_s2t_count.save_map(0,m,logger.appendDir("seeds_to_"+tmpname));
      }	
    }

    if(opts.s2tastext.value()){
      write_ascii_matrix(m_s2tastext,logger.appendDir("matrix_seeds_to_all_targets"));
    }

    
  }
    
  // the following is a helper function for save_matrix*
  //  to convert between nifti coords (external) and newimage coord (internal)
  void applycoordchange(volume<int>& coordvol, const Matrix& old2new_mat)
  {
    for (int n=0; n<=coordvol.maxx(); n++) {
      ColumnVector v(4);
      v << coordvol(n,0,0) << coordvol(n,1,0) << coordvol(n,2,0) << 1.0;
      v = old2new_mat * v;
      coordvol(n,0,0) = MISCMATHS::round(v(1));
      coordvol(n,1,0) = MISCMATHS::round(v(2));
      coordvol(n,2,0) = MISCMATHS::round(v(3));
    }
  }
  void applycoordchange(Matrix& coordvol, const Matrix& old2new_mat)
  {
    for (int n=1; n<=coordvol.Nrows(); n++) {
      ColumnVector v(4);
      v << coordvol(n,1) << coordvol(n,2) << coordvol(n,3) << 1.0;
      v = old2new_mat * v;
      coordvol(n,1) = v(1);
      coordvol(n,2) = v(2);
      coordvol(n,3) = v(3);
    }
  }

  void Counter::save_matrix1(){
    m_ConMat1->Print(logger.appendDir("fdt_matrix1.dot"));
  }

  void Counter::save_matrix2(){
    m_ConMat2->Print(logger.appendDir("fdt_matrix2.dot")); 
  }


  void Counter::save_matrix3(){
    m_ConMat3->Print(logger.appendDir("fdt_matrix3.dot"));
  }

  int Seedmanager::run(const float& x,const float& y,const float& z,
		       bool onewayonly, int fibst,bool sampvox){
    int seed_ind=-1;
    // run without forcing the tracking direction
    return run(x,y,z,onewayonly,fibst,sampvox,seed_ind);
  }

  // this function returns the total number of pathways that survived a streamlining 
  int Seedmanager::run(const float& x,const float& y,const float& z,
		       bool onewayonly, int fibst,bool sampvox,const int& seed_ind){
    //Tracer_Plus tr("Seedmanager::run");

    //onewayonly for mesh things..
    if(opts.fibst.set()){
      fibst=opts.fibst.value()-1;      
    }
    else{
      if(fibst == -1){
	fibst=0;
      }
      //TB moved randfib option inside tractvols.h 28/10/2009
      // This means that we have access to fsamples when figuring out fibst
      // so we can choose to seed in proportion to f in that voxel. 
    }
  
  
    int nlines=0;
    for(int p=0;p<opts.nparticles.value();p++){
      if(opts.randfib.value()>2){ 
	//This bit of code just does random sampling from all fibre populations - even those whose f value is less than f-thresh. 
	//3 other possibilities - randfib==0 -> use fibst (default first fibre but can be set)
	// randfib==1 - random sampling of fibres bigger than fthresh
	// randfib==2 random sampling of fibres bigger than fthresh in proporthion to their f-values. 
	float tmp=rand()/float(RAND_MAX) * float(m_counter.get_stline().nfibres()-1);
	fibst = (int)round(tmp);
      }
    
      // random sampling within a seed voxel
      float newx=x,newy=y,newz=z;    
      if(sampvox){
	newx+=(float)rand()/float(RAND_MAX)-0.5;
	newy+=(float)rand()/float(RAND_MAX)-0.5;
	newz+=(float)rand()/float(RAND_MAX)-0.5;
      }
    
      if(opts.verbose.value()>1)
	logger.setLogFile("particle"+num2str(p));
    
      m_counter.get_stline().reset(); //This now includes a vols.reset() in order to get fibst right. 

      bool forwardflag=false,backwardflag=false;    
      int  rejflag1=1,rejflag2=1; // 0:accept, 1:reject, 2:wait
      m_counter.clear_path();
    
      // track in one direction
      if(!onewayonly || opts.matrix3out.value()){//always go both ways in matrix3 mode
	rejflag1 = m_counter.get_stline().streamline(newx,newy,newz,m_seeddims,fibst,seed_ind);

	if(rejflag1==0 || rejflag1==2){ 
	  forwardflag=true;
	  m_counter.append_path();
	  if(opts.save_paths.value())
	    m_counter.add_path();
	}
	m_counter.get_stline().reverse();
      }
      

      // track in the other direction
      // if rotdir!=0 then this tracks again in the same direction (doubles the number of samples...)
      rejflag2=m_counter.get_stline().streamline(newx,newy,newz,m_seeddims,fibst,seed_ind);

      if(rejflag2==0){	
	backwardflag=true;
      }
      if(rejflag2>0){
	backwardflag=false;
	if(rejflag1>0)
	  forwardflag=false;
      }
    
      if(!forwardflag)
	m_counter.clear_path();
      if(backwardflag){
	m_counter.append_path();   	
	if(opts.save_paths.value())
	  m_counter.add_path();
      }
      if(forwardflag || backwardflag){
	nlines++; 
	m_counter.count_streamline();
      
	if(opts.matrix3out.value()||opts.matrix1out.value()){
	  float pathlength;
	  if( !forwardflag || !backwardflag)
	    pathlength = m_counter.calc_pathlength(0);
	  else // if the paths are appended, one step has zero length
	    pathlength = m_counter.calc_pathlength(1);
	
	  if(opts.matrix3out.value() && pathlength>=opts.distthresh3.value())
	    m_counter.update_matrix3();
	  
	  if(opts.matrix1out.value() && pathlength>=opts.distthresh1.value())
	    m_counter.update_matrix1();
	}
      }

      m_counter.clear_streamline(); 
    }
  
    m_counter.count_seed();
  
  
    return nlines;
    
  }

}

  
  

