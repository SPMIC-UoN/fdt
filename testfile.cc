/*  Copyright (C) 2004 University of Oxford  */

/*  CCOPYRIGHT  */

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;
int read_avg_file (vector<vector<int> >& avgs,const string fname){
  avgs.clear();
  ifstream avg_file(fname.c_str());
  if(!avg_file){return -1;}
  avno=0;
  int start,length;
  string myline;
  else{
      avg_file.getline(myline);
      istrstream mylinestr(myline);
      mylinestr>>start;
      mylinestr>>length;
      cerr<<start<<" "<<length<<endl;
      mylinestr>>start;
      mylinestr>>length;
      cerr<<start<<" "<<length<<endl;

    }
  }
  
  return 0;
}
 

int main(int argc, char** argv){
  vector<vector<int> > poo;
  string af="avg_file";
  read_avg_file(poo,af);
}
