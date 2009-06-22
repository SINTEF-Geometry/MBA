//================================================================
//
// Created: September 18. 2000
//                                                                           
// Author: �yvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description:
//                                                                           
//================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//================================================================

#include "checkWIN32andSGI.h"
#include <MBAdata.h>

#ifdef WIN32
#include <minmax.h>
#endif

#ifdef SGI
#include <math.h>
#else
#include <cmath>
#endif
#include <iostream>
#include <fstream>
using namespace std;


MBAdata::MBAdata() {
  baseType_ = MBA_CONSTLS;
  offset_ = 0.0;
  umin_=vmin_=umax_=vmax_=MBA_UNDEFREAL;
  urange_inv_ = vrange_inv_ = MBA_UNDEFREAL;

}

void MBAdata::init(boost::shared_ptr<dVec> U, boost::shared_ptr<dVec> V, boost::shared_ptr<dVec> Z) {

#ifdef MBA_DEBUG  
  if (U->size() != V->size() || Z->size() != U->size()) {
    cout << "ERROR, U, V and Z must have equal size" << endl;
    return;
  }
#endif

  U_ = U;
  V_ = V;
  Zorig_ = Z;
  Z_.assign(Z->begin(), Z->end()); // Copy

  //initDefaultDomain();
}

void MBAdata::initDefaultDomain() {

  if (U_->size() == 0)
    return;

  umin_ = *std::min_element(U_->begin(),U_->end());
  vmin_ = *std::min_element(V_->begin(),V_->end());
  umax_ = *std::max_element(U_->begin(),U_->end());
  vmax_ = *std::max_element(V_->begin(),V_->end());

  urange_inv_ = double(1) / (umax_ - umin_);
  vrange_inv_ = double(1) / (vmax_ - vmin_);

}


void MBAdata::buildOffset() {
  int noPoints = Z_.size();
  for (int ip = 0; ip < noPoints; ip++)
    Z_[ip] = (*Zorig_)[ip] - offset_;  
}

static double average(const std::vector<double>& vec) {
  int no = vec.size();
  double sum = 0.0;
  for (int ip = 0; ip < no; ip++)
     sum += vec[ip];

  return sum/(double)no;
}


void MBAdata::buildBaseSurface() {
  if (baseType_ == MBA_CONSTLS) {
    offset_ = average(Z_);
#ifdef MBA_DEBUG
    cout << "averageZ = " << offset_ << endl;
#endif
    buildOffset();
  }
  else if (baseType_ == MBA_CONSTVAL) {
    buildOffset();
  }
}

void MBAdata::readScatteredData(char filename[]) {
#ifdef MBA_DEBUG
  cout << "Reading data from file..." << endl;
#endif
  
  ifstream ifile(filename);
  double u,v,z,zmin,zmax;
  int no=0;
  umin_ = vmin_ = zmin =  1.0e+20;
  umax_ = vmax_ = zmax = -1.0e+20;
	
  do {
    ifile >> u >> v >> z;
    if (ifile) {
      no++;
      umin_ = min(u, umin_);
      vmin_ = min(v, vmin_);
      umax_ = max(u, umax_);
      vmax_ = max(v, vmax_);

      zmin = min(z,zmin);
      zmax = max(z,zmax);
    }
  } while(ifile);
  

#ifdef MBA_DEBUG
  cout << "No. of points read = " << no << " min and max: " << endl;
  cout << "(" << umin_ << ", " << vmin_ << ", " << zmin << ") " << "(" << umax_ << ", " << vmax_ << ", " << zmax << ")" << endl;
#endif

  urange_inv_ = double(1) / (umax_ - umin_);
  vrange_inv_ = double(1) / (vmax_ - vmin_);

  ifile.clear();
  ifile.seekg(0);
	
  if (U_.get() == NULL)
    U_.reset(new std::vector<double>);
  U_->resize(no);
  V_->resize(no);
  Z_.resize(no);
  Zorig_->resize(no);
	
  for (int i = 0; i < no; i++) {
    ifile >> (*U_)[i] >> (*V_)[i] >> (*Zorig_)[i];
    Z_[i] = (*Zorig_)[i];
  }
}
