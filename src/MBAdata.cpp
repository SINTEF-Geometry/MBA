/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of MBA.
 *
 * MBA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * MBA is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with MBA. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using MBA.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the MBA library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */


#include "checkWIN32andSGI.h"
#include <MBAdata.h>
#include <algorithm>

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

void MBAdata::init(std::shared_ptr<dVec> U, std::shared_ptr<dVec> V, std::shared_ptr<dVec> Z) {

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

void MBAdata::readScatteredData(const char filename[]) {
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
  if (V_.get() == NULL)
      V_.reset(new std::vector<double>);
  V_->resize(no);
  Z_.resize(no);
  if (Zorig_.get() == NULL)
      Zorig_.reset(new std::vector<double>);
  Zorig_->resize(no);
	
  for (int i = 0; i < no; i++) {
    ifile >> (*U_)[i] >> (*V_)[i] >> (*Zorig_)[i];
    Z_[i] = (*Zorig_)[i];
  }
}
