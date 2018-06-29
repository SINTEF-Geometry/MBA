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

#ifndef _MBADATA_H_
#define _MBADATA_H_

#include <MBAtypedef.h>
#include <vector>
#include <memory>

//===========================================================================
/** \brief Holds the scattered data
 * 
 * MBAdata - Holds scattered data for Multilevel B-spline approximation
 * (and interpolation) for functional surfaces.
 * Note that this is a (mathematical) lower level interface; thus there is no check
 * for the validity of given scattered data or of arguments passed to
 * member functions.
 * See further documentation in MBA and MBAadaptive
 * 
 * \author �yvind Hjelle <Oyvind.Hjelle@math.sintef.no>
 * \see MBA, MBAadaptive, MBAdataPar
 */
//===========================================================================
class MBAdata {
  friend class MBA;
  friend class MBAadaptive;

  double umin_, vmin_, umax_, vmax_; // possibly user defined (expanded)
    
    double urange_inv_;
    double vrange_inv_;

  MBAbaseType  baseType_;  
  double offset_;
  std::shared_ptr<dVec> U_;
  std::shared_ptr<dVec> V_;
  std::shared_ptr<dVec> Zorig_;
    std::vector<double> Z_;

  /// Read scattered data from file
    void readScatteredData(const char filename[]);

  void buildOffset();
  void buildBaseSurface();

  /*  Expand the rectangular domain beyond the domain of the given 
   *  scattered data; see doc in MBA.
   *  This must be done after data initialization.
   *  Note that this can only be used to expand the domain 
   *  (and not shrink it in any direction). 
   */
  //void expandDomain(double Dumin, double Dvmin, double Dumax, double Dvmax) {
    //umin_ -= Dumin; vmin_ -= Dvmin; umax_ += Dumax; vmax_ += Dvmax;}

  //  Set the rectangular domain. See also above.
  void setDomain(double umin, double vmin, double umax, double vmax)
  {
      umin_ = umin; 
      vmin_ = vmin; 
      umax_ = umax; 
      vmax_ = vmax;
      urange_inv_ = double(1) / (umax_ - umin_);
      vrange_inv_ = double(1) / (vmax_ - vmin_);
  }

  // Clear all allocated memory
  void clear() {U_->clear(); V_->clear(); Z_.clear(); Zorig_->clear();}

  // Non-const for class MBA
  std::vector<double>& Z() {return Z_;};

  // From data limits
  void initDefaultDomain();

public:
	MBAdata();
	~MBAdata() {}

  /// Initialize with scattered data
  void init(std::shared_ptr<dVec> U, std::shared_ptr<dVec> V, std::shared_ptr<dVec> Z);
  
  /// min u-value of the actual data domain
  const double& umin() const {return umin_;}
  /// min v-value of the actual data domain
  const double& vmin() const {return vmin_;}
  /// max u-value of the actual data domain
  const double& umax() const {return umax_;}
  /// max v-value of the actual data domain
  const double& vmax() const {return vmax_;}

  const double& rangeUInv() const {return urange_inv_;}
  const double& rangeVInv() const {return vrange_inv_;}

  // Access to scattered data
  const std::shared_ptr<dVec>& U() const {return U_;};
  const std::shared_ptr<dVec>& V() const {return V_;};
  const std::shared_ptr<dVec>& Zorig() const {return Zorig_;};
  const std::vector<double>& Z() const {return Z_;}; // also non-const private members
  double U(int i) const {return (*U_)[i];};
  double V(int i) const {return (*V_)[i];};

  /// Number of scattered data points
    int size() const {return U_->size();}


  // Evaluator.
  // Assumes that the base surface is a constant
  inline double f() const {return offset_;}
};

#endif //_MBADATA_H_
