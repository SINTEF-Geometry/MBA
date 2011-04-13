//================================================================
//
// Created: September 18. 2000
//                                                                           
// Author: Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description:
//                                                                           
//================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//================================================================

#ifndef _MBADATA_H_
#define _MBADATA_H_

#include <MBAtypedef.h>
#include <vector>
#include <boost/shared_ptr.hpp>

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
 * \author Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
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
  boost::shared_ptr<dVec> U_;
  boost::shared_ptr<dVec> V_;
  boost::shared_ptr<dVec> Zorig_;
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
  void init(boost::shared_ptr<dVec> U, boost::shared_ptr<dVec> V, boost::shared_ptr<dVec> Z);
  
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
  const boost::shared_ptr<dVec>& U() const {return U_;};
  const boost::shared_ptr<dVec>& V() const {return V_;};
  const boost::shared_ptr<dVec>& Zorig() const {return Zorig_;};
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
