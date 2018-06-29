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

#include <UCBsplineSurface.h>
#include <UCBsplines.h>

#ifdef WIN32
#define WIN32ORSGI
#endif

#ifdef SGI
#define WIN32ORSGI
#endif

#ifdef SGI6
#define WIN32ORSGI
#endif

using namespace UCBspl;

// NOTE: namespace UCBspl:

SplineSurface::SplineSurface(std::shared_ptr<GenMatrixType> PHI,
                             double umin, double vmin, 
                             double umax, double vmax) {
  PHI_  = PHI;
  umin_ = umin;
  vmin_ = vmin;
  umax_ = umax;
  vmax_ = vmax;
}


SplineSurface::SplineSurface(const SplineSurface& surf) {
  PHI_  = surf.PHI_;
  umin_ = surf.umin_;
  vmin_ = surf.vmin_;
  umax_ = surf.umax_;
  vmax_ = surf.vmax_;
}

void SplineSurface::init(std::shared_ptr<GenMatrixType> PHI,
                         double umin, double vmin, 
                         double umax, double vmax) {
  PHI_  = PHI;
  umin_ = umin;
  vmin_ = vmin;
  umax_ = umax;
  vmax_ = vmax;
}


double SplineSurface::f(double u, double v) const { // original domain
  
  // NOTE for maintanance : This is the same algorithm as MBA::f_pure(u,v) apart from
  // adding base surface

  // Evaluation without base surface

  // Map to the half open domain Omega = [0,m) x [0,n)
  int m_ = PHI_->noX()-3;
  int n_ = PHI_->noY()-3;
  double uc = (u - umin_)/(umax_-umin_) * (double)m_;
  double vc = (v - vmin_)/(vmax_-vmin_) * (double)n_;
  
  int i, j;
  double s, t;
  UCBspl::ijst(m_, n_, uc, vc, i, j, s, t);

  double Bks[4];
  double Blt[4];

  Bks[0] = UCBspl::B_0(s); Blt[0] = UCBspl::B_0(t);
  Bks[1] = UCBspl::B_1(s); Blt[1] = UCBspl::B_1(t);
  Bks[2] = UCBspl::B_2(s); Blt[2] = UCBspl::B_2(t);
  Bks[3] = UCBspl::B_3(s); Blt[3] = UCBspl::B_3(t);

    // following commented-out part replaced by the four above lines for 
    // optimality reasons (Odd Andersen 15. dec. 2003)

//   for (int kl = 0; kl < 4; kl++) {
//     Bks[kl] = UCBspl::B(kl,s);
//     Blt[kl] = UCBspl::B(kl,t);
//   }

  double val = 0.0;
  for (int k = 0; k <= 3; k++)
    for (int l = 0; l <= 3; l++)
      //val += PHI_(i+k,j+l)*B(k,s)*B(l,t);
      val += (*PHI_)(i+k,j+l)*Bks[k]*Blt[l];
        
  return val;
}


double SplineSurface::f(int ii, int jj) const {
  
  // Size of matrix in UnifBsplines for pre-evaluated bases,
  // and index of lower left coefficient
#ifdef  UNIFORM_CUBIC_C1_SPLINES
  int i = 2*ii - 1;
  int j = 2*jj - 1;
  int m_size = 2;
#else
  int i = ii - 1;
  int j = jj - 1;
  int m_size = 3;
#endif
      
  double val = 0.0;
  // Only over 3x3 coefficients in the C2 case and 2x2 in the C1 case
  for (int k = 0; k < m_size; k++) {
    for (int l = 0; l < m_size; l++) {
      val += (*PHI_)(i+k,j+l)*UCBspl::tensor_BB[k][l];
    }
  }

  return val;
}


void SplineSurface::normalVector(int ii, int jj, double& gx, double& gy, double& gz) const { // original domain

#ifdef  UNIFORM_CUBIC_C1_SPLINES
  int i = 2*ii - 1;
  int j = 2*jj - 1;
  int m_size = 2;
  int incr = 1;
#else
  int i = ii - 1;
  int j = jj - 1;
  int m_size = 3; // size of matrices with pre-evaluate tensors
  int incr = 2;
#endif

  int k,l;
  double val1 = 0.0;
  // note for C2: k += 2 to skip second row of tensor
  for (k = 0; k < m_size; k+=incr) {
    for (l = 0; l < m_size; l++) {
      val1 += (*PHI_)(i+k,j+l) * UCBspl::tensor_dBB[k][l];
    }
  }
  double val2 = 0.0;
  // note for C2: l += 2 to skip second column of tensor
  for (k = 0; k < m_size; k++) {
    for (l = 0; l < m_size; l+=incr) {
      val2 += (*PHI_)(i+k,j+l) * UCBspl::tensor_BdB[k][l];
    }
  }

  // find cross product and normalize
  // (-df/du, -df/dv, 1)
  int m_ = PHI_->noX()-3;
  int n_ = PHI_->noY()-3;
  val1 *= ((double)m_/(umax_-umin_));
  val2 *= ((double)n_/(vmax_-vmin_));

  double len = sqrt(val1*val1 + val2*val2 + 1.0);

  gx = -(val1/len);
  gy = -(val2/len);
  gz = 1.0/len;
}

// Added by Vegard Stenerud July 2001, but not properly tested
// -----------------------------------------------------------
void SplineSurface::derivatives(double u, double v, double& dx, double& dy) const { // original domain
  
  // Map to the half open domain Omega = [0,m) x [0,n)
  int m_ = PHI_->noX()-3;
  int n_ = PHI_->noY()-3;
  double uc = (u - umin_)/(umax_-umin_) * (double)m_;
  double vc = (v - vmin_)/(vmax_-vmin_) * (double)n_;
  
  int i, j;
  double s, t;
  UCBspl::ijst(m_, n_, uc, vc, i, j, s, t);
  
  double  Bks[4];
  double  Blt[4];
  double dBks[4];
  double dBlt[4];
  
    Bks[0] = UCBspl::B_0(s); Blt[0] = UCBspl::B_0(t);
    Bks[1] = UCBspl::B_1(s); Blt[1] = UCBspl::B_1(t);
    Bks[2] = UCBspl::B_2(s); Blt[2] = UCBspl::B_2(t);
    Bks[3] = UCBspl::B_3(s); Blt[3] = UCBspl::B_3(t);


    dBks[0] = UCBspl::dB_0(s); dBlt[0] = UCBspl::dB_0(t);
    dBks[1] = UCBspl::dB_1(s); dBlt[1] = UCBspl::dB_1(t);
    dBks[2] = UCBspl::dB_2(s); dBlt[2] = UCBspl::dB_2(t);
    dBks[3] = UCBspl::dB_3(s); dBlt[3] = UCBspl::dB_3(t);

    // following commented-out part replaced by the eight above lines for 
    // optimality reasons (Odd Andersen 15. dec. 2003)

//   for (int kl = 0; kl < 4; kl++) {
//     Bks[kl] =  UCBspl::B(kl,s);
//     Blt[kl] =  UCBspl::B(kl,t);
//     dBks[kl] = UCBspl::dB(kl,s);
//     dBlt[kl] = UCBspl::dB(kl,t);
//   }
  
  double val1 = 0.0;
  double val2 = 0.0;
  
  for (int k = 0; k <= 3; k++) {
    for (int l = 0; l <= 3; l++) {
      val1 += (*PHI_)(i+k,j+l) * dBks[k]*  Blt[l];
      val2 += (*PHI_)(i+k,j+l) *  Bks[k] * dBlt[l];
    }
  }
  
  // find cross product and normalize
  // (-df/du, -df/dv, 1)
  val1 *= ((double)m_/(umax_-umin_));
  val2 *= ((double)n_/(vmax_-vmin_));
  
  dx = val1;
  dy = val2;
  
}

// Added by Vegard Stenerud July 2001, but not properly tested
// -----------------------------------------------------------
void SplineSurface::secondDerivatives(double u, double v, double& ddx, double& ddy, double& dxdy) const { // original domain
  
  // Map to the half open domain Omega = [0,m) x [0,n)
  int m_ = PHI_->noX()-3;
  int n_ = PHI_->noY()-3;
  double uc = (u - umin_)/(umax_-umin_) * (double)m_;
  double vc = (v - vmin_)/(vmax_-vmin_) * (double)n_;
  
  int i, j;
  double s, t;
  UCBspl::ijst(m_, n_, uc, vc, i, j, s, t);
  
  double   Bks[4];
  double   Blt[4];
  double  dBks[4];
  double  dBlt[4];
  double ddBks[4];
  double ddBlt[4];
  
  Bks[0] = UCBspl::B_0(s); Blt[0] = UCBspl::B_0(t);
  Bks[1] = UCBspl::B_1(s); Blt[1] = UCBspl::B_1(t);
  Bks[2] = UCBspl::B_2(s); Blt[2] = UCBspl::B_2(t);
  Bks[3] = UCBspl::B_3(s); Blt[3] = UCBspl::B_3(t);

  dBks[0] = UCBspl::dB_0(s); dBlt[0] = UCBspl::dB_0(t);
  dBks[1] = UCBspl::dB_1(s); dBlt[1] = UCBspl::dB_1(t);
  dBks[2] = UCBspl::dB_2(s); dBlt[2] = UCBspl::dB_2(t);
  dBks[3] = UCBspl::dB_3(s); dBlt[3] = UCBspl::dB_3(t);

  ddBks[0] = UCBspl::ddB_0(s); ddBlt[0] = UCBspl::ddB_0(t);
  ddBks[1] = UCBspl::ddB_1(s); ddBlt[1] = UCBspl::ddB_1(t);
  ddBks[2] = UCBspl::ddB_2(s); ddBlt[2] = UCBspl::ddB_2(t);
  ddBks[3] = UCBspl::ddB_3(s); ddBlt[3] = UCBspl::ddB_3(t);

  // following commented-out part replaced by the twelve above lines for 
  // optimality reasons (Odd Andersen 15. dec. 2003)


//   for (int kl = 0; kl < 4; kl++) {
//     Bks[kl] =   UCBspl::B(kl,s);
//     Blt[kl] =   UCBspl::B(kl,t);
//     dBks[kl] =  UCBspl::dB(kl,s);
//     dBlt[kl] =  UCBspl::dB(kl,t);
//     ddBks[kl] = UCBspl::ddB(kl,s);
//     ddBlt[kl] = UCBspl::ddB(kl,t);
//   }
  
  double val1 = 0.0;
  double val2 = 0.0;
  double val3 = 0.0;
  
  for (int k = 0; k <= 3; k++) {
    for (int l = 0; l <= 3; l++) {
      val1 += (*PHI_)(i+k,j+l) * ddBks[k] *   Blt[l];
      val2 += (*PHI_)(i+k,j+l) *   Bks[k] * ddBlt[l];
      val3 += (*PHI_)(i+k,j+l) *  dBks[k] *  dBlt[l];
    }
  }
  
  // Har ikke sett p� hvordan det blir med skaleringa her.
  
  // find cross product and normalize
  // (-df/du, -df/dv, 1)
  val1 *= ((double)m_/(umax_-umin_));
  val2 *= ((double)n_/(vmax_-vmin_));
  //val3 *= ((double)n_/(vmax_-vmin_));
  
  ddx = val1;
  ddy = val2;
  dxdy = val3;
  
}

// Added by Vegard Stenerud July 2001, but not properly tested
// -----------------------------------------------------------
void SplineSurface::curvatures(double u, double v, double& profC, double& planC) const { // original domain
  
  double ddx;
  double ddy;
  double dxdy;
  
  secondDerivatives(u,v,ddx,ddy,dxdy);
  
  double gx;
  double gy;
  double gz;
  
  normalVector(u,v,gx,gy,gz);
  
  // gx^2 + gy^2 + gz^2 = 1
  double len = sqrt(1-gz*gz);
  
  double u1 = gx/len;
  double u2 = gy/len;
  
  double v1 = -gy/len;
  double v2 = gx/len; 
  
  //Henviser til notat om retnings andre deriverte og
  //Earth surface processes and landforms, vol 12, 47-56 (1987) (Fikk av J. R. Sulebakk)
  
  double dirDDProfC = u1*u1*ddx + 2*u1*u2*dxdy + u2*u2*ddy;
  double dirDDPlanC = v1*v1*ddx + 2*v1*v2*dxdy + v2*v2*ddy;
  
  
  profC = dirDDProfC/((1+(dirDDProfC)*(dirDDProfC))*sqrt(1+(dirDDProfC)*(dirDDProfC)));
  planC = dirDDPlanC/((1+(dirDDPlanC)*(dirDDPlanC))*sqrt(1+(dirDDPlanC)*(dirDDPlanC)));
  
}
// -------------------- Vegard -------------------

void SplineSurface::eval(int i, int j, double& z, double& gx, double& gy, double& gz) const {
  // Just to complete the interface also for eval
  // No optimization compared to calling f(i,j) and grad(i,j,...) separately

  z = f(i,j);
  normalVector(i, j, gx, gy, gz);
}

void SplineSurface::normalVector(double u, double v, double& gx, double& gy, double& gz) const { // original domain
  
  // Map to the half open domain Omega = [0,m) x [0,n)
  int m_ = PHI_->noX()-3;
  int n_ = PHI_->noY()-3;
  double uc = (u - umin_)/(umax_-umin_) * (double)m_;
  double vc = (v - vmin_)/(vmax_-vmin_) * (double)n_;
  
  int i, j;
  double s, t;
  UCBspl::ijst(m_, n_, uc, vc, i, j, s, t);
  
  double  Bks[4];
  double  Blt[4];
  double dBks[4];
  double dBlt[4];

  Bks[0] = UCBspl::B_0(s); Blt[0] = UCBspl::B_0(t);
  Bks[1] = UCBspl::B_1(s); Blt[1] = UCBspl::B_1(t);
  Bks[2] = UCBspl::B_2(s); Blt[2] = UCBspl::B_2(t);
  Bks[3] = UCBspl::B_3(s); Blt[3] = UCBspl::B_3(t);

  dBks[0] = UCBspl::dB_0(s); dBlt[0] = UCBspl::dB_0(t);
  dBks[1] = UCBspl::dB_1(s); dBlt[1] = UCBspl::dB_1(t);
  dBks[2] = UCBspl::dB_2(s); dBlt[2] = UCBspl::dB_2(t);
  dBks[3] = UCBspl::dB_3(s); dBlt[3] = UCBspl::dB_3(t);


  // following commented-out part replaced by the eight above lines for 
  // optimality reasons (Odd Andersen 15. dec. 2003)

//   for (int kl = 0; kl < 4; kl++) {
//      Bks[kl] = UCBspl::B(kl,s);
//      Blt[kl] = UCBspl::B(kl,t);
//     dBks[kl] = UCBspl::dB(kl,s);
//     dBlt[kl] = UCBspl::dB(kl,t);
//   }
  
  double val1 = 0.0;
  double val2 = 0.0;
  for (int k = 0; k <= 3; k++) {
    for (int l = 0; l <= 3; l++) {
      val1 += (*PHI_)(i+k,j+l) * dBks[k]*  Blt[l];
      val2 += (*PHI_)(i+k,j+l) *  Bks[k] * dBlt[l];
    }
  }
  
  // find cross product and normalize
  // (-df/du, -df/dv, 1)
  val1 *= ((double)m_/(umax_-umin_));
  val2 *= ((double)n_/(vmax_-vmin_));

  double len = sqrt(val1*val1 + val2*val2 + 1.0);

  gx = -(val1/len);
  gy = -(val2/len);
  gz = 1.0/len;
}

void SplineSurface::eval(double u, double v, double& z, double& gx, double& gy, double& gz) const {

  // Map to the half open domain Omega = [0,m) x [0,n)
  int m_ = PHI_->noX()-3;
  int n_ = PHI_->noY()-3;
  double uc = (u - umin_)/(umax_-umin_) * (double)m_;
  double vc = (v - vmin_)/(vmax_-vmin_) * (double)n_;
  
  int i, j;
  double s, t;
  UCBspl::ijst(m_, n_, uc, vc, i, j, s, t);
  
  double Bks[4];
  double Blt[4];
  double dBks[4];
  double dBlt[4];

  Bks[0] = UCBspl::B_0(s); Blt[0] = UCBspl::B_0(t);
  Bks[1] = UCBspl::B_1(s); Blt[1] = UCBspl::B_1(t);
  Bks[2] = UCBspl::B_2(s); Blt[2] = UCBspl::B_2(t);
  Bks[3] = UCBspl::B_3(s); Blt[3] = UCBspl::B_3(t);

  dBks[0] = UCBspl::dB_0(s); dBlt[0] = UCBspl::dB_0(t);
  dBks[1] = UCBspl::dB_1(s); dBlt[1] = UCBspl::dB_1(t);
  dBks[2] = UCBspl::dB_2(s); dBlt[2] = UCBspl::dB_2(t);
  dBks[3] = UCBspl::dB_3(s); dBlt[3] = UCBspl::dB_3(t);

  // following commented-out part replaced by the eight above lines for 
  // optimality reasons (Odd Andersen 15. dec. 2003)

//   for (int kl = 0; kl < 4; kl++) {
//     Bks[kl] = UCBspl::B(kl,s);
//     Blt[kl] = UCBspl::B(kl,t);
//     dBks[kl] = UCBspl::dB(kl,s);
//     dBlt[kl] = UCBspl::dB(kl,t);
//   }

  double val1 = 0.0;
  double val2 = 0.0;

  z = 0.0;
  for (int k = 0; k <= 3; k++) {
    for (int l = 0; l <= 3; l++) {
      z    += (*PHI_)(i+k,j+l) *  Bks[k] *  Blt[l];
      val1 += (*PHI_)(i+k,j+l) * dBks[k] *  Blt[l];
      val2 += (*PHI_)(i+k,j+l) *  Bks[k] * dBlt[l];
    }
  }

  // find cross product and normalize
  // (-df/du, -df/dv, 1)
  val1 *= ((double)m_/(umax_-umin_));
  val2 *= ((double)n_/(vmax_-vmin_));
  
  double len = sqrt(val1*val1 + val2*val2 + 1.0);

  gx = -(val1/len);
  gy = -(val2/len);
  gz = 1.0/len;
}

// ---------------------------------------------------------
void SplineSurface::refineCoeffs() {

  GenMatrix<UCBspl_real>* PHIrefined = new GenMatrix<UCBspl_real>();
  UCBspl::refineCoeffsC2(*PHI_, *PHIrefined);	
  PHI_.reset(PHIrefined);
}

// ---------------------------------------------------------
bool SplineSurface::restrictCoeffs() {

  GenMatrix<UCBspl_real>* PHIrestricted = new GenMatrix<UCBspl_real>();
  bool status = UCBspl::restrictCoeffsC2(*PHI_, *PHIrestricted);	
  if (!status)
    return status;
  PHI_.reset(PHIrestricted);

  return true;
}    
