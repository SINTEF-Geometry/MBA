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

#ifndef _UNIFBSPLINES_H_
#define _UNIFBSPLINES_H_

#include <GenMatrix.h>

#include <cmath>
using namespace std;

//#define MBA_CUBIC_C1 1 

#ifdef  MBA_CUBIC_C1
// Basis functions
inline double B_0(double t) {double tmp = 1.0-t; return tmp*tmp*tmp/2.0;}
inline double B_1(double t) {return (2.5*t*t*t - 4.5*t*t + 1.5*t + 0.5);}
inline double B_2(double t) {return (-2.5*t*t*t + 3.0*t*t);}
inline double B_3(double t) {return t*t*t/2.0;}

// And their derivatives
inline double dB_0(double t) {return 1.5*(1.0-t)*(t-1.0);}
inline double dB_1(double t) {return (7.5*t*t - 9.0*t + 1.5);}
inline double dB_2(double t) {return (-7.5*t*t + 6.0*t);}
inline double dB_3(double t) {return 1.5*t*t;}

// PRE EVALUATTED tensor products of B-splines, in a grid index
// -------------------------------------------------------------

// Pre-evaluated tensors equivalent to B(k,s)*B(l,t)
static double tensor_BB[2][2] = {
  {1./4., 1./4.},
  {1./4., 1./4.}
};

// ???? only temp. for compilation
// equivalent to dB(k,s)*B(l,t).
static double tensor_dBB[2][2] = {
  {-3./4., -3./4.},
  { 3./4.,  3./4.}
};

// The transpose of the matrix above, for the y-components,
// equivalent to B(k,s)*dB(l,t)
static double tensor_BdB[2][2] = {
  {-3./4., 3./4.},
  {-3./4., 3./4.}
};

#else
// The uniform (cardinal), univariate cubic C2 basis functions
// 0 <= t < 1;
inline double B_0(double t) {double tmp = 1.0-t; return tmp*tmp*tmp/6.0;}

//static inline double B_1(double t) {return (3.0*t*t*t - 6.0*t*t+4.0)/6.0;}
static double div46 = 4.0/6.0;
inline double B_1(double t) {return (0.5*t*t*t - t*t + div46);}

//static inline double B_2(double t) {return (-3.0*t*t*t + 3.0*t*t + 3.0*t + 1.0)/6.0;}
static double div16 = 1.0/6.0;
inline double B_2(double t) {return (-0.5*t*t*t + 0.5*t*t + 0.5*t + div16);}

inline double B_3(double t) {return t*t*t/6.0;}

// FIRST DERIVATIVES to the uniform cubic C2 basis functions
inline double dB_0(double t) {return 0.5*(1.0-t)*(t-1.0);}
inline double dB_1(double t) {return (1.5*t*t - 2.0*t);}
inline double dB_2(double t) {return (-(1.5*t*t) + t + 0.5);}
inline double dB_3(double t) {return t*t/2.0;}

// Lagt til, Vegard
// SECOND DERIVATIVES to the uniform cubic C2 basis functions
inline double ddB_0(double t) {return (-t+1.0);}
inline double ddB_1(double t) {return (3.0*t-2.0);}
inline double ddB_2(double t) {return (-3.0*t+1.0);}
inline double ddB_3(double t) {return t;}

// PRE EVALUATTED tensor products of B-splines, in a grid index
// -------------------------------------------------------------

// Pre-evaluated tensors equivalent to B(k,s)*B(l,t)
static double tensor_BB[3][3] = {
  {1./36., 1./9., 1./36.},
  {1./9., 4./9., 1./9.},
  {1./36., 1./9., 1./36.}
};


// Pre-evaluated tensor products of B-splines, in a grid index
// Mixed derivative an positional tensor for the x-component.
// equivalent to dB(k,s)*B(l,t).
// (The second row with zeros is not used in evaluations.)
static double tensor_dBB[3][3] = {
  {-1./12., -1./3., -1./12.},
  {     0.,     0.,      0.},
  { 1./12.,  1./3.,  1./12.}
};

// The transpose of the matrix above, for the y-components,
// equivalent to B(k,s)*dB(l,t)
static double tensor_BdB[3][3] = {
  {-1./12., 0., 1./12.},
  { -1./3., 0., 1./3.},
  {-1./12., 0., 1./12.}
};

#endif


// Generic interface for evaluation to all univariate basis functions
inline double B(int k, double t) {
  if (k == 0)
    return B_0(t);
  else if (k == 1)
    return B_1(t);
  else if (k == 2)
    return B_2(t);
  else
    return B_3(t);
}


// Generic interface to first derivatives
inline double dB(int k, double t) {
  if (k == 0)
    return dB_0(t);
  else if (k == 1)
    return dB_1(t);
  else if (k == 2)
    return dB_2(t);
  else
    return dB_3(t);
}

// Lagt til, Vegard
// Generic interface to second derivatives
inline double ddB(int k, double t) {
  if (k == 0)
    return ddB_0(t);
  else if (k == 1)
    return ddB_1(t);
  else if (k == 2)
    return ddB_2(t);
  else
    return ddB_3(t);
}

// Evaluation of a bivariate tensor product basis function
inline double w(int k, int l, double s, double t) {
  return B(k,s)*B(l,t);
}

// s,t \in [0,1) (but special on gridlines m and n)
// i,j \in [-1, ???
inline void ijst(int m, int n, double uc, double vc, int& i, int& j, double& s, double& t) {
	//int i = std::min((int)uc - 1, m-2);
	//int j = std::min((int)vc - 1, n-2);

#ifdef  MBA_CUBIC_C1
	i = 2*((int)uc) - 1;
	j = 2*((int)vc) - 1;
#else
	i = (int)uc - 1;
	j = (int)vc - 1;
#endif

#ifdef WIN32
  s = uc - floor(uc);
	t = vc - floor(vc);
#else
  s = uc - std::floor(uc);
  t = vc - std::floor(vc);
#endif
	
	// adjust for x or y on gridlines m and n (since impl. has 0 <= x <= m and 0 <= y <= n 
#ifdef  MBA_CUBIC_C1
	if (i == 2*m-1) {
		i-=2;
		s = 1;
	}
	if (j == 2*n-1) {
		j-=2;
		t = 1;
	}
#else
	if (i == m-1) {
		i--;
		s = 1;
	}
	if (j == n-1) {
		j--;
		t = 1;
	}
#endif
}

inline void WKLandSum2(double s, double t, double w_kl[4][4], double& sum_w_ab2) {
  int k,l;
  sum_w_ab2 = 0.0; 
  for (k = 0; k <= 3; k++) {
    for (l = 0; l <=3; l++) {
      double tmp = w(k, l, s, t);
      w_kl[k][l] = tmp;
      sum_w_ab2 += (tmp*tmp); 
    }
  }
}

// for check, should give unity
// inline double sumWKL(double s, double t) {
//  int k,l;
//  double sum_w_ab = 0.0; 
//  for (k = 0; k <= 3; k++) {
//    for (l = 0; l <=3; l++) {
//      sum_w_ab += w(k, l, s, t);
//    }
//  }
//  return sum_w_ab;
//}

// -------------------------------------------------------------------------
// Refinement (The Oslo algorithm)
// Generic implementation for functional and parametric case
// There are 4 different configurations of vertices in the refined lattice:

// For the uniform bicubic C2 case

// 1) On an existing vertex
template <class Type>
inline Type phi_2i_2j(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/64.0*( mat(i-1,j-1) + mat(i-1,j+1) + mat(i+1,j-1) + mat(i+1,j+1)
    +6.0*(mat(i-1,j) + mat(i,j-1) + mat(i,j+1) + mat(i+1,j)) + 36.0*mat(i,j) );
}

// 2) On an existing vertical line
template <class Type>
inline Type phi_2i_2jPluss1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/16.0*( mat(i-1,j) + mat(i-1,j+1) + mat(i+1,j) + mat(i+1,j+1)
    +6.0*(mat(i,j) + mat(i,j+1)) );
}

// 3) On an existing horizontal line
template <class Type>
inline Type phi_2iPlus1_2j(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/16.0*( mat(i,j-1) + mat(i,j+1) + mat(i+1,j-1) + mat(i+1,j+1)
    +6.0*(mat(i,j) + mat(i+1,j)) );
}

// 4) a new crossing in both directions
template <class Type>
inline Type phi_2iPlus1_2jPlus1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/4.0*( mat(i,j) + mat(i,j+1) + mat(i+1,j) + mat(i+1,j+1) );
}


// The Oslo Algorithm For the uniform bicubic C1 case 
// There are 16 configurations:
// ----------------------------
//
//      |       :
//      |       :
//      |       :
//.....9|13...15:4...........
//     7|11    3:16
//      |       :
//      |       :
//_____5|2____12:14__________
//     1|6     8:10
//      |       :
//      |       :
//      |       :

// !!!! int i, int j denotes (2i) and (2j) in the coefficient matrix PHI !!!
// 1) SW of an existing grid point
template <class Type>
inline Type phi_4iMin1_4jMin1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/16.0*( 9.0*mat(i-1,j-1) + 3.0*(mat(i-1,j) + mat(i,j-1)) + mat(i,j) );
}

// 2) NE of an existing grid point
template <class Type>
inline Type phi_4i_4j(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/16.0*( mat(i-1,j-1) + 3.0*(mat(i-1,j) + mat(i,j-1)) + 9.0*mat(i,j) );
}

// 3) SW of a new grid point
template <class Type>
inline Type phi_4iPlus1_4jPlus1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/64.0*( mat(i-1,j-1) + 5.0*(mat(i-1,j) + mat(i,j-1)) + 25.0*mat(i,j) + 2.0*(mat(i-1,j+1) + mat(i+1,j-1)) + 10.0*(mat(i+1,j) + mat(i,j+1)) + 4.0*mat(i+1,j+1) );
}

// 4) NE of a new grid point
template <class Type>
inline Type phi_4iPlus2_4jPlus2(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/64.0*( 4.0*mat(i,j) + 5.0*mat(i+2,j+1) + 2*(mat(i+2,j) + mat(i,j+2)) + mat(i+2,j+2) + 25*mat(i+1,j+1) + 10.0*(mat(i+1,j) + mat(i,j+1)) + 5.0*mat(i+1,j+2) );
}

// 5) NW of an existing grid point
template <class Type>
inline Type phi_4iMin1_4j(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/16.0*( 3.0*(mat(i-1,j-1) + mat(i,j)) + 9.0*mat(i-1,j) + mat(i,j-1) );
}

// 6) SE of an existing grid point
template <class Type>
inline Type phi_4i_4jMin1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/16.0*( 3.0*(mat(i-1,j-1) + mat(i,j)) + mat(i-1,j) + 9.0*mat(i,j-1) );
}

// 7) SW of a new grid point on an existing vertical line
template <class Type>
inline Type phi_4iMin1_4jPlus1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( 3.0*mat(i-1,j-1) + 15.0*mat(i-1,j) + mat(i,j-1) + 5.0*mat(i,j) + 2.0*mat(i,j+1) + 6.0*mat(i-1,j+1) );
}

// 8) SW of a new grid point on an existing horizontal line
template <class Type>
inline Type phi_4iPlus1_4jMin1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( 3.0*mat(i-1,j-1) + mat(i-1,j) + 15*mat(i,j-1) + 5.0*mat(i,j) + 2.0*mat(i+1,j) + 6.0*mat(i+1,j-1) );
}

// 9) NW of a new grid point on an existing vertical line
template <class Type>
inline Type phi_4iMin1_4jPlus2(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( 6.0*mat(i-1,j) + 2.0*mat(i,j) + 15*mat(i-1,j+1) + 5.0*mat(i,j+1) + mat(i,j+2) + 3.0*mat(i-1,j+2) );
}

// 10) SE of a new grid point on an existing horizontal line
template <class Type>
inline Type phi_4iPlus2_4jMin1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( 6.0*mat(i,j-1) + 2.0*mat(i,j) + mat(i+2,j) + 15.0*mat(i+1,j-1) + 5.0*mat(i+1,j) + 3.0*mat(i+2,j-1) );
}

// 11) SE of a new grid point on an existing vertical line
template <class Type>
inline Type phi_4i_4jPlus1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( mat(i-1,j-1) + 5.0*mat(i-1,j) + 3.0*mat(i,j-1) + 15.0*mat(i,j) + 2.0*mat(i-1,j+1) + 6.0*mat(i,j+1) );
}

// 12) NW of a new grid point on an existing horizontal line
template <class Type>
inline Type phi_4iPlus1_4j(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( mat(i-1,j-1) + 3.0*mat(i-1,j) + 5.0*mat(i,j-1) + 15.0*mat(i,j) + 2.0*mat(i+1,j-1) + 6.0*mat(i+1,j) );
}

// 13) NE of a new grid point on an existing vertical line
template <class Type>
inline Type phi_4i_4jPlus2(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( 2.0*mat(i-1,j) + 6.0*mat(i,j) + 5.0*mat(i-1,j+1) + 15.0*mat(i,j+1) + mat(i-1,j+2) + 3.0*mat(i,j+2) );
}

// 14) NE of a new grid point on an existing horizontal line
template <class Type>
inline Type phi_4iPlus2_4j(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/32.0*( 2.0*mat(i,j-1) + 6.0*mat(i,j) + 5.0*mat(i+1,j-1) + 15.0*mat(i+1,j) + mat(i+2,j-1) + 3.0*mat(i+2,j) );
}

// 15) NW of a new grid point
template <class Type>
inline Type phi_4iPlus1_4jPlus2(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/64.0*( 2.0*(mat(i-1,j) + mat(i+1,j+2)) + 10.0*(mat(i,j) + mat(i+1,j+1)) + 5.0*(mat(i-1,j+1) + mat(i,j+2)) + 25.0*mat(i,j+1) + 4.0*mat(i+1,j) + mat(i-1,j+2) );
}

// 16) SE of a new grid point
template <class Type>
inline Type phi_4iPlus2_4jPlus1(GenMatrix<Type>& mat, int i, int j) {
  return 1.0/64.0*( 2.0*(mat(i,j-1) + mat(i+2,j+1)) + 10.0*(mat(i,j) + mat(i+1,j+1)) + 5.0*(mat(i+1,j-1) + mat(i+2,j)) + 4.0*mat(i,j+1) + 25.0*mat(i+1,j) + mat(i+2,j-1) );
}

#endif // _UNIFBSPLINES_H_
