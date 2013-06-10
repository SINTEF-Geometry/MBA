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

#include <UCBsplines.h>

#ifdef DEBUG_UCBspl
#include <iostream>
#endif

// Refine lattice for C2
// =====================
void UCBspl::refineCoeffsC2(const GenMatrix<UCBspl_real>& PSI, GenMatrix<UCBspl_real>& PSIprime) {
  
  int mm = PSI.noX()-3; // corresponds to current m_ in MBA
  int nn = PSI.noY()-3; // corresponds to current n_ in MBA
  
  PSIprime.resize(2*mm+3, 2*nn+3);
  
  int i,j;
  
  for (j = 0; j <= nn; j++) {
    // calculate in first column, sw and left of current
    PSIprime(-1,2*j-1) = phi_2iPlus1_2jPlus1(PSI,-1,j-1);   // 4) a new crossing in both directions
    PSIprime(-1,2*j)   = phi_2iPlus1_2j(PSI,-1,j);          // 3) On an existing horizontal line
    
    for (i = 0; i <= mm; i++) {
      PSIprime(2*i,2*j)     = phi_2i_2j(PSI,i,j);           // 1) On an existing vertex
      PSIprime(2*i,2*j+1)   = phi_2i_2jPluss1(PSI,i,j);     // 2) On an existing vertical line
      PSIprime(2*i+1,2*j)   = phi_2iPlus1_2j(PSI,i,j);      // 3) On an existing horizontal line
      PSIprime(2*i+1,2*j+1) = phi_2iPlus1_2jPlus1(PSI,i,j); // 4) a new crossing in both directions
    }
  }
  
  // Upper left corner
  PSIprime(-1,2*nn+1) = phi_2iPlus1_2jPlus1(PSI,-1,nn);     // 4) a new crossing in both directions
  
  // Lower row except lower left corner
  for (i = 0; i <= mm; i++) {
    PSIprime(2*i,-1)   = phi_2i_2jPluss1(PSI,i,-1);      // 2) On an existing vertical line
    PSIprime(2*i+1,-1) = phi_2iPlus1_2jPlus1(PSI,i,-1);  // 4) a new crossing in both directions
  }
}


// Refine lattice for C1
// =====================
void UCBspl::refineCoeffsC1(const GenMatrix<UCBspl_real>& PSI, GenMatrix<UCBspl_real>& PSIprime) {
  
  int mm = (PSI.noX()-2)/2; // corresponds to current m_ in MBA
  int nn = (PSI.noY()-2)/2; // corresponds to current n_ in MBA
  
  PSIprime.resize(4*mm+2, 4*nn+2);
  
  int i,j;
  for (j = 0; j < nn; j++) {
    
    int jj = 2*j;
    for (i = 0; i < mm; i++) {
      int ii = 2*i;
      
      // Around existing grid point
      PSIprime(4*i-1,4*j-1) = phi_4iMin1_4jMin1(PSI,ii,jj);
      PSIprime(4*i-1,4*j) = phi_4iMin1_4j(PSI,ii,jj);
      PSIprime(4*i,4*j-1) = phi_4i_4jMin1(PSI,ii,jj);
      PSIprime(4*i,4*j) = phi_4i_4j(PSI,ii,jj);
      
      // Around midpoint of existing vertical line
      PSIprime(4*i-1,4*j+1) = phi_4iMin1_4jPlus1(PSI,ii,jj);
      PSIprime(4*i-1,4*j+2) = phi_4iMin1_4jPlus2(PSI,ii,jj);
      PSIprime(4*i,4*j+1) = phi_4i_4jPlus1(PSI,ii,jj);
      PSIprime(4*i,4*j+2) = phi_4i_4jPlus2(PSI,ii,jj);
      
      // Around midpoint of existing horizontal line
      PSIprime(4*i+1,4*j-1) = phi_4iPlus1_4jMin1(PSI,ii,jj);
      PSIprime(4*i+1,4*j) = phi_4iPlus1_4j(PSI,ii,jj);
      PSIprime(4*i+2,4*j-1) = phi_4iPlus2_4jMin1(PSI,ii,jj);
      PSIprime(4*i+2,4*j) = phi_4iPlus2_4j(PSI,ii,jj);
      
      // Around new grid point
      PSIprime(4*i+1,4*j+1) = phi_4iPlus1_4jPlus1(PSI,ii,jj);
      PSIprime(4*i+1,4*j+2) = phi_4iPlus1_4jPlus2(PSI,ii,jj);
      PSIprime(4*i+2,4*j+1) = phi_4iPlus2_4jPlus1(PSI,ii,jj);
      PSIprime(4*i+2,4*j+2) = phi_4iPlus2_4jPlus2(PSI,ii,jj);
    }
    
    // Around existing grid point of last vertical line
    int ii = 2*mm;
    PSIprime(4*mm-1,4*j-1) = phi_4iMin1_4jMin1(PSI,ii,jj);
    PSIprime(4*mm-1,4*j) = phi_4iMin1_4j(PSI,ii,jj);
    PSIprime(4*mm,4*j-1) = phi_4i_4jMin1(PSI,ii,jj);
    PSIprime(4*mm,4*j) = phi_4i_4j(PSI,ii,jj);
    
    // Around midpoint of last existing vertical line
    PSIprime(4*mm-1,4*j+1) = phi_4iMin1_4jPlus1(PSI,ii,jj);
    PSIprime(4*mm-1,4*j+2) = phi_4iMin1_4jPlus2(PSI,ii,jj);
    PSIprime(4*mm,4*j+1) = phi_4i_4jPlus1(PSI,ii,jj);
    PSIprime(4*mm,4*j+2) = phi_4i_4jPlus2(PSI,ii,jj);
  }
  
  // Around last existing horizontal line
  // ------------------------------------
  int jj = 2*nn;
  for (i = 0; i < mm; i++) {
    // Around existing grid point 
    int ii = 2*i;
    PSIprime(4*i-1,4*nn-1) = phi_4iMin1_4jMin1(PSI,ii,jj);
    PSIprime(4*i-1,4*nn)   = phi_4iMin1_4j(PSI,ii,jj);
    PSIprime(4*i,4*nn-1)   = phi_4i_4jMin1(PSI,ii,jj);
    PSIprime(4*i,4*nn)     = phi_4i_4j(PSI,ii,jj);
    
    // Around midpoint
    PSIprime(4*i+1,4*nn-1) = phi_4iPlus1_4jMin1(PSI,ii,jj);
    PSIprime(4*i+1,4*nn)   = phi_4iPlus1_4j(PSI,ii,jj);
    PSIprime(4*i+2,4*nn-1) = phi_4iPlus2_4jMin1(PSI,ii,jj);
    PSIprime(4*i+2,4*nn)   = phi_4iPlus2_4j(PSI,ii,jj);
  }
  
  // Upper right corner (around existing grid point)
  int ii = 2*mm;
  // int jj = 2*nn from above
  PSIprime(4*mm-1,4*nn-1) = phi_4iMin1_4jMin1(PSI,ii,jj);
  PSIprime(4*mm-1,4*nn)   = phi_4iMin1_4j(PSI,ii,jj);
  PSIprime(4*mm,4*nn-1)   = phi_4i_4jMin1(PSI,ii,jj);
  PSIprime(4*mm,4*nn)     = phi_4i_4j(PSI,ii,jj);
}


// Restriction operator for uniform cubic C2 splines
// (Designed as the "transposed" of the corresponding refinement operator)
// =======================================================================
bool UCBspl::restrictCoeffsC2(const GenMatrix<UCBspl_real>& rr, GenMatrix<UCBspl_real>& r) {
  
  int old_noX = rr.noX();
  int old_noY = rr.noY();
  
  
  if ((old_noX-3)%2 != 0 || (old_noY-3)%2 != 0) {
#ifdef DEBUG_UCBspl
    cout << "ERROR, invalid grid for projection; see comments in the code above: " << old_noX << " " << old_noY << endl;
    exit(-1);
#endif
    return false;
  }
  
  // Size of grid at the new coarser level
  int noX = (old_noX-3)/2 + 3;
  int noY = (old_noY-3)/2 + 3;
  r.resize(noX, noY);
  // r.fill(1.0e+9f); just to test if all elements are set
  
  double denominatorFull   = 256.0; // 64.
  double denominatorCorner =  25.0; // (op3)
  double denominatorOp2    = 225.0;
  double denominatorOp4    =  80.0;
  double denominatorOp5    = 240.0;
  double denominatorOp6    =  75.0;
  
  int kk,ll,kk_new,ll_new;
  
  // The interior full operator (1,...,31)
  for (kk_new = 1; kk_new <= noX-4; kk_new++) {
    for (ll_new = 1; ll_new <= noY-4; ll_new++) {
      kk = 2*kk_new; // corresponding indices in old matrix
      ll = 2*ll_new;      
      double val = rr(kk-2,ll-2) + rr(kk-2,ll+2) + rr(kk+2,ll-2) + rr(kk+2,ll+2) // outer corners
        + 16.0 * (rr(kk-1,ll-1) + rr(kk-1,ll+1) + rr(kk+1,ll-1) + rr(kk+1,ll+1)) // inner corners
        +  6.0 * (rr(kk,ll-2) + rr(kk,ll+2) + rr(kk-2,ll) + rr(kk+2,ll)) // outer sides
        + 24.0 * (rr(kk,ll-1) + rr(kk,ll+1) + rr(kk-1,ll) + rr(kk+1,ll)) // inner sides
        +  4.0 * (rr(kk-2,ll-1) + rr(kk-1,ll-2) + rr(kk+1,ll-2) + rr(kk+2,ll-1)) // remains llower hallf
        +  4.0 * (rr(kk-2,ll+1) + rr(kk-1,ll+2) + rr(kk+1,ll+2) + rr(kk+2,ll+1)) // remains uppe hallf
        + 36.0 * rr(kk,ll);
      r(kk_new,ll_new) = val/denominatorFull;
    }
  }
  
  
  // Upper and lower boundary where full width of operator
  // Here kk = 2*kk_new
  // -----------------------------------------------------
  for (kk_new = 1; kk_new <= noX-4; kk_new++) {
    
    double val;
    kk = 2*kk_new;
    
    // Upper boundary where full width of operator (=op4)
    ll_new = noY-2;
    ll = 2*ll_new-1;
    val = rr(kk-2,ll-1) + rr(kk+2,ll-1) // outer corners
      + 16.0 * (rr(kk-1,ll) + rr(kk+1,ll)) // inner corners
      +  6.0 *  rr(kk,ll-1) // outer sides
      + 24.0 *  rr(kk,ll)   // inner sides (but now center
      +  4.0 * (rr(kk-2,ll) + rr(kk+2,ll)) // remains
      +  4.0 * (rr(kk-1,ll-1) + rr(kk+1,ll-1)); // remains
    r(kk_new,ll_new) = val/denominatorOp4;
    
    // Next to upper boundary where full width of operator (=op5)
    ll_new = noY-3;
    ll = 2*ll_new;
    val = rr(kk-2,ll-2) + rr(kk+2,ll-2) // outer corners
      + 16.0 * (rr(kk-1,ll-1) + rr(kk-1,ll+1) + rr(kk+1,ll-1) + rr(kk+1,ll+1)) // inner corners
      +  6.0 * (rr(kk,ll-2) + rr(kk-2,ll) + rr(kk+2,ll)) // outer sides
      + 24.0 * (rr(kk,ll-1) + rr(kk,ll+1) + rr(kk-1,ll) + rr(kk+1,ll)) // inner sides
      +  4.0 * (rr(kk-2,ll-1) + rr(kk-1,ll-2) + rr(kk+1,ll-2) + rr(kk+2,ll-1)) // remains llower hallf
      +  4.0 * (rr(kk-2,ll+1) + rr(kk+2,ll+1)) // remains uppe hallf
      + 36.0 * rr(kk,ll);
    r(kk_new,ll_new) = val/denominatorOp5;
    
    // Lower boundary where full width of operator (=op4, with the two rows interchanged)
    ll_new = ll = -1;
    val = rr(kk-2,ll+1) + rr(kk+2,ll+1) // outer corners
      + 16.0 * (rr(kk-1,ll) + rr(kk+1,ll)) // inner corners
      +  6.0 *  rr(kk,ll+1) // outer sides
      + 24.0 *  rr(kk,ll)   // inner sides (but now center
      +  4.0 * (rr(kk-2,ll) + rr(kk+2,ll)) // remains
      +  4.0 * (rr(kk-1,ll+1) + rr(kk+1,ll+1)); // remains
    r(kk_new,ll_new) = val/denominatorOp4;
    
    // Next to lower boundary where full width of operator (=op5, with the lower row moved to top)
    ll_new = ll = 0;
    val = rr(kk-2,ll+2) + rr(kk+2,ll+2) // outer corners
      + 16.0 * (rr(kk-1,ll-1) + rr(kk-1,ll+1) + rr(kk+1,ll-1) + rr(kk+1,ll+1)) // inner corners
      +  6.0 * (rr(kk,ll+2) + rr(kk-2,ll) + rr(kk+2,ll)) // outer sides
      + 24.0 * (rr(kk,ll-1) + rr(kk,ll+1) + rr(kk-1,ll) + rr(kk+1,ll)) // inner sides
      +  4.0 * (rr(kk-2,ll+1) + rr(kk-1,ll+2) + rr(kk+1,ll+2) + rr(kk+2,ll+1)) // remains upper hallf
      +  4.0 * (rr(kk-2,ll-1) + rr(kk+2,ll-1)) // remains llower
      + 36.0 * rr(kk,ll);
    r(kk_new,ll_new) = val/denominatorOp5;
  } 
  
  
  
  // Right and left boundary where full width of operator
  // Here ll = 2*ll_new
  // -----------------------------------------------------
  for (ll_new = 1; ll_new <= noY-4; ll_new++) {
    
    double val;
    ll = 2*ll_new;
    
    // Right boundary where full width of operator (=op4 rotated clockwise)
    kk_new = noX-2;
    kk = 2*kk_new-1;
    val = rr(kk-1,ll+2) + rr(kk-1,ll-2) // outer corners
      + 16.0 * (rr(kk,ll+1) + rr(kk,ll-1)) // inner corners
      +  6.0 *  rr(kk-1,ll) // outer sides
      + 24.0 *  rr(kk,ll)   // inner sides (but now center
      +  4.0 * (rr(kk,ll+2) + rr(kk,ll-2)) // remains
      +  4.0 * (rr(kk-1,ll+1) + rr(kk-1,ll-1)); // remains
    r(kk_new,ll_new) = val/denominatorOp4;
    
    // Next to right boundary where full width of operator (=op5 rotated clockwise)
    kk_new = noX-3;
    kk = 2*kk_new;
    val = rr(kk-2,ll-2) + rr(kk-2,ll+2) // outer corners
      + 16.0 * (rr(kk-1,ll-1) + rr(kk-1,ll+1) + rr(kk+1,ll-1) + rr(kk+1,ll+1)) // inner corners
      +  6.0 * (rr(kk,ll+2) + rr(kk,ll-2) + rr(kk-2,ll)) // outer sides
      + 24.0 * (rr(kk,ll-1) + rr(kk,ll+1) + rr(kk-1,ll) + rr(kk+1,ll)) // inner sides
      +  4.0 * (rr(kk-1,ll+2) + rr(kk-1,ll-2) + rr(kk-2,ll+1) + rr(kk-2,ll-1)) // remains lleft
      +  4.0 * (rr(kk+1,ll+2) + rr(kk+1,ll-2)) // remains right
      + 36.0 * rr(kk,ll);
    r(kk_new,ll_new) = val/denominatorOp5;
    
    // Left boundary where full width of operator (=op4, transposed)
    kk_new = kk = -1;
    val = rr(kk+1,ll+2) + rr(kk+1,ll-2) // outer corners
      + 16.0 * (rr(kk,ll+1) + rr(kk,ll-1)) // inner corners
      +  6.0 *  rr(kk+1,ll) // outer sides
      + 24.0 *  rr(kk,ll)   // inner sides (but now center
      +  4.0 * (rr(kk,ll+2) + rr(kk,ll-2)) // remains
      +  4.0 * (rr(kk+1,ll+1) + rr(kk+1,ll-1)); // remains
    r(kk_new,ll_new) = val/denominatorOp4;
    
    // Next to left boundary where full width of operator (=op5, transposed)
    kk_new = kk = 0;
    val = rr(kk+2,ll+2) + rr(kk+2,ll-2) // outer corners
      + 16.0 * (rr(kk-1,ll-1) + rr(kk-1,ll+1) + rr(kk+1,ll-1) + rr(kk+1,ll+1)) // inner corners
      +  6.0 * (rr(kk,ll+2) + rr(kk,ll-2) + rr(kk+2,ll)) // outer sides
      + 24.0 * (rr(kk,ll-1) + rr(kk,ll+1) + rr(kk-1,ll) + rr(kk+1,ll)) // inner sides
      +  4.0 * (rr(kk+1,ll+2) + rr(kk+1,ll-2) + rr(kk+2,ll+1) + rr(kk+2,ll-1)) // remains right
      +  4.0 * (rr(kk-1,ll+2) + rr(kk-1,ll-2)) // remains lleft
      + 36.0 * rr(kk,ll);
    r(kk_new,ll_new) = val/denominatorOp5;
  }
  
  // Set the four nodes in each corner.
  // ==================================
  
  // For example, the lower left nodes have indices:
  // have indices (kk,ll,kk_new and ll_new):
  // (-1, 0), (0, 0)
  // (-1,-1), (0,-1)
  
  // See the report for calculation of the different operators
  
  // CORNERS (op3)
  // -------------
  // lower left corner
  double val = 16.0*rr(-1,-1) + 4.0*(rr(-1,0)+rr(0,-1)) + rr(0,0);
  r(-1,-1) = val/denominatorCorner;
  
  // upper right corner
  val = 16.0*rr(old_noX-2,old_noY-2) + 4.0*(rr(old_noX-3,old_noY-2)+rr(old_noX-2,old_noY-3)) 
    + rr(old_noX-3,old_noY-3);
  r(noX-2,noY-2) = val/denominatorCorner;
  
  // lower right corner
  val = 16.0*rr(old_noX-2,-1) + 4.0*(rr(old_noX-2,0)+rr(old_noX-3,-1)) + rr(old_noX-3,0);
  r(noX-2,-1) = val/denominatorCorner;
  
  // upper left corner
  val = 16.0*rr(-1,old_noY-2) + 4.0*(rr(-1,old_noY-3)+rr(0,old_noY-2)) + rr(0,old_noY-3);
  r(-1,noY-2) = val/denominatorCorner;
  
  
  // NEXT to the corner, symmetric (op2)
  // -----------------------------------
  // upper left
  val = 36.0 * rr(0,old_noY-3)
    +16.0 * (rr(-1,old_noY-2) + rr(-1,old_noY-4) + rr(1,old_noY-2) + rr(1,old_noY-4))
    +24.0 * (rr(-1,old_noY-3) + rr(1,old_noY-3) + rr(0,old_noY-2) + rr(0,old_noY-4))
    + 6.0 * (rr(0,old_noY-5)  + rr(2,old_noY-3))
    + 4.0 * (rr(-1,old_noY-5) + rr(1,old_noY-5) + rr(2,old_noY-4) + rr(2,old_noY-2));
  r(0,noY-3) = val/denominatorOp2;
  
  // lower left
  val = 36.0 * rr(0,0)
    +16.0 * (rr(-1,-1) + rr(-1,1) + rr(1,-1) + rr(1,1))
    +24.0 * (rr(-1,0) + rr(1,0) + rr(0,-1) + rr(0,1))
    + 6.0 * (rr(0,2)  + rr(2,0))
    + 4.0 * (rr(2,-1) + rr(2,1) + rr(1,2) + rr(-1,2));
  r(0,0) = val/denominatorOp2;
  
  // upper right
  val = 36.0 * rr(old_noX-3,old_noY-3)
    +16.0 * (rr(old_noX-4,old_noY-4) + rr(old_noX-4,old_noY-2) + rr(old_noX-2,old_noY-4) + rr(old_noX-2,old_noY-2))
    +24.0 * (rr(old_noX-4,old_noY-3) + rr(old_noX-3,old_noY-4) + rr(old_noX-2,old_noY-3) + rr(old_noX-3,old_noY-2))
    + 6.0 * (rr(old_noX-5,old_noY-3)  + rr(old_noX-3,old_noY-5))
    + 4.0 * (rr(old_noX-5,old_noY-2) + rr(old_noX-5,old_noY-4) + rr(old_noX-4,old_noY-5) + rr(old_noX-2,old_noY-5));
  r(noX-3,noY-3) = val/denominatorOp2;
  
  // lower right
  val = 36.0 * rr(old_noX-3,0)
    +16.0 * (rr(old_noX-4,-1) + rr(old_noX-2,-1) + rr(old_noX-2,1) + rr(old_noX-4,1))
    +24.0 * (rr(old_noX-4,0) + rr(old_noX-3,-1) + rr(old_noX-2,0) + rr(old_noX-3,1))
    + 6.0 * (rr(old_noX-5,0)  + rr(old_noX-3,2))
    + 4.0 * (rr(old_noX-5,-1) + rr(old_noX-5,1) + rr(old_noX-4,2) + rr(old_noX-2,2));
  r(noX-3,0) = val/denominatorOp2;
  
  
  // The unsymmetric operator next to the corner (op6)
  // --------------------------------------------------
  //
  // upper left-up:
  // 16 <24> 16 4
  //  4   6   4 1
  val = 24.0 * rr(0,old_noY-2)
    +16.0 * (rr(-1,old_noY-2) + rr(1,old_noY-2))
    + 4.0 * (rr(-1,old_noY-3) + rr(1,old_noY-3) + rr(2,old_noY-2))
    + 6.0 *  rr(0,old_noY-3)
    +       rr(2,old_noY-3);
  r(0,noY-2) = val/denominatorOp6;
  
  // upper left-left:
  //  16  4
  // <24> 6
  //  16  4
  //   4  1
  val = 24.0 * rr(-1,old_noY-3)
    +16.0 * (rr(-1,old_noY-2) + rr(-1,old_noY-4))
    + 4.0 * (rr(-1,old_noY-5) + rr(0,old_noY-4) + rr(0,old_noY-2))
    + 6.0 *  rr(0,old_noY-3)
    +       rr(0,old_noY-5);
  r(-1,noY-3) = val/denominatorOp6;
  
  // lower left-left
  //   4  1
  //  16  4
  // <24> 6
  //  16  4
  val = 24.0 * rr(-1,0)
    +16.0 * (rr(-1,-1) + rr(-1,1))
    + 4.0 * (rr(0,-1) + rr(0,1) + rr(-1,2))
    + 6.0 *  rr(0,0)
    +       rr(0,2);
  r(-1,0) = val/denominatorOp6;
  
  // lower left-down
  //  4   6   4 1
  // 16 <24> 16 4
  val = 24.0 * rr(0,-1)
    +16.0 * (rr(-1,-1) + rr(1,-1))
    + 4.0 * (rr(-1,0) + rr(1,0) + rr(2,-1))
    + 6.0 *  rr(0,0)
    +       rr(2,0);
  r(0,-1) = val/denominatorOp6;
  
  
  // upper right-up:
  // 4 16 <24> 16
  // 1  4   6   4
  val = 24.0 * rr(old_noX-3,old_noY-2)
    +16.0 * (rr(old_noX-4,old_noY-2) + rr(old_noX-2,old_noY-2))
    + 4.0 * (rr(old_noX-5,old_noY-2) + rr(old_noX-4,old_noY-3) + rr(old_noX-2,old_noY-3))
    + 6.0 *  rr(old_noX-3,old_noY-3)
    +       rr(old_noX-5,old_noY-3);
  r(noX-3,noY-2) = val/denominatorOp6;
  
  
  // upper right-right
  // 4  16
  // 6 <24>
  // 4  16
  // 1   4
  val = 24.0 * rr(old_noX-2,old_noY-3)
    +16.0 * (rr(old_noX-2,old_noY-2) + rr(old_noX-2,old_noY-4))
    + 4.0 * (rr(old_noX-3,old_noY-2) + rr(old_noX-3,old_noY-4) + rr(old_noX-2,old_noY-5))
    + 6.0 *  rr(old_noX-3,old_noY-3)
    +       rr(old_noX-3,old_noY-5);
  r(noX-2,noY-3) = val/denominatorOp6;
  
  
  // lower right-down
  // 1  4   6   4
  // 4 16 <24> 16
  val = 24.0 * rr(old_noX-3,-1)
    +16.0 * (rr(old_noX-4,-1) + rr(old_noX-2,-1))
    + 4.0 * (rr(old_noX-5,-1) + rr(old_noX-4,0) + rr(old_noX-2,0))
    + 6.0 *  rr(old_noX-3,0)
    +       rr(old_noX-5,0);
  r(noX-3,-1) = val/denominatorOp6;
  
  // lower right-right
  // 1   4
  // 4  16
  // 6 <24>
  // 4  16
  val = 24.0 * rr(old_noX-2,0)
    +16.0 * (rr(old_noX-2,-1) + rr(old_noX-2,1))
    + 4.0 * (rr(old_noX-3,-1) + rr(old_noX-3,1) + rr(old_noX-2,2))
    + 6.0 *  rr(old_noX-3,0)
    +       rr(old_noX-3,2);
  r(noX-2,0) = val/denominatorOp6;

  return true;
}
