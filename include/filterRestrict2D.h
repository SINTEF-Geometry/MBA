//===========================================================================
//                                                                           
// File: filterRestrict2D.h                                                  
//                                                                           
// Created: Thu Dec  4 17:17:28 2003                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: filterRestrict2D.h,v 1.1 2004-02-02 11:34:41 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FILTERRESTRICT2D_H
#define _FILTERRESTRICT2D_H

#include "UCBtypedef.h"

void filterRestrict2D(int m,                // width and height of filter (2m+1)
		      const UCBspl_real* F, // complete 2D filter F(i,j) = F[j * (2m+1) + i]
		      bool left,            // true = left boundary (index i)
		      bool lower,           // true = lower boundary, (index j)
		      int k,                // number of depassing coefficients in i direction
		      int l,                // number of depassing coefficients in j direction
		      UCBspl_real* result); // pointer to memory area where the result is written
		      
// NB: k and l must both be in [0, m].  If both are equal to 0, then the unrestricted operator
// is returned.  In general, the result operator will have the resolution (2m+1-k, 2m+1-l).
// The original filter and the result are both stored row-wise (i has lowest stride).

#endif // _FILTERRESTRICT2D_H

