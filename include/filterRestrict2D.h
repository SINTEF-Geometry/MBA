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

