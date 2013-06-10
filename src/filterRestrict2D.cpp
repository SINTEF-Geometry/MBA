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


#include <stdexcept>
#include <vector>
#include "filterRestrict2D.h"

using namespace std;

namespace {
    void argument_check(int m, const UCBspl_real* F, int k, int l, const UCBspl_real* result);
}; // end anonymous namespace


//===========================================================================
void filterRestrict2D(int m,                // width and height of filter (2m+1)
		      const UCBspl_real* F, // complete 2D filter F(i,j) = F[j * (2m+1) + i]
		      bool left,            // true = left boundary (index i)
		      bool lower,           // true = lower boundary, (index j)
		      int k,                // number of depassing coefficients in i direction
		      int l,                // number of depassing coefficients in j direction
		      UCBspl_real* result)  // pointer to memory area where the result is written
//===========================================================================
{
    argument_check(m, F, k, l, result);
    // we now suppose that the arguments are valid

    // making temporary copy of filter
    int rlen = 2 * m + 1;
    std::vector<UCBspl_real> Ftemp(F, F + rlen * rlen);

    // modifying borders
    if (left) {
	for (int j = 0; j < rlen; ++j) {
	    for (int i = 1; i <= k; ++i) {
		UCBspl_real temp = Ftemp[j * rlen + k - i];
		Ftemp[j * rlen + k] += (i+1) * temp;
		Ftemp[j * rlen + k + 1] += (-i) * temp;
	    }
	} 
    } else {
	for (int j = 0; j < rlen; ++j) {
	    for (int i = 1; i <= k; ++i) {
		UCBspl_real temp = Ftemp[j * rlen + (rlen - 1) - k + i];
		Ftemp[j * rlen + (rlen - 1) - k] += (i+1) * temp;
		Ftemp[j * rlen + (rlen - 1) - k - 1] += (-i) * temp;
	    }
	}
    }
    if (lower) {
	for (int i = 0; i < rlen; ++i) {
	    for (int j = 1; j <= l; ++j) {
		UCBspl_real temp = Ftemp[(l - j) * rlen + i];
		Ftemp[l * rlen + i] += (j+1) * temp;
		Ftemp[(l+1) * rlen + i] += (-j) * temp;
	    }
	}
    } else {
	for (int i = 0; i < rlen; ++i) {
	    for (int j = 1; j <= l; ++j) {
		UCBspl_real temp = Ftemp[(rlen - 1 - l + j) * rlen + i];
		Ftemp[(rlen - 1 - l) * rlen + i] += (j+1) * temp;
		Ftemp[(rlen - 1 - 1 - l) * rlen + i] += (-j) * temp;
	    }
	}
    }


    // copy the relevant parts of the filter to the result array
    int lower_bound_i = 0;
    int upper_bound_i = rlen;
    int lower_bound_j = 0;
    int upper_bound_j = rlen;
    if (left) {
	lower_bound_i += k;
    } else {
	// right
	upper_bound_i -= k;
    }
    if (lower) {
	lower_bound_j += l;
    } else {
	// upper
	upper_bound_j -= l;
    }
    
    int pos = 0;
    for (int j = lower_bound_j; j < upper_bound_j; ++j) {
	for (int i = lower_bound_i; i < upper_bound_i; ++i) {
	    result[pos++] = Ftemp[j * rlen + i];
	}
    }
}

namespace {

//===========================================================================
void argument_check(int m, const UCBspl_real* F, int k, int l, const UCBspl_real* result)
//===========================================================================
{
    if (m < 1) {
	throw runtime_error("Invalid size of filter.");
    }
    if (!F || !result) {
	throw runtime_error("Null pointers detected.");
    }
    if (k < 0 || k > m) {
	throw runtime_error("Invalid value for k (number of steps outside boundary).  "
			    "Should be between 0 and m.");
    }
    if (k < 0 || k > m) {
	throw runtime_error("Invalid value for l (number of steps outside boundary).  "
			    "Should be between 0 and m.");
    }
}

}; //end anonymous namespace 
