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

#include <vector>
#include <iostream>
#include "filterRestrict2D.h"

using namespace std;

int main()
{
    // making bogus 5x5 filter
    UCBspl_real data[] = {1,  2, 3,  2, 1,
			  2, -4, 6, -4, 2,
			  3,  6, 0,  6, 3,
			  2, -4, 6, -4, 2,
			  1,  2, 3,  2, 1};
    
    vector<UCBspl_real> result(25, 0);

    bool left = true;
    bool lower = true;
    int k = 2;
    int l = 2;
    filterRestrict2D(2, data, left, lower, k, l, &result[0]);

    int lower_bound_i = 0;
    int upper_bound_i = 5;
    int lower_bound_j = 0;
    int upper_bound_j = 5;

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
	    cout << result[pos++] << '\t';
	}
	cout << '\n';
    }
    
    return 1;
};
