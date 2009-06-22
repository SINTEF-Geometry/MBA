//===========================================================================
//                                                                           
// File: filterRestrictTest.cpp                                              
//                                                                           
// Created: Fri Dec  5 11:20:11 2003                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: filterRestrictTest.cpp,v 1.1 2004-02-02 11:34:40 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
