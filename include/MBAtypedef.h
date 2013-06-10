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

#ifndef _MBATYPEDEF_H
#define _MBATYPEDEF_H

/// For printing "debug" information to standard output
//#define MBA_DEBUG 1


#define MBA_UNDEFREAL 1572312345624422229996576879160.0

#include <GenMatrix.h>
#include <UCBtypedef.h>

#include <vector>

/** \interface dVec
* \brief typedef std::vector<double> dVec
*
*/
typedef std::vector<double> dVec;

/** \interface MBAbaseType
 * \brief enum MBAbaseType {MBA_ZERO, MBA_CONSTLS, MBA_CONSTVAL}
 *
 * This is the surface on which the spline surface is built over.
 * \see Details in MBA::setBaseType
 */
enum MBAbaseType {MBA_ZERO, MBA_CONSTLS, MBA_CONSTVAL}; // default is none

#endif
