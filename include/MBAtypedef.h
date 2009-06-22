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
