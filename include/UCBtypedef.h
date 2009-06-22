//================================================================
//
// Created: May 2002
//                                                                           
// Author: Øyvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description:
//                                                                           
//================================================================
// Copyright (c) 2002 SINTEF Applied Mathematics
//================================================================
#ifndef _UCB_TYPEDEF_H_
#define _UCB_TYPEDEF_H_


/** \interface UCBspl_real
 *  \brief \#define UCBspl_real float; //Real type for large data objects (matrices)
 *
 *  Defines the real type for for large data objects such as vector and matrices.
 *  (Should be float or double) 
 */
#ifndef UCBspl_real
#define UCBspl_real float
#endif

template <class Type>
class GenMatrix; //<Type>;
/** \interface GenMatrixType 
 * \brief typedef GenMatrix<UCBspl_real> GenMatrixType;
 */
typedef GenMatrix<UCBspl_real> GenMatrixType;


#endif
