//===========================================================================
//                                                                           
// File: checkWIN32andSGI.h                                                  
//                                                                           
// Created: Thu Nov 13 10:20:21 2003                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: checkWIN32andSGI.h,v 1.1 2004-02-02 11:34:40 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _CHECKWIN32ANDSGI_H
#define _CHECKWIN32ANDSGI_H

#ifdef WIN32
#define WIN32ORSGI
#endif

#ifdef SGI
#define WIN32ORSGI
#endif



#endif // _CHECKWIN32ANDSGI_H

