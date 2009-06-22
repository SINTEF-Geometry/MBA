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
#ifndef _UCB_UTILS_H_
#define _UCB_UTILS_H_

#include <UCBsplineSurface.h>
#include <vector>


#include <PointAccessUtils.h> // ??? Temporary here now so we dodn't need to include it in app.


/** \brief Misc. utilities, e.g., reading scattered data from file and
 *  printing surface to different plotting formats.
 *
 *  The functions are not documented in detail as many are self-explanatory.
 *  Consult the source code for details.
 */
namespace UCBspl {
  
  // Operations on UCBspl::SplineSurface, file access
  // ------------------------------------------------
  void printVRMLgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV, double scale = 1.0);
  void printVTKgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV, double scale = 1.0);
  void printVTKtriangleStrips(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV, double scale = 1.0);
  void printGNUgrid (const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV);  
  void printIRAPgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV);
  void printGLgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV);
  void printGLgridBin(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV,
                      const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z, // scattered data
                      double scale = 1.0);
  void saveSplineSurface(const char filename[], const UCBspl::SplineSurface& surf);
  void readSplineSurface(const char filename[], UCBspl::SplineSurface& surf);
  void saveSplineSurfaceBin(const char filename[], const UCBspl::SplineSurface& surf);
  void readSplineSurfaceBin(const char filename[], UCBspl::SplineSurface& surf);
  
}; // end namespace

#endif
