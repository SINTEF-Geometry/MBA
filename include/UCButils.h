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
