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

#ifndef POINT_ACCESS_UTILS_H
#define POINT_ACCESS_UTILS_H

#include <vector>

namespace UCBspl { // ??? temporary

// Temporary namespace scope here now. Later we may:
//   i) Skip namespace here and in cpp file
//  ii) include this file and GenMatrix in a file inside namespace scope
// iii) Also compile the cpp file inside the namespace scope


  // Operations on scattered data
  // ----------------------------
  // ascii file access
  int  numberOfPoints(const char filename[]);
  void readScatteredData(const char filename[], std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, int incr=1);
  void readScatteredDataBin(const char filename[], std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z, bool invZ = false);
  void readScatteredDataBin2(const char filename[], std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z);
  void printVRMLpoints(const char filename[], const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z,
		       double xmin, double ymin,
                       int incr=1,  // print every incr point (default=1)
                       double scale = 1.0);
  void printVTKpoints(const char filename[], const std::vector<double>& X, const std::vector<double>& Y, const std::vector<double>& Z,
                      double scale = 1.0);

  // binary file accesss
  void asciiXYZ2Bin(const char infile[], const char outfile[], int incr = 1);
  void asciiXYZ2Bin2(const char infile[], const char outfile[], int incr = 1);
  void readScatteredDataFileSetBin(const char metafile[], std::vector<double>& X, std::vector<double>& Y, std::vector<double>& Z);
  void grid2scat(const char infile[], const char outfile[]);

  //void printVTKtriangleStrips(const char filename[], const GenMatrix<UCBspl_real>& mat, double scale = 1.0);  

  void gisGeo2Planar(std::vector<double>& X, std::vector<double>& Y);

  //void exploreData(const char metafile[]); // Misc. on explore files

}; // end namespace

#endif

