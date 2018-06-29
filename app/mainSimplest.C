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

//! \page mainsimplest A complete but simple main program
/// \include mainSimplest.C
#include <MBA.h>
#include <UCButils.h>
#include <PointAccessUtils.h>


int main() {

  // Read scattered data from the file Data/rygg1.dat.
  // The format is assumed to be:
  // x y z
  // x y z
  // x y z
  // etc.

  typedef std::vector<double> dVec;
  std::shared_ptr<dVec> x_arr = std::make_shared<std::vector<double> >();
  std::shared_ptr<dVec> y_arr = std::make_shared<std::vector<double> >();
  std::shared_ptr<dVec> z_arr = std::make_shared<std::vector<double> >();
  UCBspl::readScatteredData("Data/rygg1.dat", *x_arr, *y_arr, *z_arr);

  MBA mba(x_arr, y_arr, z_arr);

  // Create spline surface.
  mba.MBAalg(1,1,7);

  UCBspl::SplineSurface surf = mba.getSplineSurface();

  // Find height and normal vector of surface in (x,y).
  double x = 5.0, y = 5.0;
  double z = surf.f(x,y);         
  double nx, ny, nz;
  surf.normalVector(x, y, nx, ny, nz);
  std::cout << "z-value in (5,5) = " << z << std::endl;
  std::cout << "Normal in (5,5) = (" << nx << "," << ny << "," << nz << ")"
	    << std::endl;

  // Sample surface and print to VRML file.
  UCBspl::printVRMLgrid("qwe.wrl", surf, 50, 50);  

  return 0;
}
