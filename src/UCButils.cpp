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

#include <UCButils.h>

#ifdef MBA_DEBUG
#include <MBAclock.h>
#endif

#ifdef WIN32
#include <minmax.h>
#endif

#include <iostream>
#include <fstream>
//#include <cmath>
#include <math.h>
#include <algorithm>
using namespace std;

#include <stdio.h>
#include <string>



void UCBspl::printVRMLgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV, double scale) {
  // Print out a sample grid of points of an MBA surface to VRML format.
  // The dimensions of the grid are noU and noV.
  // Number of faces = noU * noV.

  
#ifdef MBA_DEBUG
  cout << "printSampleVRML to: " << filename << endl; 
#endif
  
  bool bscale = false;
  if (scale != 1.0)
    bscale = true;
  
  ofstream os(filename);
  
  os << "#VRML V1.0 ascii" << endl;
  os << endl;
  os << "Separator {" << endl; // Separator 1
  os << "    ShapeHints {" << endl;
  os << "       vertexOrdering  COUNTERCLOCKWISE" << endl;
  os << "       shapeType       SOLID" << endl;
  os << "       faceType        CONVEX" << endl;
  os << "       creaseAngle     30.0" << endl;
  os << "    }" << endl; // ShapeHints
  
  os << "    Separator {" << endl; // Separator 2
  os << "       Coordinate3 {" << endl;
  os << "           point       [" << endl;
  
  double umin = surf.umin();
  double vmin = surf.vmin();
  double du = (surf.umax() - umin)/(double)(noU-1);
  double dv = (surf.vmax() - vmin)/(double)(noV-1);  
    
  // VERTICES
  int i,j;
  for (j=0; j<noV; j++) {
    double v = vmin + (double)j*dv;
    for(i=0; i<noU; i++) {
      double u = umin + (double)i*du;
      double z = surf.f(u,v);
      if (bscale)
        z *= scale;
      
      os << u - umin << " " << v - vmin << " " << z; // Transform to origin
      
      if(i < noU-1 || j < noV-1) 
        os << ",\n";
    }
  }
  os << "]\n";
  os << "}\n"; // Coordinate3
  
  os << " IndexedFaceSet {\n";
  os << "    coordIndex [\n" << flush;
  
  //FACES
  for(j=0; j<noV-1; j++) {
    for(i=0; i<noU-1; i++) {
      os << j*noU+i << ", " << j*noU+i+1 << ", " << (j+1)*noU+i+1 << ", " << (j+1)*noU+i << ", " << -1 << "," << endl;
    }
    os << "\n";
  }
  os << "       ]\n"; // coordIndex
  os << "    }\n"; // IndexedFaceSet
  os << "  }\n"; // Separator 2
  
  /*
  if (displayScatteredData) {
    // The scattered data
    os << "  Separator {" << endl;
    // R G B color
    os << " Material { emissiveColor " << 1.0 << " " << 0.0 << " " << 0.0 << " }" << endl;
    
    os << "     Coordinate3 {" << endl;
    os << "        point [" << endl;
    
    int noPoints = mba.getData().size();
    const vector<double>& U = *mba.getData().U(); 
    const vector<double>& V = *mba.getData().V();
    const vector<double>& Zorig = *mba.getData().Zorig();
    for (int i = 0; i < noPoints; i++) {
      double z_orig = Zorig[i];
      if (bscale)
        z_orig *= scale;
      
      os << U[i]-umin << " " << V[i]-vmin << " " << z_orig;
      if (i == noPoints-1)
        os << endl;
      else
        os << "," << endl;
    }
    os << "          ]" << endl; // point
    os << "     }" << endl; // Coordinate3
    os << "     PointSet {" << endl;
    os << "        startIndex 0" << endl;
    os << "        numPoints  -1" << endl;
    os << "     }" << endl; // PointSet
    os << "  }" << endl; // Separator
  } // displayScatteredData
  */
  
  os << "}\n"; // Separator 1
}

void UCBspl::printGNUgrid (const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV) {
  
#ifdef MBA_DEBUG
  cout << "Printing grid to: " << filename << endl;
#endif
  
  double du = (surf.umax() - surf.umin())/(double)(noU-1);
  double dv = (surf.vmax() - surf.vmin())/(double)(noV-1);  
  ofstream ofile(filename);
  
  for (int j = 0; j < noV; j++) {      
    double v = surf.vmin() + j*dv;
    for (int i = 0; i < noU; i++) {
      double u = surf.umin() + i*du;
      
      ofile << surf.f(u,v) << "\n";
    }
    ofile << endl;
  }  
}


void UCBspl::printVTKtriangleStrips(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV, double scale) {
  
  cout << "Printing grid to vtk poly data file with triangle strips...." << endl;

  double umin = surf.umin();
  double vmin = surf.vmin();
  double umax = surf.umax();
  double vmax = surf.vmax();
  double du = (umax - umin)/(double)(noU-1);
  double dv = (vmax - vmin)/(double)(noV-1);  
  int noPoints = noU*noV;

  ofstream os(filename);

  os << "# vtk DataFile Version 2.0" << endl;
  os << "This file was generated by class UCButils (triangle strips)" << endl;
  os << "ASCII" << endl;
  os << "DATASET POLYDATA" << endl;
  os << "POINTS " << noPoints << " float\n";

  // The points row by row
  int i,j;
  for (j = 0; j < noV; j++) {
    double v = vmin + j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + i*du;      
      os << u << " " << v << " " << surf.f(u,v)*scale << '\n';
    }
  }

  int noStrips = noV-1;
  int size = (1 + noU*2)*(noV - 1);
  os << "TRIANGLE_STRIPS " << noStrips << " " << size << endl;

  const int izoff = noU + 1; // offset in 1D array to north-east
  int indz = 0;
  for (j = 0; j < noV-1; j++) {
    os << 2*noU << " ";  // *ptr_tris = 2*dim[0]; // As in the file format, no. of nodes 
    //ptr_tris++;
    os << indz+noU << " "; // *ptr_tris = indz + dim[0]; // first point in strip
    //ptr_tris++;
    for (int i = 0; i < noU-1; i++) {
			
      // lower left and upper right point of each diagonal
      os << indz << " "; // *ptr_tris = indz;

      //ptr_tris++;
      os << indz + izoff << " ";// *ptr_tris = indz + izoff;
      //ptr_tris++;
			indz++;
    }
		os << indz << '\n';// *ptr_tris = indz;  // last point in strip
    indz++;
    //ptr_tris++;
  }

  // And then the normals
  double gx,gy,gz;
  os << "POINT_DATA " << noPoints << '\n';
  os << "NORMALS normals float" << '\n';
  for (j = 0; j < noV; j++) {
    double v = vmin + j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + i*du;
      surf.normalVector(u,v,gx,gy,gz);
      os << gx << " " << gy << " " << gz << '\n';
    }
  }
}




void UCBspl::printVTKgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV, double scale) {

  cout << "Printing grid to vtk structured points file......" << endl;

  double umin = surf.umin();
  double vmin = surf.vmin();
  double umax = surf.umax();
  double vmax = surf.vmax();
  double du = (umax - umin)/(double)(noU-1);
  double dv = (vmax - vmin)/(double)(noV-1);  

  ofstream os(filename);

  os << "# vtk DataFile Version 2.0" << endl;
  os << "This file was generated by class UCButils (structured points)" << endl;
  os << "ASCII" << endl;
  os << "DATASET STRUCTURED_POINTS" << endl;
  os << "DIMENSIONS " << noU << " " << noV << " " << 1 << endl;
  os << "ORIGIN " << umin << " " << vmin << " " << 999 << endl;
  os << "SPACING " << du << " " << dv << " " << 999 << endl;
  os << "POINT_DATA " << noU*noV*1 << endl;
  os << "SCALARS volume_scalars float" << endl;
  os << "LOOKUP_TABLE default" << endl;

  // The points:
  //int noPoints = noU*noV;
  
  // size parameters: noX, noY, xmin, ymin, dx, dy
  //os << noU << " " << noV << " " << umin << " " << vmin << " " << du << " " << dv << endl;
  
  // z-arr row by row
  int i,j;
  for (j = 0; j < noV; j++) {
    double v = vmin + j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + i*du;      
      os << surf.f(u,v)*scale << '\n';
    }
  }
}


void UCBspl::printIRAPgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV) {
  
#ifdef MBA_DEBUG
  cout << "Printing grid to Irap grid file (format may be wrong): " << filename << endl;
#endif
  
  double du = (surf.umax() - surf.umin())/(double)(noU-1);
  double dv = (surf.vmax() - surf.vmin())/(double)(noV-1);  
  ofstream ofile(filename);
  ofile << noU << " " << noV << " " << du << " " << dv << endl;
  ofile << surf.umin() << " " << surf.umax() << " " << surf.vmin() << " " << surf.vmax() << endl;
  
  
  for (int j = 0; j < noV; j++) {
    double v = surf.vmin() + j*dv;
    for (int i = 0; i < noU; i++) {
      double u = surf.umin() + i*du;
      
      ofile << surf.f(u,v) << "\n";
    }
    //ofile << endl;
  }
}

void UCBspl::printGLgrid(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV) {
  
#ifdef MBA_DEBUG
  cout << "Printing GridStrips......" << endl;
#endif
  
  double umin = surf.umin();
  double vmin = surf.vmin();
  double umax = surf.umax();
  double vmax = surf.vmax();
  double du = (umax - umin)/(double)(noU-1);
  double dv = (vmax - vmin)/(double)(noV-1);  
  
  // Assumes uniform grid
  ofstream os(filename);
  
  // size parameters: noX, noY, xmin, ymin, dx, dy
  os << noU << " " << noV << " " << umin << " " << vmin << " " << du << " " << dv << endl;
  
  // z-arr row by row
  int i,j;
  for (j = 0; j < noV; j++) {
    double v = vmin + j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + i*du;      
      os << surf.f(u,v) << "\n";
    }
  }
  
  // normals row by row
  double gx,gy,gz;
  
  for (j = 0; j < noV; j++) {
    double v = vmin + (double)j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + (double)i*du;
      surf.normalVector(u,v,gx,gy,gz);
      os << gx << " " << gy << " " << gz << endl;
    }
  }
}


void UCBspl::printGLgridBin(const char filename[], const UCBspl::SplineSurface& surf, int noU, int noV,
                              const std::vector<double>& U, const std::vector<double>& V, const std::vector<double>& Z, // scattered data
                              double scale) {
  
#ifdef MBA_DEBUG
  cout << "Printing GridStrips (binary)......" << endl;
#endif
  

  double umin = surf.umin();
  double vmin = surf.vmin();
  double umax = surf.umax();
  double vmax = surf.vmax();
  double du = (umax - umin)/(double)(noU-1);
  double dv = (vmax - vmin)/(double)(noV-1);  
  
  // Assumes uniform grid
  FILE* fp = fopen(filename,"wb");
  
  fwrite(&noU,sizeof(int),1,fp);
  fwrite(&noV,sizeof(int),1,fp);
  fwrite(&umin,sizeof(double),1,fp);
  fwrite(&vmin,sizeof(double),1,fp);
  fwrite(&du,sizeof(double),1,fp);
  fwrite(&dv,sizeof(double),1,fp);
  
  // z-arr row by row
  int i,j;
  for (j = 0; j < noV; j++) {
    double v = vmin + j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + i*du;      
      double z = surf.f(u,v) * scale;
      fwrite(&z,sizeof(double),1,fp);
    }
  }
  
  // normals row by row
  double gx,gy,gz;
  
  for (j = 0; j < noV; j++) {
    double v = vmin + (double)j*dv;
    for (i = 0; i < noU; i++) {
      double u = umin + (double)i*du;
      surf.normalVector(u,v,gx,gy,gz);
      fwrite(&gx,sizeof(double),1,fp);
      fwrite(&gy,sizeof(double),1,fp);
      fwrite(&gz,sizeof(double),1,fp);
    }
  }
  
  // And the scattered data
  //const MBAdata data = mba.getData();
  int noPoints = U.size();
  fwrite(&noPoints,sizeof(int),1,fp);
  
#ifdef MBA_DEBUG
  cout << "No. points = " << noPoints << endl;
#endif

  for (i = 0; i < noPoints; i++) {
    double u = U[i];
    double v = V[i];
    double z = Z[i] * scale;
    
    fwrite(&u,sizeof(double),1,fp);
    fwrite(&v,sizeof(double),1,fp);
    fwrite(&z,sizeof(double),1,fp);
  }
  
  fclose(fp);
}


void UCBspl::saveSplineSurface(const char filename[], const UCBspl::SplineSurface& surf) {
  
#ifdef MBA_DEBUG
  cout << "Writing spline surface to ascii file: " << filename << endl;
#endif
  
  // store to file
  ofstream os(filename);
  
  // limits
  os << surf.umin() << '\n';
  os << surf.vmin() << '\n';
  os << surf.umax() << '\n';
  os << surf.vmax() << '\n';
    
  // continuity
  /*
  int cont = 1;
  if (!surf.C2_)
    cont = 0;
  os << cont << '\n';
  */
  
  
  // spline coefficients
  const GenMatrix<UCBspl_real>& PHI = *surf.getCoefficients();
  int noX = PHI.noX();
  int noY = PHI.noY(); 
  
  os << noX << " " << noY << '\n';
  for (int i = 0; i < noX; i++) {
    for (int j = 0; j < noY; j++) {
      os << PHI(i-1,j-1) << '\n';
    }
  }
}

void UCBspl::saveSplineSurfaceBin(const char filename[], const UCBspl::SplineSurface& surf) {
  
#ifdef MBA_DEBUG
  cout << "Writing spline surface to binary file: " << filename << endl;
#endif
  
  // store to binary file
  
  FILE* fp = fopen(filename,"wb");
  // limits
  
  double umin,vmin,umax,vmax;
  surf.getDomain(umin,vmin,umax,vmax);
  fwrite(&umin,sizeof(double),1,fp);
  fwrite(&vmin,sizeof(double),1,fp);
  fwrite(&umax,sizeof(double),1,fp);
  fwrite(&vmax,sizeof(double),1,fp);
  
  
  // continuity
  /*
  int cont = 1;
  if (!surf.C2_)
    cont = 0;
  fwrite(&cont,sizeof(int),1,fp);
  */
  
  // spline coefficients
  const GenMatrix<UCBspl_real>& PHI = *surf.getCoefficients();
  int noX = PHI.noX();
  int noY = PHI.noY(); 
  
  fwrite(&noX,sizeof(int),1,fp);
  fwrite(&noY,sizeof(int),1,fp);
  
  for (int i = 0; i < noX; i++) {
    for (int j = 0; j < noY; j++) {
      UCBspl_real coeff = PHI(i-1,j-1);
      fwrite(&coeff,sizeof(UCBspl_real),1,fp);
    }
  }
  
  fclose(fp);
}


void UCBspl::readSplineSurface(const char filename[], UCBspl::SplineSurface& surf) {
  
#ifdef MBA_DEBUG
  cout << "Reading spline surface surface from ascii file: " << filename << endl;
#endif
  
  ifstream is(filename);
  if (!is) {
    cout << "ERROR: File does not exist (or empty?)\n" << endl;
    exit(-1);
  }
  
  // limits
  double umin, vmin, umax, vmax;
  is >> umin;
  is >> vmin;
  is >> umax;
  is >> vmax;
  
  // offset
  double offset;
  is >> offset;
  
  // continuity
  int cont;
  is >> cont;
//   bool C2 = true;
//   if (cont == 0)
//     C2 = false;
  
  // spline coefficients
  int noX, noY;
  is >> noX;
  is >> noY;
#ifdef MBA_DEBUG
  cout << "Size of surface = " << noX << " X " << noY << endl;
#endif
  
  std::shared_ptr<GenMatrixType> PHI(new GenMatrix<UCBspl_real>(noX, noY));
  for (int i = 0; i < noX; i++) {
    for (int j = 0; j < noY; j++) {
      is >> (*PHI)(i-1,j-1);
    }
  }
  
  
  surf.init(PHI, umin, vmin, umax, vmax);
}


void UCBspl::readSplineSurfaceBin(const char filename[], UCBspl::SplineSurface& surf) {

#ifdef MBA_DEBUG
  cout << "Reading spline surface surface from binary file: " << filename << endl;
#endif
  
  FILE* fp = fopen(filename,"rb");
  
  if (fp == NULL) {
    cout << "ERROR: File does not exist (or empty?)\n" << endl;
    exit(-1);
  }
  
#ifdef MBA_DEBUG
  MBAclock rolex;
#endif
  // limits
  double umin, vmin, umax, vmax;
  fread(&umin,sizeof(double),1,fp);
  fread(&vmin,sizeof(double),1,fp);
  fread(&umax,sizeof(double),1,fp);
  fread(&vmax,sizeof(double),1,fp);
  
  // offset
  double offset;
  fread(&offset,sizeof(double),1,fp);
  
  // continuity
  int cont;
  fread(&cont,sizeof(int),1,fp);
  
//   bool C2 = true;
//   if (cont == 0)
//     C2 = false;
  
  // spline coefficients
  int noX, noY;
  fread(&noX,sizeof(int),1,fp);
  fread(&noY,sizeof(int),1,fp);
  
#ifdef MBA_DEBUG
  cout << "Size of surface = " << noX << " X " << noY << endl;
#endif
  
  std::shared_ptr<GenMatrixType> PHI(new GenMatrix<UCBspl_real>(noX, noY));
  UCBspl_real coef;
  for (int i = 0; i < noX; i++) {
    for (int j = 0; j < noY; j++) {
      fread(&coef,sizeof(UCBspl_real),1,fp);
      (*PHI)(i-1,j-1) = coef;
    }
  }
  fclose(fp);
#ifdef MBA_DEBUG
  cout << "Time used on reading data = " << rolex.getInterval() << endl;
#endif
  
  surf.init(PHI, umin, vmin, umax, vmax);
}

