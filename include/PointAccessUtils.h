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

