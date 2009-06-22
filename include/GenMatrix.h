//================================================================
//
// Created: September 18. 2000
//                                                                           
// Author: �yvind Hjelle <Oyvind.Hjelle@math.sintef.no>
//                                                                           
// Revised:
//                                                                           
// Description:
//                                                                           
//================================================================
// Copyright (c) 2000 SINTEF Applied Mathematics
//================================================================

#ifndef _GenMatrix_h_
#define _GenMatrix_h_

#include <iostream>

//#define ARRAY_CHECK 1

#ifdef ARRAY_CHECK
static void myMessage(int i, int j) {
  std::cout << "Matrix out of range: (" << i << "," << j << ")" << std::endl; 
  std::cout << " press [return]" << std::endl;
  getchar();
}
#endif

//===========================================================================
/** \brief A customized matrix class  
 *
 * GenMatrix - Matrix for storing B-spline coefficients in 
 * Multilevel B-spline approximation.
 * The indices goes from -1.
 * The class has reserve (and capasity) functionality similar to std::vector.
 * 
 * \author �yvind Hjelle <Oyvind.Hjelle@math.sintef.no>
 */
//===========================================================================

template <class Type>
class GenMatrix {
  Type** arr_;
  int noX_, noY_;
  int rX_, rY_; // reserved space
public:
  GenMatrix() {arr_ = NULL; noX_ = noY_ = rX_ = rY_ = 0;}
  GenMatrix(int noX, int  noY)
  {noX_ = noY_ = rX_ = rY_ = 0; arr_ = NULL; resize(noX, noY);}
  ~GenMatrix() {clear();}


  GenMatrix(const GenMatrix& G) { // copy constructor
    // By purpose, the copy constructor is not implemented.
    // Use the init function instead.
    std::cout << "GenMatrix:: copy constructor not implemented, use init()" << std::endl;
    std::cout << "EXIT FROM GenMatrix:: copy constructor" << std::endl;

    exit(-1);
  }
  
  
  void init(const GenMatrix& G) { // instead of copy constructor
    clear();
    reserve(G.rX_,G.rY_);
    resize(G.noX_,G.noY_);
    for (int j = 0; j < noY_; j++) {
      for (int i = 0; i < noX_; i++) {
        arr_[j][i] = G.arr_[j][i];
      }
    }
  }


  void swap(GenMatrix& mat) {
    std::swap(arr_, mat.arr_);

    std::swap(noX_, mat.noX_);
    std::swap(noY_, mat.noY_);
    std::swap(rX_ , mat.rX_);
    std::swap(rY_ , mat.rY_);
  }

  void operator += (GenMatrix& mat) {
    for (int j = 0; j < noY_; j++) {
      for (int i = 0; i < noX_; i++) {
        arr_[j][i] += mat.arr_[j][i];
      }
    }
  }

  Type norm_2() const
  {
      Type temp;
      Type res = 0;
      for (int j = 0; j < noY_; ++j) {
	  for (int i = 0; i < noX_; ++i) {
	      temp = arr_[j][i];
	      res += temp*temp;
	  }
      }
      return res;
  }


  // Temporary ?
  void operator += (double offset) {
    for (int j = 0; j < noY_; j++) {
      for (int i = 0; i < noX_; i++) {
        arr_[j][i] += offset;
      }
    }
  }
  
  // Note that resize does not preserve content (and no fill with zeros is done)
  void resize(int noX, int noY) {
    if (noX > rX_ || noY > rY_) {
      clear();
      rX_ = noX; rY_ = noY; // reserved = allocated
      arr_ = new Type*[rY_];
      for (int j = 0; j < rY_; j++)
        arr_[j] = new Type[rX_];
    }
    noX_ = noX; noY_ = noY;
  }
  
  void reserve(int rX, int rY) {resize(rX, rY); noX_ = 0; noY_ = 0;}

  void clear() {
    if (arr_) {
      for (int j = 0; j < rY_; j++)
        delete [] arr_[j];
      delete [] arr_;
      arr_ = NULL;
    }
    rX_ = rY_ = noX_ = noY_ = 0;
  }
  
  inline int noX() const {return noX_;}
  inline int noY() const {return noY_;}
  void capacity(int& rX, int& rY) const {rX = rX_; rY = rY_;}
  
  inline const Type& operator()(int ii, int jj) const {
    
#ifdef ARRAY_CHECK
    int i = ii+1;
    int j = jj+1;
    if (i < 0 || i >= noX_ || j < 0 || j >= noY_)
      myMessage(ii,jj);
#endif
    return arr_[jj+1][ii+1];
  }
  
  inline       Type& operator()(int ii, int jj) {
    
#ifdef ARRAY_CHECK
    int i = ii+1;
    int j = jj+1;
    if (i < 0 || i >= noX_ || j < 0 || j >= noY_)
      myMessage(ii,jj);
#endif
    return arr_[jj+1][ii+1];
  }
  
  
  void fill(Type val) {
    for (int j = 0; j < noY_; j++)
      for (int i = 0; i < noX_; i++)
        arr_[j][i] = val;
  }
  
  
  void print() const {
    std::cout << "Matrix..." << std::endl;
    for (int j = 0; j < noY_; j++) {
      for (int i = 0; i < noX_; i++) {
        std::cout << arr_[j][i] << " ";
      }
      std::cout << std::endl;
    }
  }
  
  
  void printBitmap() const {
    //cout << "GenMatrix::printBitmap... " << arr_[noY_-1][0] << endl;
    for (int j = 0; j < noY_; j++) {
      for (int i = 0; i < noX_; i++) {
        if (arr_[j][i] == 0)
          std::cout << 0;
        else
          std::cout << 1;
      }
      std::cout << std::endl;
    }  
  }
};

#endif
