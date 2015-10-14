#ifndef __TRAIN_H__
#define __TRAIN_H__

#include <Rcpp.h>
#include "bwpmf.h"

struct Phi {
  
  double *data;
  
  Phi() : data(new double[Param::K]) 
  { }
  
  ~Phi() { delete [] data; }
  
};

typedef ListOfList<Phi> PhiList;

#endif // __TRAIN_H__