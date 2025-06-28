//
// ======================================================================================
// Copyright © 2002-2008 Peter Jäckel. Distributed with "Monte Carlo methods in finance"
// published by John Wiley & Sons Ltd. All rights reserved.
// 
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// ======================================================================================
//

#include "Halton.h"

const std::vector<double>& Halton::nextUniformVector() {
  ++mySequenceCounter;
  unsigned long b,i,k;
  double f, h;
  for (i=0;i<myDimensionality;++i){
    h = 0;
    b = myPrimeNumbers[i];
    f = 1.;
    k = mySequenceCounter;
    while (k){
      f /= b;
      h += (k%b)*f;
      k /= b;
    }
    mySequenceVector[i] = h;
  }
  return mySequenceVector;
}

