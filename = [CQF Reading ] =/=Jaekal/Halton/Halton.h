//
// ======================================================================================
// Copyright © 2002-2008 Peter Jäckel. Distributed with "Monte Carlo methods in finance"
// published by John Wiley & Sons Ltd. All rights reserved.
// 
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// ======================================================================================
//

#include <cmath>
#include <vector>

#include "PrimeNumbers.h"

class Halton {
public:

  Halton(unsigned long dimensionality)
    : myDimensionality(dimensionality), mySequenceCounter(0), mySequenceVector(dimensionality) {}

  const std::vector<double>& nextUniformVector();

private:

  unsigned long myDimensionality, mySequenceCounter;

  std::vector<double> mySequenceVector;

  PrimeNumbers myPrimeNumbers;

};

