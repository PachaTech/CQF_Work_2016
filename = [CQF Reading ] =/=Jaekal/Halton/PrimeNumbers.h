//
// ======================================================================================
// Copyright © 2002-2008 Peter Jäckel. Distributed with "Monte Carlo methods in finance"
// published by John Wiley & Sons Ltd. All rights reserved.
// 
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// ======================================================================================
//

#include <vector>

class PrimeNumbers {
public:

  PrimeNumbers(){}

  // Get and store one after another.
  unsigned long nextPrimeNumber();

  unsigned long indexedPrimeNumber(unsigned long absoluteIndex /* starting at 0 in C convention */){
    while (primeNumbers.size()<=absoluteIndex)
      nextPrimeNumber();
    return primeNumbers[absoluteIndex];
  }

  unsigned long operator[](unsigned long absoluteIndex /* starting at 0 in C convention */){
    return indexedPrimeNumber(absoluteIndex);
  }

private:
  std::vector<unsigned long> primeNumbers;
};

