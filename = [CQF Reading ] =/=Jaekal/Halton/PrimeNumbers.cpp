//
// ======================================================================================
// Copyright © 2002-2008 Peter Jäckel. Distributed with "Monte Carlo methods in finance"
// published by John Wiley & Sons Ltd. All rights reserved.
// 
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// ======================================================================================
//

#include "PrimeNumbers.h"
#include <cmath>

unsigned long PrimeNumbers::nextPrimeNumber() {
  unsigned long m;
  if (primeNumbers.size()) {
    unsigned long i, n, p;
    m = primeNumbers.back();
    do {
      i = 0;
      ++m;
      n = static_cast<unsigned long>(sqrt(m));
      do {
        p = primeNumbers[i];
        ++i;
      }
      while ( m%p && p<=n );
    } while ( p<=n );
  } else {
    m=2;
  }
  primeNumbers.push_back(m);
  return m;
}

