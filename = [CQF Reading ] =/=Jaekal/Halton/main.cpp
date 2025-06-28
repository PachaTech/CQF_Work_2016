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
#include <iostream>

//
//
// Please note that I do *not* advocate the use of Halton numbers in finanical applications.
// The code here is only provided for comparative purposes.
//
// Instead of Halton numbers, I recommend the use of Sobol' numbers as they can be constructed
// using the algorithms published in "Numerical Recipes in C", W. H. Press, S. A. Teukolsky,
// W. t. Vetterling, B. P. Flannery, second edition, section 7.7, page 309, ISBN 0-521-43108-5.
// The dimensionality of the Sobol' sequence can be extended to virtually any size you ever might
// need by the aid of the table of primitive polynomials modulo two as provided ready for
// compilation in the directory "PrimitivePolynomialsModuloTwo" on the CD accompanying the book
// "Monte Carlo Methods in Finance" by Peter Jäckel. You should initialise them as described in
// the section "The initialisation of Sobol' numbers". Note that you should always use Sobol'
// sequences in conjunction with the Brownian Bridge.
//
//

int main(){
  unsigned long i, j, d=6, n=17;
  Halton theHaltonSequence(d);
  for (i=0;i<n;++i){
    const std::vector<double>& uniformVector( theHaltonSequence.nextUniformVector() );
    std::cout << i << " : ";
    for (j=0;j<d;++j)
      std::cout << " " << uniformVector[j];
    std::cout << std::endl;
  }
  return 0;
}

// This program should produce the following output:
//
//  0 :  0.5 0.333333 0.2 0.142857 0.0909091 0.0769231
//  1 :  0.25 0.666667 0.4 0.285714 0.181818 0.153846
//  2 :  0.75 0.111111 0.6 0.428571 0.272727 0.230769
//  3 :  0.125 0.444444 0.8 0.571429 0.363636 0.307692
//  4 :  0.625 0.777778 0.04 0.714286 0.454545 0.384615
//  5 :  0.375 0.222222 0.24 0.857143 0.545455 0.461538
//  6 :  0.875 0.555556 0.44 0.0204082 0.636364 0.538462
//  7 :  0.0625 0.888889 0.64 0.163265 0.727273 0.615385
//  8 :  0.5625 0.037037 0.84 0.306122 0.818182 0.692308
//  9 :  0.3125 0.37037 0.08 0.44898 0.909091 0.769231
//  10 :  0.8125 0.703704 0.28 0.591837 0.00826446 0.846154
//  11 :  0.1875 0.148148 0.48 0.734694 0.0991736 0.923077
//  12 :  0.6875 0.481481 0.68 0.877551 0.190083 0.00591716
//  13 :  0.4375 0.814815 0.88 0.0408163 0.280992 0.0828402
//  14 :  0.9375 0.259259 0.12 0.183673 0.371901 0.159763
//  15 :  0.03125 0.592593 0.32 0.326531 0.46281 0.236686
//  16 :  0.53125 0.925926 0.52 0.469388 0.553719 0.313609

