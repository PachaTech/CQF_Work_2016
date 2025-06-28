//
// ======================================================================================
// Copyright © 2002-2008 Peter Jäckel. Distributed with "Monte Carlo methods in finance"
// published by John Wiley & Sons Ltd. All rights reserved.
// 
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// ======================================================================================
//

#ifndef   BROWNIAN_BRIDGE_H
#define   BROWNIAN_BRIDGE_H

#include <vector>
#include <assert.h>

class BrownianBridge {    //    Builds Wiener process paths of equal time steps of delta t = 1.
public:
  BrownianBridge( unsigned long numberOfSteps );
  void buildPath( std::vector <double> &theWienerProcessPath, const std::vector <double> &gaussianVariates );
private:
  unsigned long numberOfSteps;
  std::vector <unsigned long> leftIndex, rightIndex, bridgeIndex;
  std::vector <double> leftWeight, rightWeight, sigma;
};

#endif // BROWNIAN_BRIDGE_H

