//
// ======================================================================================
// Copyright © 2002-2008 Peter Jäckel. Distributed with "Monte Carlo methods in finance"
// published by John Wiley & Sons Ltd. All rights reserved.
// 
// Permission to use, copy, modify, and distribute this software is freely granted,
// provided that this notice is preserved.
// ======================================================================================
//

#include "BrownianBridge.h"
#include <cmath>

BrownianBridge::BrownianBridge( unsigned long numberOfSteps_ ) :
  numberOfSteps(numberOfSteps_),
  leftIndex(numberOfSteps), rightIndex(numberOfSteps), bridgeIndex(numberOfSteps),
  leftWeight(numberOfSteps), rightWeight(numberOfSteps), sigma(numberOfSteps)
{
  assert(numberOfSteps);               //  There must be at least one step.
  std::vector <unsigned long> map(numberOfSteps);
  // map is used to indicate which points are already constructed. If map[i] is zero, path point i
  // is yet unconstructed.  map[i]-1 is the index of the variate that constructs the path point # i.
  unsigned long i,j,k,l;
  map[numberOfSteps-1] = 1;            //  The first point in the construction is the global step.
  bridgeIndex[0] = numberOfSteps-1;    //  The global step is constructed from the first variate.
  sigma[0] = sqrt(numberOfSteps);      //  The variance of the global step is numberOfSteps*1.0.
  leftWeight[0] = rightWeight[0] = 0.; //  The global step to the last point in time is special.
  for (j=0,i=1;i<numberOfSteps;++i){
    while (map[j]) ++j;                //  Find the next unpopulated entry in the map.
    k=j;
    while ((!map[k])) ++k;             //  Find the next populated entry in the map from there.
    l=j+((k-1-j)>>1);                  //  l-1 is now the index of the point to be constructed next.
    map[l]=i;
    bridgeIndex[i] = l;                //  The i-th Gaussian variate will be used to set point l-1.
    leftIndex[i]   = j;
    rightIndex[i]  = k;
    leftWeight[i]  = (k-l)/(k+1.-j);
    rightWeight[i] = (l+1.-j)/(k+1.-j);
    sigma[i]       = sqrt(((l+1.-j)*(k-l))/(k+1.-j));
    j=k+1;
    if (j>=numberOfSteps) j=0;	//	Wrap around.
  }
}

void BrownianBridge::buildPath( std::vector <double> &path, const std::vector <double> &gaussianVariates ){
  assert( gaussianVariates.size() == numberOfSteps && path.size() == numberOfSteps );
  unsigned long i,j,k,l;
  path[numberOfSteps-1] = sigma[0]*gaussianVariates[0];       //  The global step.
  for (i=1;i<numberOfSteps;i++){
    j = leftIndex[i];
    k = rightIndex[i];
    l = bridgeIndex[i];
    if (j) path[l] = leftWeight[i]*path[j-1] + rightWeight[i]*path[k] + sigma[i]*gaussianVariates[i];
    else   path[l] = rightWeight[i]*path[k] + sigma[i]*gaussianVariates[i];
  }
}

