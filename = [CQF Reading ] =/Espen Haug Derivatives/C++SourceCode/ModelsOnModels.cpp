/*
 *  ModelsOnModels.cpp
 *  OpsjonerC++
 *
 *  Created by Espen Gaarder Haug on 10/5/06.
 *  Copyright 2007 Espen Gaarder Haug. All rights reserved.
 *
 */

#include "ModelsOnModels.h"


//CRR TREE USED IN NEGATIVE PROBABILITY CHAPTER

// Cox-Ross-Rubinstein binomial tree
double CRRBinomial(char OutputFLag, char AmeEurFlag, char CallPutFlag, double S, double X, double T,
                double r, double b, double v, int n)
{

// By Espen Gaarder Haug

	std::vector<double> OptionValue(n + 1);
    double ReturnValue[4];
    double u, d, p, dt, Df;
    int i, j, z;

    if (CallPutFlag == 'c')
    {
    z = 1;
    }
    else if (CallPutFlag == 'p' )
    {
    z = -1;
    };

    dt = T / n;
    u = exp(v * sqrt(dt));
    d = 1.0 / u;
    p = (exp(b * dt) - d) / (u - d);
    Df = exp(-r * dt);

    for (i = 0; i<= n; i ++)
    {
         OptionValue[i] = esmax(0, z * (S * pow(u, i) * pow(d, n - i) - X));
    };

    for (j = n - 1 ; j>= 0 ; j--)
    {
        for (i = 0; i<= j; i ++)
        {
                OptionValue[i] = (p * OptionValue[i + 1] + (1 - p) * OptionValue[i]) * Df;
            if (AmeEurFlag == 'a' )
             {
              OptionValue[i] = esmax(z * (S * pow(u, i) * pow(d,abs(i - j)) - X), OptionValue[i]);
            };
        };
     if (j == 2)
     {
           ReturnValue[2] = ((OptionValue[2] - OptionValue[1]) / (S * u * u - S)
                - (OptionValue[1] - OptionValue[0]) / (S - S * d * d)) / (0.5 * (S * u * u - S * d * d));
                ReturnValue[3] = OptionValue[1];
      }
      else if (j == 1)
      {
           ReturnValue[1] = (OptionValue[1] - OptionValue[0]) / (S * u - S * d);
      }
    };

   	ReturnValue[3] = (ReturnValue[3] - OptionValue[0]) / (2.0 * dt) / 365.0;
   	ReturnValue[0] = OptionValue[0];
    if (OutputFLag == 'p')  //Option value
        return ReturnValue[0];
    else if (OutputFLag == 'd')  //Delta
        return ReturnValue[1];
    else if (OutputFLag == 'g')  //Gamma
        return ReturnValue[2];
    else if (OutputFLag == 't')  //Theta
        return ReturnValue[3];

    return 0;
};


//ENERGY OPTIONS USED IN CHAPTER POWER DERIVATIVES

// Energy Swaption
double EnergySwaption(char CallPutFlag, double F, double X, double T, double Tb, double rj, double rb, int j, int n, double v)
{

// By Espen Gaarder Haug

    double d1, d2, Df, OptionValue;
	
    // T: years to option expiry
    // Tb: years to start of swap delivery period (Tb>=T)
    // rb: zero coupon rate from now to start of swap delivery period
    // rj: swap rate covering delivery period with j compoundings per year
    // j: number of compoundings per year
    // n: number of days in delivery period
    
     Df = (1.0 - 1.0 / pow(1.0 + rj / j, n) ) / rj * j / n;
  
     d1 = (log(F / X) + v * v / 2.0 * T) / (v * sqrt(T));
     d2 = d1 - v * sqrt(T);

    if (CallPutFlag == 'c' )
	{
        OptionValue = Df * exp(-rb * Tb) * (F * CND(d1) - X * CND(d2));
	}
    else if (CallPutFlag == 'p' )
	{
        OptionValue = Df * exp(-rb * Tb) * (X * CND(-d2) - F * CND(-d1));
	}
    return OptionValue;
    
};

// Energy swaption approximation
double EnergySwaptionApproximation(char CallPutFlag, double F, double X , double T, double Tm, double re, double v) 
{

// By Espen Gaarder Haug

    double d1, d2, OptionValue;
    // T: years to option expiry
    // Tm: years to start of swap delivery period (Tm>=T)
    // re: zero coupon rate from now to start of swap delivery period
    
    d1 = (log(F / X) + v * v / 2.0 * T) / (v * sqrt(T));
	    d2 = d1 - v * sqrt(T);

    if (CallPutFlag == 'c' )
	{
        OptionValue = exp(-re * Tm) * (F * CND(d1) - X * CND(d2));
	}
    else if ( CallPutFlag == 'p' )
	{
        OptionValue = exp(-re * Tm) * (X * CND(-d2) - F * CND(-d1));
    }
	
	return OptionValue;
    
};



// This is the generlaized Black-Scholes-Merton formula including all greeeks
// This function is simply calling all the other functions
double inCGBlackScholes(char OutputFlag[], char CallPutFlag, double S, double X,double T, double r, double b, double v)
{

// By Espen Gaarder Haug

                    double output;
                  
                    output = 0.0;
                    
                if (strcmp(OutputFlag,"p") == 0) //Value
				{
                    output = GBlackScholes(CallPutFlag, S, X, T, r, b, v);
				}
                // DELTA GREEKS
                 else if (strcmp(OutputFlag,"d") == 0) // Delta
				{
                    output = GDelta(CallPutFlag, S, X, T, r, b, v);
				}
				else if (strcmp(OutputFlag,"dddv") == 0) //DDeltaDvol
                {
					    output = GDdeltaDvol(S, X, T, r, b, v) / 100.0;
				}
                else if (strcmp(OutputFlag,"dt") == 0) //DDeltaDtime/Charm
				{
                    output =  GDdeltaDtime(CallPutFlag, S, X, T, r, b, v) / 365.0;
				}
                else if (strcmp(OutputFlag,"dmx") == 0)
				{
                    output = S * S / X * exp((2.0 * b + v * v) * T);
				}
                 else if (strcmp(OutputFlag,"e") == 0) // Elasticity
				{
                     output = GElasticity(CallPutFlag, S, X, T, r, b, v);
				}
                // GAMMA GREEKS
                else if (strcmp(OutputFlag,"sg") == 0) //SaddleGamma
                {
					output = GSaddleGamma(X, T, r, b, v);
				}
				else if (strcmp(OutputFlag,"g") == 0) //Gamma
                {
					    output = GGamma(S, X, T, r, b, v);
				}
				else if (strcmp(OutputFlag,"s") == 0) //DgammaDspot/speed
				{
					   output = GDgammaDspot(S, X, T, r, b, v);
				}
				else if (strcmp(OutputFlag,"gv") == 0) //DgammaDvol/Zomma
                {
					output = GDgammaDvol(S, X, T, r, b, v) / 100.0;
				}
				else if (strcmp(OutputFlag,"gt") == 0) //DgammaDtime
                {    
					output = GDgammaDtime(S, X, T, r, b, v) / 365.0;
                }
				else if (strcmp(OutputFlag,"gp") == 0) //GammaP
                {
					    output = GGammaP(S, X, T, r, b, v);
                }
				else if (strcmp(OutputFlag,"gps") == 0) //DgammaPDspot
                {
					    output = GDgammaPDspot(S, X, T, r, b, v);
                }
				else if (strcmp(OutputFlag,"gpv") == 0) //DgammaDvol/Zomma
                {
					    output = GDgammaPDvol(S, X, T, r, b, v) / 100.0;
                }
				else if (strcmp(OutputFlag,"gpt") == 0) //DgammaPDtime
                {
					output = GDgammaPDtime(S, X, T, r, b, v) / 365.0;
                }
				// VEGA GREEKS
				else if (strcmp(OutputFlag,"v") == 0) //Vega
                {
					output = GVega(S, X, T, r, b, v) / 100.0;
				}
				else if (strcmp(OutputFlag,"vt") == 0) //DvegaDtime
                {
					    output = GDvegaDtime(S, X, T, r, b, v) / 36500.0;
                }
				else if (strcmp(OutputFlag,"dvdv") == 0) //DvegaDvol/Vomma
                {
					    output = GDvegaDvol(S, X, T, r, b, v) / 10000.0;
                }
				else if (strcmp(OutputFlag,"vp") == 0) //VegaP
                {
					    output = GVegaP(S, X, T, r, b, v);
                }
				else if (strcmp(OutputFlag,"vpv") == 0) //DvegaPDvol/VommaP
                {
					    output = GDvegaPDvol(S, X, T, r, b, v) / 100.0;
                }
				else if (strcmp(OutputFlag,"vl") == 0) //Vega Leverage
                {
					    output = GVegaLeverage(CallPutFlag, S, X, T, r, b, v);
                }
				// THETA GREEKS
                else if (strcmp(OutputFlag,"t") == 0) //Theta
                {
					    output = GTheta(CallPutFlag, S, X, T, r, b, v) / 365.0;
                }
				else if (strcmp(OutputFlag,"Dlt") == 0) //Drift-less Theta
                {
					    output = GThetaDriftLess(S, X, T, r, b, v) / 365.0;
                }
				// RATE/CARRY GREEKS
                else if (strcmp(OutputFlag,"r") == 0) //Rho
                {
					    output = GRho(CallPutFlag, S, X, T, r, b, v) / 100.0;
				}
				else if (strcmp(OutputFlag,"fr") == 0) //Rho futures option
                {
					    output = GRhoFO(CallPutFlag, S, X, T, r, b, v) / 100.0;
                }
				else if (strcmp(OutputFlag,"b") == 0) //Carry Rho
                {
					    output = GCarry(CallPutFlag, S, X, T, r, b, v) / 100.0;
                }
				else if (strcmp(OutputFlag,"f") == 0) //Phi/Rho2
                {
					    output = GPhi(CallPutFlag, S, X, T, r, b, v) / 100.0;
                }
				// PROB GREEKS
                else if (strcmp(OutputFlag,"z") == 0) //Zeta/In-the-money risk neutral probability
                {
					    output = GInTheMoneyProbability(CallPutFlag, S, X, T, b, v);
                }
				else if (strcmp(OutputFlag,"zv") == 0) //DzetaDvol
                {
					    output = GDzetaDvol(CallPutFlag, S, X, T, r, b, v) / 100.0;
                }
				else if (strcmp(OutputFlag,"zt") == 0) //DzetaDtime
                {
					    output = GDzetaDtime(CallPutFlag, S, X, T, r, b, v) / 365.0;
                }
				else if (strcmp(OutputFlag,"bp") == 0) // Break even probability
                {
					    output = GBreakEvenProbability(CallPutFlag, S, X, T, r, b, v);
				}
				else if (strcmp(OutputFlag,"dx") == 0) //StrikeDelta
                {
					    output = GStrikeDelta(CallPutFlag, S, X, T, r, b, v);
                }
				else if (strcmp(OutputFlag,"dxdx") == 0)//Risk Neutral Density
                {
					    output = GRiskNeutralDensity(S, X, T, r, b, v);
                 }   
				//CALCULATIONS
                else if (strcmp(OutputFlag,"d1") == 0) //d1
                {
					    output = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
                }
				else if (strcmp(OutputFlag,"d2") == 0) // d2
                {
					    output = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));
                }
				else if (strcmp(OutputFlag,"nd1") == 0) // n(d1)
                {
					    output = ND((log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T)));
                }
				else if (strcmp(OutputFlag,"nd2") == 0) // n(d2)
                {
					    output = ND((log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T)));
                }
				else if (strcmp(OutputFlag,"CNDd1") == 0)  // N(d1)
                {
					    output = CND((log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T)));
                }
				else if (strcmp(OutputFlag,"CNDd2") == 0) //N(d2)
				{
                    output = CND((log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T)));
               }
			   
			   return output;
};


// Converting interest rate to continuous compounding rate
double ConvertingToCCRate( double r , double Compoundings) 
{
    if (Compoundings == 0) 
        return  r;
    else
        return  Compoundings * log(1.0 + r / Compoundings) ;  
}

// sgn function (PRIVATE FUNCTION)
double sgn( double x )
{
 if(x>=0.0) return 1.0;
 return -1.0;
}

// Generalized Black-Scholes Option Pricing Formula
double GBlackScholes(char CallPutFlag, double S, double X, double T, double r, double b, double v)
{
	 // By Espen Gaarder Haug
	 // CallPutFlag = 'c' for call and 'p' for put option, Default 'c'
    
 double d1, d2;

 d1=(log(S/X)+(b+v*v/2.0)*T)/(v*sqrt(T));
 d2=d1-v*sqrt(T);

 if(CallPutFlag == 'c')
  return  S*exp((b-r)*T)*CND(d1)-X*exp(-r*T)*CND(d2);
 else if(CallPutFlag == 'p')
  return  X*exp(-r*T)*CND(-d2)-S*exp((b-r)*T)*CND(-d1);
};

// *************************************************************************************
//										DELTA GREEKS
// *************************************************************************************

// Delta for the generalized Black-Scholes formula (PRIVATE FUNCTION)
double GDelta(double CallPutFlag, double S,double X,double T ,double r, double b, double v)
// By Espen Gaarder Haug
{
    double d1;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));

    if (CallPutFlag == 'c') {
        return exp((b - r) * T) * CND(d1);
    }
    else if (CallPutFlag == 'p'){
        return exp((b - r) * T) * (CND(d1) - 1.0);
    }
};

// DDeltaDvol also known as vanna
double  GDdeltaDvol(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug

	double d1, d2;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
   return -exp((b - r) * T) * d2 / v * ND(d1);
}

// DdeltaDtime/Charm for the generalized Black and Scholes formula
double  GDdeltaDtime(char CallPutFlag , double  S, double X , double T, double r, double b, double v) 
{
// By Espen Gaarder Haug
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
    
    if (CallPutFlag == 'c' )
	{
          return  -exp((b - r) * T) * (ND(d1) * (b / (v * sqrt(T)) - d2 / (2.0 * T)) + (b - r) * CND(d1));
    }
	else if ( CallPutFlag == 'p' )
	{
        return -exp((b - r) * T) * (ND(d1) * (b / (v * sqrt(T)) - d2 / (2.0 * T)) - (b - r) * CND(-d1));
	}
   
    
};

// Elasticity for the generalized Black and Scholes formula
double GElasticity(char CallPutFlag, double S, double X, double T, double r, double b, double v)
 {   
 // By Espen Gaarder Haug
        return GDelta(CallPutFlag, S, X, T, r, b, v) * S / GBlackScholes(CallPutFlag, S, X, T, r, b, v);
    
}


// Closed form solution to find strike given the delta
double GStrikeFromDelta(char CallPutFlag, double S, double T, double r, double b, double v, double delta)
{
// By Espen Gaarder Haug
        
       if (CallPutFlag == 'c' )
	   {
          return S * exp(-CNDEV(delta * exp((r - b) * T)) * v * sqrt(T) + (b + v * v / 2.0) * T);
		}
        {
            return S * exp(CNDEV(-delta * exp((r - b) * T)) * v * sqrt(T) + (b + v * v / 2.0) * T);
        }
        
};

// *************************************************************************************
//										GAMMA GREEKS
// *************************************************************************************


// Gamma for the generalized Black-Scholes formula (PRIVATE FUNCTION)
double GGamma(double S,double X,double T ,double r, double b, double v)
{
// By Espen Gaarder Haug
    double d1;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));

    return ND(d1)*exp((b - r) * T) /(S*v*sqrt(T));

};

// GammaP for the generalized Black and Scholes formula
double GGammaP(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    return  S * GGamma(S, X, T, r, b, v) / 100.0;
    
};


// DgammaDvol/Zomma for the generalized Black and Scholes formula
double GDgammaDvol(double S, double X, double T, double r, double b, double v) 
{
// By Espen Gaarder Haug
    
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
	
   return GGamma(S, X, T, r, b, v) * ((d1 * d2 - 1.0) / v);

};

// DgammaPDvol for the generalized Black and Scholes formula
double GDgammaPDvol(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
    
	 return S / 100.0 * GGamma(S, X, T, r, b, v) * ((d1 * d2 - 1.0) / v);

};

// DgammaDspot/Speed for the generalized Black and Scholes formula
double GDgammaDspot(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double d1;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    
   return -GGamma(S, X, T, r, b, v) * (1.0 + d1 / (v * sqrt(T))) / S;

};


// DgammaPDspot/SpeedP for the generalized Black and Scholes formula
double GDgammaPDspot(double S, double X ,double T ,double r,double b,double v )
{
// By Espen Gaarder Haug
    
    double d1 ;
    
    d1 = (log(S / X) + (b + v * v/ 2.0) * T) / (v * sqrt(T));
    
    return -GGamma(S, X, T, r, b, v) * (d1) / (100.0 * v * sqrt(T));

};


// GGammaDtime for the generalized Black and Scholes formula
double GDgammaDtime(double S, double X, double T, double r, double b, double v) 
{
// By Espen Gaarder Haug
    
    double d1 , d2 ;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
	
	return GGamma(S, X, T, r, b, v) * (r - b + b * d1 / (v * sqrt(T)) + (1.0 - d1 * d2) / (2.0 * T));

};

// GGammaPDtime for the generalized Black and Scholes formula
double GDgammaPDtime(double S, double X, double T, double r, double b, double v) 
{
// By Espen Gaarder Haug
    
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
    
	return GGammaP(S, X, T, r, b, v) * (r - b + b * d1 / (v * sqrt(T)) + (1.0 - d1 * d2) / (2.0 * T));

};


// SaddleGamma for the generalized Black and Scholes formula
double  GSaddleGamma(double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    return sqrt(exp(1.0) / Pi) * sqrt((2.0 * b - r) / (v * v) + 1.0) / X;
    
};




// *************************************************************************************
//										VEGA GREEKS
// *************************************************************************************



// Vega for the generalized Black-Scholes formula (PRIVATE FUNCTION)
double GVega(double S,double X,double T ,double r,
                double b, double v)
{
// By Espen Gaarder Haug

    double d1;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));

    return S* exp((b - r) * T) * ND(d1) *sqrt(T);

};

// VegaP for the generalized Black and Scholes formula
double GVegaP(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug

    return v / 10.0 * GVega(S, X, T, r, b, v);

};


// Vega for the generalized Black and Scholes formula
double GVegaLeverage(char CallPutFlag, double S, double X , double T , double r , double b , double v ) 
{
// By Espen Gaarder Haug
    
    return GVega(S, X, T, r, b, v) * v / GBlackScholes(CallPutFlag, S, X, T, r, b, v);

};

// DvegaDvol/Vomma for the generalized Black and Scholes formula
double GDvegaDvol(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
    
	return GVega(S, X, T, r, b, v) * d1 * d2 / v;

};


// DvegaPDvol/VommaP for the generalized Black and Scholes formula
double GDvegaPDvol(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
    
	return GVegaP(S, X, T, r, b, v) * d1 * d2 / v;

};


// DvegaDtime for the generalized Black and Scholes formula
double GDvegaDtime(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double  d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
	
    return GVega(S, X, T, r, b, v) * (r - b + b * d1 / (v * sqrt(T)) - (1.0 + d1 * d2) / (2.0 * T));

};


// *************************************************************************************
//										THETA GREEKS
// *************************************************************************************


// Theta for the generalized Black-Scholes formula (PRIVATE FUNCTION)
double GTheta(double CallPutFlag, double S,double X,double T ,double r, double b, double v)
{
// By Espen Gaarder Haug

    double d1, d2;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1- v * sqrt(T);

    if (CallPutFlag == 'c') {
        return -S*exp((b - r) * T) * ND(d1)*v/(2.0 *sqrt(T))
        -(b-r)*S*exp((b - r) * T)*CND(d1)-r*X*exp(-r*T)*CND(d2);
    }
    else if (CallPutFlag == 'p')
	{
        return -S*exp((b - r) * T) * ND(d1)*v/(2.0 *sqrt(T))
        +(b-r)*S*exp((b - r) * T)*CND(-d1)+r*X*exp(-r*T)*CND(-d2);
    }
};


// Drift-less Theta for the generalized Black and Scholes formula
double GThetaDriftLess(double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double d1;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    
	return -S * exp((b - r) * T) * ND(d1) * v / (2.0 * sqrt(T));
    
};

// *************************************************************************************
//										RHO GREEKS
// *************************************************************************************


// Rho for the generalized Black-Scholes formula 
// Do not work for options on futures
double GRho(double CallPutFlag, double S,double X,double T ,double r,
                double b, double v)
{
// By Espen Gaarder Haug

    double d2;

    d2 = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));

    if (CallPutFlag == 'c')
    {
          return T*X*exp(-r*T)*CND(d2);
    }
    else if (CallPutFlag == 'p')
    {
          return -T*X*exp(-r*T)*CND(-d2);
	}

};

// Rho for the generalized Black and Scholes formula for Futures option
double GRhoFO(char CallPutFlag, double S, double X, double T, double r, double b, double v) 
{
// By Espen Gaarder Haug
    
            return  -T * GBlackScholes(CallPutFlag, S, X, T, r, 0.0, v);
   
};

// Phi/Rho2 for the generalized Black-Scholes formula  
double GPhi(double CallPutFlag, double S,double X,double T ,double r, double b, double v)
{
// By Espen Gaarder Haug

    double d1;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));

    if (CallPutFlag == 'c')
    {
        return -T*S*exp((b-r)*T)*CND(d1);
    }
    else if (CallPutFlag == 'p'){
        return T*S*exp((b-r)*T)*CND(-d1);
    }

};

// Carry for the generalized Black-Scholes formula  (PRIVATE FUNCTION)
double GCarry(double CallPutFlag, double S,double X,double T ,double r,
                double b, double v)
{
// By Espen Gaarder Haug

    double d1;

    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));

    if (CallPutFlag == 'c') {
        return T*S*exp((b-r)*T)*CND(d1);
    }
    else if (CallPutFlag == 'p'){
        return -T*S*exp((b-r)*T)*CND(-d1);
    }

};


// *************************************************************************************
//										PROBABILITY GREEKS
// *************************************************************************************




// Risk neutral probability of ending up in-the-money for the generalized Black and Scholes formula
double GInTheMoneyProbability(char CallPutFlag, double S, double X, double T, double b, double v) 
{
// By Espen Gaarder Haug
                
    double d2;
    
    d2 = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));
    
    if ( CallPutFlag == 'c' )
	{
        return CND(d2);
	}
    else if (CallPutFlag == 'p' )
	{
        return CND(-d2);
	};
    
};

// StrikeDelta for the generalized Black and Scholes formula
double GStrikeDelta(char CallPutFlag, double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    
    double d2;
    
    d2 = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));
    if (CallPutFlag == 'c' )
	{
        return -exp(-r * T) * CND(d2);
	}
	else
    {
        return exp(-r * T) * CND(-d2);
    }
    
};


// DZetaDvol for the generalized Black and Scholes formula
double GDzetaDvol(char CallPutFlag, double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
                
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
    if (CallPutFlag == 'c') 
	{
        return -ND(d2) * d1 / v;
	}
    else
	{
        return ND(d2) * d1 / v;
	}
	
};

// DZetaDtime for the generalized Black and Scholes formula
double GDzetaDtime(char CallPutFlag, double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
                
    double d1, d2;
    
    d1 = (log(S / X) + (b + v * v / 2.0) * T) / (v * sqrt(T));
    d2 = d1 - v * sqrt(T);
	
    if (CallPutFlag == 'c' )
	{
        return ND(d2) * (b / (v * sqrt(T)) - d1 / (2.0 * T));
	}
    else
	{
       return -ND(d2) * (b / (v * sqrt(T)) - d1 / (2.0 * T));
     }

};

// Risk neutral break even probability for the generalized Black and Scholes formula
double GBreakEvenProbability(char CallPutFlag, double S, double X, double T, double r, double b , double v)
{
// By Espen Gaarder Haug
                
    double d2;
    
    if (CallPutFlag = 'c') 
	{
        X = X + GBlackScholes('c', S, X, T, r, b, v) * exp(r * T);
        d2 = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));
        return CND(d2);
	}
    else if (CallPutFlag == 'p' )
	{
        X = X - GBlackScholes('p', S, X, T, r, b, v) * exp(r * T);
        d2 = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));
        return  CND(-d2);
    }
    
};

// Risk Neutral Denisty for the generalized Black and Scholes formula
double GRiskNeutralDensity(double S,double X,double T ,double r, double b, double v)
{
// By Espen Gaarder Haug

    double d2;

    d2 = (log(S / X) + (b - v * v / 2.0) * T) / (v * sqrt(T));

    return ND(d2)*exp(-r * T) /(X*v*sqrt(T));

};


// *************************************************************************************
//										EXOTIC OPTIONS
// *************************************************************************************


// American Knock-in Barrier options, call down-and-in and put up-and-in
// Here using the Bjerksund-Stensland 2002 approximation
double AmericanKnockInBarriers(char CallPutFlag, double S, double X, 
        double H, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
		
		double OptionValue;
		
        if (CallPutFlag == 'c' ) // call down-and-in
		{
            if (H <= X )
			{
                OptionValue = pow(S / H, 1.0 - 2.0 * b / (v*v)) 
                * BSAmericanApprox2002('c', H * H / S, X, T, r, b, v);
			}
           else if ( H <= esmax(X, r / (r - b) * X) )
		   {
                OptionValue = pow(S / H, 1.0 - 2.0 * b / (v * v)) 
                * (BSAmericanApprox2002('c', H * H / S, X, T, r, b, v) 
                    - GBlackScholes('c', H * H / S, X, T, r, b, v)) 
                    + StandardBarrier("cdi", S, X, H, 0, T, r, b, v);
           }
		}
         else if (CallPutFlag == 'p' )
		 {
            if (H >= X )
			{
                OptionValue = pow(S / H, 1.0 - 2.0 * b /  (v * v) ) 
                * BSAmericanApprox2002('p', H * H / S, X, T, r, b, v);
			}
           else if( H >= esmax(X, r / (r - b) * X) )
		   {
                OptionValue = pow(S / H, 1.0 - 2.0 * b / (v * v)) 
                * (BSAmericanApprox2002(CallPutFlag, H * H / S, X, T, r, b, v) 
                    - GBlackScholes('p', H * H / S, X, T, r, b, v)) 
                    + StandardBarrier("pui", S, X, H, 0, T, r, b, v);
          }
        
       }
	   return OptionValue;
};




// The Bjerksund and Stensland (2002) American approximation
// This approximaion option function is used by American Barrier Option
// For more inormation about this option formula see "The Complete Guide To Option Pricing Formulas" 2006
double BSAmericanApprox2002(char CallPutFlag, double S, double X, double T, double r, double b, double v)
{
// By Espen Gaarder Haug
    double OptionValue;
	
    if( CallPutFlag == 'c' )
    {
		    OptionValue = BSAmericanCallApprox2002(S, X, T, r, b, v);
    }
	else if (CallPutFlag == 'p' )  // Use the Bjerksund and Stensland put-call transformation
   {
		     OptionValue = BSAmericanCallApprox2002(X, S, T, r - b, -b, v);
   }
	
	return OptionValue;
    
};

double BSAmericanCallApprox2002(double S, double X, double T, double r, double b, double v)
{    
    double BInfinity, B0, OptionValue;
    double ht1, ht2, I1, I2;
    double alfa1, alfa2, Beta, t1;
    
    t1 = 0.5 * (sqrt(5.0) - 1.0) * T;
    
    if( b >= r )  // Never optimal to exersice before maturity
	{
            OptionValue = GBlackScholes('c', S, X, T, r, b, v);
	}
    else
	{
        
        Beta = (0.5 - b / ( v* v)) + sqrt(pow(b / (v * v) - 0.5, 2) + 2.0 * r / (v * v));
        BInfinity = Beta / (Beta - 1.0) * X;
        B0 = esmax(X, r / (r - b) * X);
        
        ht1 = -(b * t1 + 2.0 * v * sqrt(t1)) * X * X / ((BInfinity - B0) * B0);
        ht2 = -(b * T + 2.0 * v * sqrt(T)) * X * X / ((BInfinity - B0) * B0);
        I1 = B0 + (BInfinity - B0) * (1.0 - exp(ht1));
        I2 = B0 + (BInfinity - B0) * (1.0 - exp(ht2));
        alfa1 = (I1 - X) * pow(I1, -Beta);
        alfa2 = (I2 - X) * pow(I2,-Beta);
    
        if (S >= I2 )
		{
            OptionValue = S - X;
		}
        else
		{
            OptionValue = alfa2 * pow(S , Beta) - alfa2 * phi(S, t1, Beta, I2, I2, r, b, v) 
                + phi(S, t1, 1, I2, I2, r, b, v) - phi(S, t1, 1, I1, I2, r, b, v) 
                - X * phi(S, t1, 0, I2, I2, r, b, v) + X * phi(S, t1, 0, I1, I2, r, b, v) 
                + alfa1 * phi(S, t1, Beta, I1, I2, r, b, v) - alfa1 * ksi(S, T, Beta, I1, I2, I1, t1, r, b, v) 
                + ksi(S, T, 1, I1, I2, I1, t1, r, b, v) - ksi(S, T, 1, X, I2, I1, t1, r, b, v) 
                - X * ksi(S, T, 0, I1, I2, I1, t1, r, b, v) + X * ksi(S, T, 0, X, I2, I1, t1, r, b, v);
           
        }
   }
   return OptionValue;
   
};

double phi(double S, double T, double gamma, double H, double i, double r, double b, double v) 
{

	double lambda, kappa, d;
    
    lambda = (-r + gamma * b + 0.5 * gamma * (gamma - 1.0) * v * v) * T;
    d = -(log(S / H) + (b + (gamma - 0.5) * v * v) * T) / (v * sqrt(T));
    kappa = 2.0 * b / (v * v) + 2.0 * gamma - 1.0;
	
    return exp(lambda) * pow(S,gamma) * (CND(d) - pow(i / S, kappa) * CND(d - 2 * log(i / S) / (v * sqrt(T))));

};


double ksi(double S, double T2, double gamma, double H, double I2, double I1, double t1, double r, double b, double v)
{

    double e1, e2, e3, e4;
    double f1, f2, f3, f4;
    double rho, kappa, lambda;
    
    e1 = (log(S / I1) + (b + (gamma - 0.5) * v * v) * t1) / (v * sqrt(t1));
    e2 = (log(I2 * I2 / (S * I1)) + (b + (gamma - 0.5) * v * v) * t1) / (v * sqrt(t1));
    e3 = (log(S / I1) - (b + (gamma - 0.5) * v * v) * t1) / (v * sqrt(t1));
    e4 = (log(I2 * I2 / (S * I1)) - (b + (gamma - 0.5) * v * v) * t1) / (v * sqrt(t1));
    
    f1 = (log(S / H) + (b + (gamma - 0.5) * v * v) * T2) / (v * sqrt(T2));
    f2 = (log(I2 * I2 / (S * H)) + (b + (gamma - 0.5) * v * v) * T2) / (v * sqrt(T2));
    f3 = (log(I1 * I1 / (S * H)) + (b + (gamma - 0.5) * v * v) * T2) / (v * sqrt(T2));
    f4 = (log((S * I1 * I1) / (H * I2 * I2)) + (b + (gamma - 0.5) * v * v) * T2) / (v * sqrt(T2));
    
    rho = sqrt(t1 / T2);
    lambda = -r + gamma * b + 0.5 * gamma * (gamma - 1.0) * v * v;
    kappa = 2.0 * b / (v * v) + 2 * gamma - 1.0;
    
    return exp(lambda * T2) * pow(S,gamma) * (CBND(-e1, -f1, rho) - pow(I2 / S,kappa) * CBND(-e2, -f2, rho) 
            - pow(I1 / S, kappa) * CBND(-e3, -f3, -rho) + pow(I1 / I2, kappa) * CBND(-e4, -f4, -rho));


};




// First-then-barrier options
double FirstThenBarrier(int TypeFlag, double S, double X, double L, double U, double T, double r, double v)
{            
//By Espen Gaarder Haug
            double k;
            
            k = 0.0;
               
            if (TypeFlag == 1)  // First-down-then-up-and-in call
                return  X / L * StandardBarrier("pdi", S, L * L / X, L * L / U, k, T, r, 0, v);
            else if (TypeFlag == 2)  // First-up-then-down-and-in call
                return X / U * StandardBarrier("pui", S, U  * U / X, U * U / L, k, T, r, 0, v);
            else if (TypeFlag == 3)  // First-down-then-up-and-in put
                return X / L * StandardBarrier("cdi", S, L * L / X, L * L / U, k, T, r, 0, v);
            else if (TypeFlag == 4)  // First-up-then-down-and-in put
                return X / U * StandardBarrier("cui", S, U * U / X, U * U / L, k, T, r, 0, v);
            else if (TypeFlag == 5)  // First-down-then-up-and-out call
                return GBlackScholes('c', S, X, T, r, 0, v) - X / L * StandardBarrier("pdi", S, L * L / X, L * L / U, k, T, r, 0, v);
            else if (TypeFlag == 6)  // First-up-then-down-and-out call
                return GBlackScholes('c', S, X, T, r, 0, v) - X / U * StandardBarrier("pui", S, U * U / X, U * U / L, k, T, r, 0, v);
            else if (TypeFlag == 7)  // First-down-then-up-and-out put
                return GBlackScholes('p', S, X, T, r, 0, v) - X / L * StandardBarrier("cdi", S, L * L / X, L * L / U, k, T, r, 0, v);
            else if (TypeFlag == 8)  // First-up-then-down-and-out put
                return GBlackScholes('p', S, X, T, r, 0, v) - X / U * StandardBarrier("cui", S, U * U  / X, U * U / L, k, T, r, 0, v);
};

// Dual Double barrier options formula
double DoubleBarrierHaugCP(char TypeFlag[], double S, double X, double L, double U, double T, double r, double v)
{        
// By Espen Gaarder Haug

        double sum;
        
        if (TypeFlag == "ci")  //Put down and in call up and in
		{
            sum = StandardBarrier("cui", S, X, U, 0, T, r, 0, v) + StandardBarrier("pdi", S, X, L, 0, T, r, 0, v) 
                - X / U * StandardBarrier("pui", S, U * U / X, U * U / L, 0, T, r, 0, v) - X / L * StandardBarrier("cdi", S, L * L / X, L * L / U, 0, T, r, 0, v) 
                + U / L * StandardBarrier("cdi", S, L * L * X / (U * U), pow(L,3) / (U * U), 0, T, r, 0, v) + L / U * StandardBarrier("pui", S, U * U * X / (L * L), pow(U,3) / (L * L), 0, T, r, 0, v) 
                - L * X / (U * U) * StandardBarrier("pui", S, pow(U,4) / (L * L * X), pow(U,4) / pow(L,3), 0, T, r, 0, v) - U * X / (L * L) * StandardBarrier("cdi", S, pow(L,4) / (U * U * X), pow(L , 4) / pow(U,3), 0, T, r, 0, v) 
                + U * U / (L * L) * StandardBarrier("cdi", S, pow(L,4) * X / pow(U,4), pow(L,5) / pow(U,4), 0, T, r, 0, v) + L * L / (U * U) * StandardBarrier("pui", S, pow(U,4) * X / pow(L,4), pow(U,5) / pow(L,4), 0, T, r, 0, v);
         }
		else if (TypeFlag == "pi" ) // Call down and in  put up and in
		{
                sum = StandardBarrier("pui", S, X, U, 0, T, r, 0, v) + StandardBarrier("cdi", S, X, L, 0, T, r, 0, v) 
                - X / U * StandardBarrier("cui", S, U * U / X, U * U / L, 0, T, r, 0, v) - X / L * StandardBarrier("pdi", S, L * L / X, L * L / U, 0, T, r, 0, v) 
                + U / L * StandardBarrier("pdi", S, L * L * X / (U * U), pow(L,3) / (U * U), 0, T, r, 0, v) + L / U * StandardBarrier("cui", S, U * U * X / (L * L), pow(U,3) / (L * L), 0, T, r, 0, v) 
                - L * X / (U * U) * StandardBarrier("cui", S, pow(U,4) / (L * L * X), pow(U,4) / pow(L,3), 0, T, r, 0, v) - U * X /(L * L) * StandardBarrier("pdi", S, pow(L,4) / (U * U * X), pow(L,4) / pow(U,3), 0, T, r, 0, v) 
                + U * U /(L * L) * StandardBarrier("pdi", S, pow(L,4) * X / pow(U,4), pow(L,5) / pow(U,4), 0, T, r, 0, v) + L * L / (U * U) * StandardBarrier("cui", S, pow(U,4) * X / pow(L,4), pow(U,5) / pow(L,4), 0, T, r, 0, v);
       }
		
          
        sum = esmin(sum, esmax(GBlackScholes('c', S, X, T, r, 0, v), GBlackScholes('p', S, X, T, r, 0, v)));
       
         if (TypeFlag == "ci" || TypeFlag == "pi" )
		 {
		return sum;
		 }
        else if( TypeFlag == "co" ) //Knock out straddle, knock out call to upside and put to down side, here you always keep one option
        {
		return GBlackScholes('c', S, X, T, r, 0, v) + GBlackScholes('p', S, X, T, r, 0, v) - DoubleBarrierHaugCP("ci", S, X, L, U, T, r, v);
        }
		else if( TypeFlag == "po" ) //Knock out straddle, knock out call to downside and put to up side side
        {
		return GBlackScholes('c', S, X, T, r, 0, v) + GBlackScholes('p', S, X, T, r, 0, v) - DoubleBarrierHaugCP("pi", S, X, L, U, T, r, v);
		}
    
};

//  Double barrier options formula using put-call barrier symmetry
double DoubleBarrierHaug(char TypeFlag[], double S, double X, double L, double U, double T, double r, double v) 
{        
// By Espen Gaarder Haug

        double sum, DoubleBarrier;
        
		if  (strcmp(TypeFlag,"ci") == 0)
		{
             sum = StandardBarrier("cui", S, X, U, 0, T, r, 0, v) + StandardBarrier("cdi", S, X, L, 0, T, r, 0, v) 
                    - X / U * StandardBarrier("pui", S, U *U / X, U *U / L, 0, T, r, 0, v) - X / L * StandardBarrier("pdi", S, L * L / X, L * L / U, 0, T, r, 0, v) 
                    + U / L * StandardBarrier("cdi", S, L * L * X / (U * U), pow(L , 3) / (U * U), 0, T, r, 0, v) + L / U * StandardBarrier("cui", S, U * U * X / (L * L), pow(U,3) / (L * L), 0, T, r, 0, v) 
                    - L * X / (U * U) * StandardBarrier("pui", S, pow(U,4) / (L * L * X), pow(U,4) / pow(L,3), 0, T, r, 0, v) - U * X /( L * L) * StandardBarrier("pdi", S, pow(L,4) / (U * U * X), pow(L,4) / pow(U,3), 0, T, r, 0, v) 
                    + (U * U) / (L * L) * StandardBarrier("cdi", S, pow(L,4) * X / pow(U,4), pow(L,5) / pow(U,4), 0, T, r, 0, v) + L * L / (U * U) * StandardBarrier("cui", S, pow(U,4) * X / pow(L,4), pow(U,5) / pow(L,4), 0, T, r, 0, v);
                DoubleBarrier = esmin(sum, GBlackScholes('c', S, X, T, r, 0, v));
		}
		else if  (strcmp(TypeFlag,"pi") == 0)
		{
                sum = StandardBarrier("pui", S, X, U, 0, T, r, 0, v) + StandardBarrier("pdi", S, X, L, 0, T, r, 0, v) 
                    - X / U * StandardBarrier("cui", S, U * U / X, U * U / L, 0, T, r, 0, v) - X / L * StandardBarrier("cdi", S, L * L / X, L * L / U, 0, T, r, 0, v) 
                    + U / L * StandardBarrier("pdi", S, L * L * X / (U * U), pow(L,3) / (U * U), 0, T, r, 0, v) + L / U * StandardBarrier("pui", S, U * U * X / (L * L), pow(U,3) / (L * L), 0, T, r, 0, v) 
                    - L * X / (U * U) * StandardBarrier("cui", S, pow(U,4) / (L * L * X), pow(U,4) / pow(L,3), 0, T, r, 0, v) - U * X / (L * L) * StandardBarrier("cdi", S, pow(L,4) / (U * U * X), pow(L,4) / pow(U,3), 0, T, r, 0, v) 
                    + U * U / (L * L) * StandardBarrier("pdi", S, pow(L,4) * X / pow(U,4), pow(L,5) / pow(U,4), 0, T, r, 0, v) + L * L / (U * U) * StandardBarrier("pui", S, pow(U,4) * X / pow(L,4), pow(U,5) / pow(L,4), 0, T, r, 0, v);
               DoubleBarrier = esmin(sum, GBlackScholes('p', S, X, T, r, 0, v));
		} 
    
    else if  (strcmp(TypeFlag,"co") == 0)
	{
        DoubleBarrier = GBlackScholes('c', S, X, T, r, 0, v) - DoubleBarrierHaug("ci", S, X, L, U, T, r, v);
	}
    else if  (strcmp(TypeFlag,"po") == 0)
	{
        DoubleBarrier = GBlackScholes('p', S, X, T, r, 0, v) - DoubleBarrierHaug("pi", S, X, L, U, T, r, v);
	}
   
   return   DoubleBarrier;       
};


// Standard barrier options
double StandardBarrier(char TypeFlag[], double S, double X, double H, double K, double T, 
            double r, double b, double v)
{    
// By Espen Gaarder Haug
    double mu, lambda, X1, X2, y1, y2, Z;
    
    int eta, phi ;     //Binary variable that can take the value of 1 or -1
     
    double f1 ;   
	double f2 ;   
	double f3 ;   
    double f4 ;   
    double f5 ;   
    double f6 ;  
	
    mu = (b - v * v / 2) / (v*v);
    lambda = sqrt(mu * mu + 2 * r / (v*v) );
    X1 = log(S / X) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T);
    X2 = log(S / H) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T);
    y1 = log(pow(H,2) / (S * X)) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T);
    y2 = log(H / S) / (v * sqrt(T)) + (1 + mu) * v * sqrt(T);
    Z = log(H / S) / (v * sqrt(T)) + lambda * v * sqrt(T);
    
    if (TypeFlag == "cdi" || TypeFlag == "cdi"){
        eta = 1;
        phi = 1;
    }    
    else if (TypeFlag == "cui" || TypeFlag == "cuo") {
        eta = -1;
        phi =  1;
    }
    else if (TypeFlag == "pdi" || TypeFlag == "pdo" ){
        eta =  1;
        phi = -1;
    }
    else if (TypeFlag == "pui" || TypeFlag == "puo" ){
        eta = -1;
        phi = -1;
    };
    
    f1 = phi * S * exp((b - r) * T) * CND(phi * X1) - phi * X * exp(-r * T) * CND(phi * X1 - phi * v * sqrt(T));
    f2 = phi * S * exp((b - r) * T) * CND(phi * X2) - phi * X * exp(-r * T) * CND(phi * X2 - phi * v * sqrt(T));
    f3 = phi * S * exp((b - r) * T) * pow(H / S , 2 * (mu + 1)) * CND(eta * y1) - phi * X * exp(-r * T) * pow(H / S, 2 * mu) * CND(eta * y1 - eta * v * sqrt(T));
    f4 = phi * S * exp((b - r) * T) * pow(H / S, 2 * (mu + 1)) * CND(eta * y2) - phi * X * exp(-r * T) * pow(H / S, 2 * mu) * CND(eta * y2 - eta * v * sqrt(T));
    f5 = K * exp(-r * T) * (CND(eta * X2 - eta * v * sqrt(T)) - pow(H / S ,2 * mu) * CND(eta * y2 - eta * v * sqrt(T)));
    f6 = K * (pow(H / S,mu + lambda) * CND(eta * Z) + pow(H / S, mu - lambda) * CND(eta * Z - 2 * eta * lambda * v * sqrt(T)));
    
    
    if (X > H)  
    {
            if (TypeFlag == "cdi")
                return f3 + f5;
            else if (TypeFlag == "cui")
                return f1 + f5;
              
            else if (TypeFlag == "pdi")
                return f2 - f3 + f4 + f5;
                
            else if (TypeFlag == "pui")
                return f1 - f2 + f4 + f5;
                
            else if (TypeFlag == "cdo")
                return f1 - f3 + f6;
                
            else if (TypeFlag == "cuo")
                return f6;
                
            else if (TypeFlag == "pdo")
                return f1 - f2 + f3 - f4 + f6;
                
            else if (TypeFlag == "puo")
                return f2 - f4 + f6;
     }      
    else if (X < H) 
    {
            if (TypeFlag ==  "cdi")
                return f1 - f2 + f4 + f5;
                
            else if (TypeFlag == "cui")
                return f2 - f3 + f4 + f5;
                
            else if (TypeFlag == "pdi")
                return f1 + f5;
                
            else if (TypeFlag == "pui")
                return f3 + f5;
                
            else if (TypeFlag == "cdo")
                return f2 + f6 - f4;
                
            else if (TypeFlag == "cuo")
                return f1 - f2 + f3 - f4 + f6;
                
            else if (TypeFlag == "pdo")
                return f6;
                
            else if (TypeFlag == "puo")
                return f1 - f3 + f6;
    };
};




// *************************************************************************************
//										DISTRIBUTION FUNCTION
// *************************************************************************************

// The normal distribution function
double ND(double x)
{
// By Espen Gaarder Haug

    return (1.0/sqrt(Pi + Pi)) * exp(-x*x/2.0);
};


// Cummulative double precision algorithm based on Hart 1968
// Based on an implementation by Graeme West
double CND(double X ) 
{
// By Espen Gaarder Haug

    double y , Exponential , SumA , SumB, CNDValue ;
    
    y = fabs(X);
    if( y > 37 )
	{
        CNDValue = 0.0;
	}
    else
	{
        Exponential = exp(-y *y / 2.0);
		
        if (y < 7.07106781186547)
		{
            SumA = 0.0352624965998911 * y + 0.700383064443688;
            SumA = SumA * y + 6.37396220353165;
            SumA = SumA * y + 33.912866078383;
            SumA = SumA * y + 112.079291497871;
            SumA = SumA * y + 221.213596169931;
            SumA = SumA * y + 220.206867912376;
            SumB = 0.0883883476483184 * y + 1.75566716318264;
            SumB = SumB * y + 16.064177579207;
            SumB = SumB * y + 86.7807322029461;
            SumB = SumB * y + 296.564248779674;
            SumB = SumB * y + 637.333633378831;
            SumB = SumB * y + 793.826512519948;
            SumB = SumB * y + 440.413735824752;
            CNDValue = Exponential * SumA / SumB;
		}
        else
		{
            SumA = y + 0.65;
            SumA = y + 4.0 / SumA;
            SumA = y + 3.0 / SumA;
            SumA = y + 2.0 / SumA;
            SumA = y + 1.0 / SumA;
            CNDValue = Exponential / (SumA * 2.506628274631);
        }
  }
  
  if( X > 0)
  {
	CNDValue = 1.0 - CNDValue;
   }  
   
   return CNDValue;
   
};



// Inverse Cummulative Normal Distriution Approximation, based on Moro 1995
double CNDEV(double u)
{

    double x, r;
	double a[4] = {2.50662823884, -18.61500062529, 41.39119773534, -25.44106049637};
	double b[4] = {-8.4735109309, 23.08336743743, -21.06224101826, 3.13082909833};
	double c[9] = {0.337475482272615, 0.976169019091719, 0.160797971491821, 2.76438810333863E-02,
					3.8405729373609E-03, 3.951896511919E-04, 3.21767881767818E-05, 2.888167364E-07, 3.960315187E-07};

    x = u - 0.5;

    if (fabs(x) < 0.42)
	{
        r = x * x;
        r = x * (((a[3] * r + a[2]) * r + a[1]) * r + a[0])
        / ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
        return  r;
    }
    r = u;
    if (x >= 0.0) r = 1.0 - u;
    r = log(-log(r));
    r = c[0] + r * (c[1] + r * (c[2] + r * (c[3] + r * (c[4] +
        r * (c[5] + r * (c[6] + r * (c[7] + r * c[8])))))));
    if (x < 0)  r = -r;
    return r;
}


// The cumulative bivariate normal distribution function
double CBND(double a, double b, double rho) 
{

    double rho1 , rho2 , delta ,a1 , b1 , Sum ;
    int i , j ;
    
    double X[5] = {0.24840615, 0.39233107, 0.21141819, 0.03324666, 0.00082485334};
    double y[5] = {0.10024215, 0.48281397, 1.0609498, 1.7797294, 2.6697604};
    a1 = a / sqrt(2.0 * (1.0 - rho * rho));
    b1 = b / sqrt(2.0 * (1.0 - rho * rho));
    
    if ((a <= 0.0) && (b <= 0.0) && (rho <= 0.0) ){
        Sum = 0;
        for (i = 0 ; i<= 4; i++){
            for (j = 0 ; j<= 4; j++){
                Sum += X[i] * X[j] * exp(a1 * (2 * y[i] - a1) 
                + b1 * (2 * y[j] - b1) + 2 * rho * (y[i] - a1) * (y[j] - b1));
            };
        };
       Sum= sqrt(1 - rho * rho) / Pi * Sum;
      return Sum;
    }
    else if ((a <= 0.0) && (b >= 0.0) && (rho >= 0.0)){
        return CND(a) - CBND(a, -b, -rho);
    }
    else if ((a >= 0.0) && (b <= 0.0) && (rho >= 0.0)){
        return CND(b) - CBND(-a, b, -rho);
    }
    else if (a >= 0.0 && b >= 0.0 && rho <= 0.0){
        return CND(a) + CND(b) - 1 + CBND(-a, -b, rho);
    }
    else if (a * b * rho > 0.0){
        rho1 = (rho * a - b) * sgn(a) / sqrt(a * a - 2 * rho * a * b + b * b);
        rho2 = (rho * b - a) * sgn(b) / sqrt(a * a - 2 * rho * a * b + b * b);
        delta = (1 - sgn(a) * sgn(b)) / (double) 4.0;
        return CBND(a, 0.0, rho1) + CBND(b, 0.0, rho2) - delta;
    };
    return -1000.0; // Just to spot errors in input or function
};



// Simple maximum function (PRIVATE FUNCTION)
double esmax(double x, double y)
{
if ( x>=y )
 return x;
else
 return y;
};


// Simple minimum function (PRIVATE FUNCTION)
double esmin(double x, double y)
{
if ( x>=y )
 return y;
else
 return x;
};