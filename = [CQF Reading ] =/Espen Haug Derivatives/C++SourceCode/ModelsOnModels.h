/*
 *  ModelsOnModels.h
 *
 *  Created by Espen2Admin on 10/5/06.
 *  Copyright 2007 Espen Gaarder Haug. All rights reserved.
 *
 */
 
// C++ Code (Progming done in Xcode 
// By Espen Gaarder Haug

#include<iostream.h>
#include<math.h>
#include <vector.h> // standard STL vector template


#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif 



double ConvertingToCCRate( double r , double Compoundings) ;
double inCGBlackScholes(char OutputFlag[], char CallPutFlag, double S, double X,double T, double r, double b, double v);
double GBlackScholes(char CallPutFlag, double S, double X, double T, double r, double b, double v);

// DELTA GREEKS
double GDelta(double CallPutFlag, double S,double X,double T ,double r,double b, double v);
double GDdeltaDvol(double S, double X, double T, double r, double b, double v);
double GDdeltaDtime(char CallPutFlag , double  S, double X , double T, double r, double b, double v);
double GElasticity(char CallPutFlag, double S, double X, double T, double r, double b, double v);
double  GStrikeFromDelta(char CallPutFlag, double S, double T, double r, double b, double v, double delta);
// GAMMA GREEKS
double GGamma(double S,double X,double T ,double r, double b, double v);
double GGammaP(double S, double X, double T, double r, double b, double v);
double GDgammaDvol(double S, double X, double T, double r, double b, double v) ;
double GDgammaPDvol(double S, double X, double T, double r, double b, double v);
double GDgammaDspot(double S, double X, double T, double r, double b, double v);
double GDgammaPDspot(double S, double X ,double T ,double r,double b,double v );
double GDgammaDtime(double S, double X, double T, double r, double b, double v);
double GDgammaPDtime(double S, double X, double T, double r, double b, double v);
double  GSaddleGamma(double X, double T, double r, double b, double v);
// VEGA GREEKS
double GVega(double S,double X,double T ,double r, double b, double v);
double GVegaP(double S, double X, double T, double r, double b, double v);
double GVegaLeverage(char CallPutFlag, double S, double X, double T, double r, double b, double v);
double GDvegaDvol(double S, double X, double T, double r, double b, double v);
double GDvegaPDvol(double S, double X, double T, double r, double b, double v);
double GDvegaDtime(double S, double X, double T, double r, double b, double v);
// THETA GREEKS
double GTheta(double CallPutFlag, double S,double X,double T ,double r, double b, double v);
double GThetaDriftLess(double S, double X, double T, double r, double b, double v);
// RATE GREEKS
double GRho(double CallPutFlag, double S,double X,double T ,double r, double b, double v);
double GRhoFO(char CallPutFlag, double S, double X, double T, double r, double b, double v) ;
double GPhi(double CallPutFlag, double S,double X,double T ,double r, double b, double v);
double GCarry(double CallPutFlag, double S,double X,double T ,double r, double b, double v);
// PROBABILITY GREEKS
double GInTheMoneyProbability(char CallPutFlag, double S, double X, double T, double b, double v);
double GStrikeDelta(char CallPutFlag, double S, double X, double T, double r, double b, double v);
double GDzetaDvol(char CallPutFlag, double S, double X, double T, double r, double b, double v);
double GDzetaDtime(char CallPutFlag, double S, double X, double T, double r, double b, double v);
double GBreakEvenProbability(char CallPutFlag, double S, double X, double T, double r, double b , double v);
double GRiskNeutralDensity(double S, double X, double T , double r, double b, double v);
// EXOTIC OPTIONS
double AmericanKnockInBarriers(char CallPutFlag, double S, double X, double H, double T, double r, double b, double v);
double StandardBarrier(char TypeFlag[], double S, double X, double H, double K, double T, double r, double b, double v);
double FirstThenBarrier(int TypeFlag, double S, double X, double L, double U, double T, double r, double v);
double DoubleBarrierHaug(char TypeFlag[], double S, double X, double L, double U, double T, double r, double v) ;
double DoubleBarrierHaugCP(char TypeFlag[], double S, double X, double L, double U, double T, double r, double v);

// ENERGY OPTIONS
double EnergySwaption(char CallPutFlag, double F, double X, double T, double Tb, double rj, double rb, int j, int n, double v);
double EnergySwaptionApproximation(char CallPutFlag, double F, double X , double T, double Tm, double re, double v) ;

//NUMERICAL
double CRRBinomial(char OutputFLag, char AmeEurFlag, char CallPutFlag, double S, double X, double T, double r, double b, double v, int n);



double BSAmericanApprox2002(char CallPutFlag, double S, double X, double T, double r, double b, double v);
double BSAmericanCallApprox2002(double S, double X, double T, double r, double b, double v);
double phi(double S, double T, double gamma, double H, double i, double r, double b, double v);
double ksi(double S, double T2, double gamma, double H, double I2, double I1, double t1, double r, double b, double v);

double inGBlackScholes(char OutputFlag, char CallPutFlag, double S, double X, double T, double r, double v);
double esmax(double x, double y);
double esmin(double x, double y);
double CNDEV(double u);
double CND(double X);
double ND(double x);
double CBND(double a, double b, double rho);


int main()
//void main()
{
	double S, X, T, r, b, v;
	char CallPutFlag;
	CallPutFlag='c';
	S=100;
	X=100;
	T=0.5;
	r=0.1;
	b=0.1;
	v=0.3;
	
		cout <<"\nHere is the CGBlackScholes value: "<< inCGBlackScholes("p",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		
		cout <<"\nHere is the delta value: "<< inCGBlackScholes("d",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DdeltaDvol value: "<< inCGBlackScholes("dddv",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DdeltaDtime value: "<< inCGBlackScholes("dt",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the delta mirror strike value: "<< inCGBlackScholes("dmx",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the elasticity value: "<< inCGBlackScholes("e",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the saddle gamma value: "<< inCGBlackScholes("sg",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		
		cout <<"\nHere is the gamma value: "<< inCGBlackScholes("g",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the speed value: "<< inCGBlackScholes("s",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DgammaDvol Zomma value: "<< inCGBlackScholes("gv",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DgammaDtime  value: "<< inCGBlackScholes("gt",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		
		cout <<"\nHere is the gammaP : "<< inCGBlackScholes("gp",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the gammaP speed: "<< inCGBlackScholes("gps",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the gammaPvol ZommaP: "<< inCGBlackScholes("gpv",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DgammaPDtime value: "<< inCGBlackScholes("gpt",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		
		cout <<"\nHere is the vega value: "<< inCGBlackScholes("v",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DvegazDtime value: "<< inCGBlackScholes("vt",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the Vomma value: "<< inCGBlackScholes("dvdv",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the vegaP value: "<< inCGBlackScholes("vp",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the VommaP value: "<< inCGBlackScholes("vpv",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the vega leverage value: "<< inCGBlackScholes("vl",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;


		cout <<"\nHere is the rho value: "<< inCGBlackScholes("r",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the futures rho value: "<< inCGBlackScholes("fr",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the carry rho value: "<< inCGBlackScholes("b",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the phi rho2 value: "<< inCGBlackScholes("f",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		
		cout <<"\nHere is the zeta in the money prob value: "<< inCGBlackScholes("z",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DzetaDvol in the money prob value: "<< inCGBlackScholes("zv",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the DzetaDtime in the money prob value: "<< inCGBlackScholes("zt",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the break even probability value: "<< inCGBlackScholes("bp",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;

		cout <<"\nHere is the RND value: "<< inCGBlackScholes("dxdx",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the strike delta value: "<< inCGBlackScholes("dx",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;

		cout <<"\nHere is the d1 value: "<< inCGBlackScholes("d1",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the d2 value: "<< inCGBlackScholes("d2",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;

		cout <<"\nHere is the nd1 value: "<< inCGBlackScholes("nd1",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the nd2 value: "<< inCGBlackScholes("nd2",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;

		cout <<"\nHere is the nd1 value: "<< inCGBlackScholes("CNDd1",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		cout <<"\nHere is the nd2 value: "<< inCGBlackScholes("CNDd2",CallPutFlag, S, X, T, r, b, v)<<"\n" << endl;
		
		//EXOTIC OPTIONS
		cout <<"\nHere is the FirstThenBarrier value: "<< FirstThenBarrier(1,  100,  100,  90,  110,  0.5,  0.1,  0.3) <<"\n" << endl;
		cout <<"\nHere is the DoubleBarrierHaug value: "<< DoubleBarrierHaug( "po",  100,  100,  90,  110,  0.5,  0.1,  0.2) <<"\n" << endl;
		cout <<"\nHere is the DoubleBarrierHaugCP value: "<< DoubleBarrierHaugCP( "ci",  100,  100,  90,  110,  0.5,  0.1,  0.2) <<"\n" << endl;
		cout <<"\nHere is the AmericanKnockInBarriers value: "<< AmericanKnockInBarriers( 'p',  100,  100, 110,  0.5,  0.05,  0.02,  0.2) <<"\n" << endl;
	
		//ENERGY OPTIONS
		cout <<"\nHere is the EnergySwaption value: "<<  EnergySwaption( 'p',  33,  35,  0.25,  0.5,  0.05,  0.04,  365,  90,  0.18) <<"\n" << endl;
		cout <<"\nHere is the EnergySwaptionApproximation value: "<<  EnergySwaptionApproximation( 'p',  33,  35 ,  0.5,  0.6726,  0.05,  0.18)  <<"\n" << endl;
	
		//CRR TREE USED IN NEGATIVE PROBABILITY CHAPTER
		cout <<"\nHere is the CRRBinomial value: "<<  CRRBinomial( 'p',  'a',  'c',  42,  40,  0.5,  0.1,  0.1,  0.2,  20)  <<"\n" << endl;
		
				
}



