#pragma once

#include "stdafx.h"
#include "CardIn.h"
#include "AuxFunc.h"
#include <eigen3/Eigen/Dense>

using namespace std;

class OnePartInt : AuxFunc
{
	// Class to evaluate one and two electron atomic integrals. 
	// It normalizes the basis functions as well.
public:
	// Relativistic instance.
	OnePartInt(CardIn & input, vector<Gspinor> & Large, vector<Gspinor> & Small);
	// Non-relativistic instance.
	OnePartInt(CardIn & input, vector<Gspinor> & Large);
	// LL and SS Overlap integral.
	double S_LL(int a, int b, Gspinor & M);
	double S_SS(int a, int b, Gspinor & M);
	//MatrixXd Coef_L;//Large expansion coefficients
	//MatrixXd Coef_S;//Large expansion coefficients
	// M*[L,a,r] and M[S,b,r] Kinetic matrix element.
	double Kinetic_SL(int a, int b, Gspinor & M_S, Gspinor & M_L);
	// LL and SS nuclear attraction integral.
	double Nuclear_LL(int a, int b, Gspinor & M);
	double Nuclear_SS(int a, int b, Gspinor & M);
	// Point-like nuclei.
	double Nuclear_LL_Point(int a, int b, Gspinor & M);
	double Nuclear_SS_Point(int a, int b, Gspinor & M);
	// Non-Relativistic Kinetic energy integral.
	double Kinetic_NR(int a, int b, Gspinor & M);
	// Non-Relativistic point like nuclear integral.
	double Nuclear_NR(int a, int b, Gspinor & M);

	~OnePartInt();
private:
	CardIn & input;
	// Parameter for nuclear attraction integral;
	double zeta;
};

