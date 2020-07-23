/*===========================================================================
This file is part of Atom.

    Atom is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Atom is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Atom.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
#pragma once

#include <vector>
#include "CardIn.h"
#include "AuxFunc.h"
#include <eigen3/Eigen/Dense>

using namespace std;

class OnePartInt : AuxFunc
{
	// Class to evaluate one electron atomic integrals. 
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

