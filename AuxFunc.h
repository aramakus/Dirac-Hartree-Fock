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
#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>

using namespace std;

namespace Constant
{
	const double Alpha = 0.0072973525698;
	const double Alpha2 = Constant::Alpha * Constant::Alpha;
	const double Fm = 1.8897261246 / 100000;
	const double Bohr = 0.529177244;
}

class AuxFunc
{
	// Class that contains auxillary functions for other calculations.
public:
	AuxFunc();
	// Integer and half integer gamma function, max x = 19.
	double Gamma(double x);
	// Odd R integral: Int_0^inf dr*r^n*exp(-lambda*r^2)*erf(zeta*r), max n=21
	double OddRINT(int n, double lambda, double zeta);
	// Incomplete Beta function: Int_0^z dt*t^(p-1/2)*(1-t)^(q-1),  
	double BetaFunc(int p, int q, double z);
	// Natural power, x^n.
	double Power(double x, int n);
	// Wigner3j for particular m1,m2,m3.
	// 3j symbol /l_a L l_b\ , la_max = lb_max = 5, L_max = 10 
	//           \ 0  0  0 /
	// 3j symbol /j_a L  j_b\ , 2ja_max = 2jb = 11, L_max = 22
	//           \0.5 0 -0.5/
	double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3);	

	~AuxFunc();
private:
	vector<double> gamma;// Gamma functions from 0 to 20 with increment of 0.5;
	double LogFF(float Num, float Denom);// Evaluates log(Num!/Denom!^k). Used in Wigner3j.
};

struct Gspinor
{
	vector<double> lambda;//exponent
	vector<double> Norm;
	int l;//angular momentum
	int k;//kappa
	vector<double> occup;//list of occupancies of orbitals with this Kappa(L)
  vector<double> energy;//single particle energy. 
	float j;//total orbital momentum
	int atom = 0;//atom to which it belongs
	bool is_open = false;
};

struct Density
{
	vector<vector<double>> cLL;
	vector<vector<double>> cSS;
	vector<vector<double>> cLS;
	vector<vector<double>> oLL;
	vector<vector<double>> oSS;
	vector<vector<double>> oLS;
};

struct Coulomb
{
	vector<vector<double>> cl; // closed.
	vector<vector<double>> op; // open.
};
