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

#include "CardIn.h"
#include <vector>
#include <eigen3/Eigen/Dense>
#include <map>

using namespace std;
using namespace Eigen;

class LinearAlgebra : AuxFunc
{
public:
	LinearAlgebra(CardIn & input);

	vector<MatrixXd> I; // Nuclear + Kinetic energy matrix.
	vector<MatrixXd> cG; // Direct and Exchange two-particle part of the Hamilotinian form.
	vector<MatrixXd> oG; // Direct and Exchange two-particle part of the Hamilotinian form.
	
  vector<Density> Dens; // Close and open shell density matrix.
	vector<Coulomb> G;
  
	vector<MatrixXd> Overlap; // Overlap matrix.
	
	void SolveFock(vector<Gspinor> & Large, vector<Gspinor> & Small); // Relativistic Hydrogenic ion. Initial guess for rel SelfConsist.
	void SolveFock(vector<Gspinor> & Large); // Non-Rel Hydrogenic ion. Initial guess for non-rel SelfConsist.
	void IterateHF(vector<Gspinor> & Large, vector<Gspinor> & Small); // DHF.
	void IterateHF(vector<Gspinor> & Large); // HF.

	int save_if_converged(vector<Gspinor> & Large, vector<Gspinor> & Small); // Save final result if convergence was achieved. TODO!!!

	~LinearAlgebra();
private:
	int max_iter = 50;
	double Etoller = 0.000000000001;// Relative energy tollerance.
	CardIn & input;
	vector<MatrixXd> EigVecs;
	vector<VectorXd> EigVals;
	vector<MatrixXd> oEigVecs;
	vector<VectorXd> oEigVals;

	void make_cG(vector<Gspinor> & Large, int K1, int K2);
	void make_oG(vector<Gspinor> & Large, int K);
  void make_cDens(double p, vector<Gspinor> & Large);
	void make_oDens(double p, vector<Gspinor> & Large);

  void make_cG(vector<Gspinor> & Large, vector<Gspinor> & Small, int K1, int K2);
	void make_oG(vector<Gspinor> & Large, vector<Gspinor> & Small, int K1, int K2);
  void make_cDens(double p, vector<Gspinor> & Large, vector<Gspinor> & Small);
	void make_oDens(double p, vector<Gspinor> & Large, vector<Gspinor> & Small);

  void sym_Matr(vector<Coulomb> & M);
  void REL_sym_Matr(vector<Coulomb> & M);
	vector<int> Dim;// List of basis set size for each Kappa(L).

	map<int, char> L_to_symb = {{0, 's'}, {1, 'p'}, {2, 'd'}, {3, 'f'}, {4, 'g'}, {5, 'h'}};
	string SymbAng(int kappa);

	bool converged = false;
};

// TODO LIST:
// 1. Parallelize make_cG and make_oG using OMP.
// 2. Add save_if_converged to save successfull calculations.