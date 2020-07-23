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

// ======================= Atom =======================
//
// Program for relativistic (DHF) and non-relativistic
// (HF) atomic calculations. Gaussian basis is used 
// for NR and G-spinor basis for DHF routines (see 
// "Handbook of Molecular Physics and Quantum 
// Chemistry" 2:6:31, p. 696–716 by H. M . Quiney). 
// Open shell atoms are treated using average over
// configuration approximation (see "Case Studies in
// Atomic Physics IV" by E. W. McDaniel and 
// M. R. C. McDowell, 1975). Finite nuclear effects
// are included in the calculation.
//
// ====================================================
#include "CardIn.h"
#include "LinearAlgebra.h"
#include <ctime>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

/*
Atomic code for atomic structure calculations using finite basis Hartree-Fock approximation.
Notations and variable names follow I. Grant "Relativistic Qunatum Theory of Atoms".
*/
int main(int argc, char *argv[])
{
  // Read input file content.
	clock_t t = clock();
	string logname = "log";
	string filename = "";
	if (argc > 1) {
		filename = argv[1];
		size_t lastdot = filename.find_last_of(".");
		if (lastdot != string::npos) filename.substr(0, lastdot);
		logname = logname + filename;
	}
	else logname = logname + ".txt";
	ofstream log(logname);
	if (argc <= 1) {
		log << "no input file found. Exiting..." << endl;
		return 1;
	}

	// Initialize the basis set.
	vector<Gspinor> Large; // Large component, radial wave function.
	vector<Gspinor> Small; // Small component (rel only), radial wave function.

	// Read an input file.
	CardIn input(filename, Large, Small, log);

  // Iterative solution of Hartree-Fock equations.
	LinearAlgebra Solver(input);
	if (input.Hamiltonian == "DHF")	{
		Solver.SolveFock(Large, Small);
		Solver.IterateHF(Large, Small);		
	}
	else {
		Solver.SolveFock(Large);
		Solver.IterateHF(Large);
	}

	log.flush();
	log.close();

	t = clock() - t; 
	printf("=============================================================================\n");
	printf("===                 total execution time : %8.2f s                     ===\n", ((float)t)/CLOCKS_PER_SEC);
	printf("===                           Normal Exit                                 ===\n");
	printf("=============================================================================\n");

	return 0;
}
