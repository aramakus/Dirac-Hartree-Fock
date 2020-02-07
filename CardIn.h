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
#include <fstream>
#include <sys/stat.h>
#include <string>
#include <vector>
#include <sstream>
#include "AuxFunc.h"

using namespace std;

class CardIn
{
	// An input reader. Reads and stores all the information needed for subsequent calculations.
	// Also initializes the basis setexponents and quantum numbers.
public:
	CardIn(string filename, vector<Gspinor> & Large, vector<Gspinor> & Small, ofstream & log);
	~CardIn();

	string label;// Name for the calculation.	
	string Hamiltonian;// Hamiltonian model.
	string BasisModel;
	string NewOld;
	int Z;
	double Amass;
private:
	int REL_assign_closed(vector<Gspinor> & Large, vector<Gspinor> & Small, string line, ofstream & log); // assign relativistic closed shell orbitals.
	int NR_assign_closed(vector<Gspinor> & Large, string line, ofstream & log);  // assign non-relativistic closed shell orbitals.
	int REL_assign_opened(vector<Gspinor> & Large, vector<Gspinor> & Small, string line, ofstream & log); // assign relativistic opened shell orbitals.
	int NR_assign_opened(vector<Gspinor> & Large, string line, ofstream & log); // assign non-relativistic opened shell orbitals.
	void read_opened_line(string line, vector<int> & n, vector<int> & l, vector<int> & occ);
};