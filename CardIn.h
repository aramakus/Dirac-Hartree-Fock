#pragma once
#include "stdafx.h"
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