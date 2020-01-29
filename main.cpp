#include "stdafx.h"
#include "CardIn.h"
#include "LinearAlgebra.h"
#include <ctime> 

using namespace std;

int main(int argc, char *argv[])
{
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
	vector<Gspinor> Large;
	vector<Gspinor> Small;

	// Read an input file.
	CardIn input(filename, Large, Small, log);

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
