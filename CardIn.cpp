#include "CardIn.h"
#include <map>

CardIn::CardIn(string filename, vector<Gspinor> & Large, vector<Gspinor> & Small, ofstream & log)
{
	if (filename == "") log << "No input file identified. Quitting..." << endl;
	else{
		log << "Input file: " << filename << endl;

		ifstream infile(filename);
		
		map<string, vector<string>> FileContent;
		string comment = "//";
		string curr_key;
		
		while (!infile.eof() && infile.is_open()) {
			string line;
			getline(infile, line);

			if (!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);

			if (!line.compare(0, 2, comment)) continue;
			if (line == "") continue;
			if (!line.compare(0, 1, "#")) {
				if ( FileContent.find(line) == FileContent.end() ) {
					FileContent[line] = vector<string>(0);
				}
				curr_key = line;
			} else {
				FileContent[curr_key].push_back(line);
			}
		}

		for (int n = 0; n < FileContent["#GENERAL"].size(); n++) {
			stringstream stream(FileContent["#GENERAL"][n]);
			if (n == 0) stream >> label;
			if (n == 1) stream >> Z;
			if (n == 2) stream >> Amass;
			if (n == 3) stream >> Hamiltonian;
			if (n == 4) stream >> NewOld;
			if (n == 5) stream >> BasisModel;
		}

		for (int n = 0; n < FileContent["#STRUCTURE"].size(); n++) {
			stringstream stream(FileContent["#STRUCTURE"][n]);
			if (n != 0) continue;

			if (Hamiltonian == "DHF") {
				REL_assign_closed(Large, Small, FileContent["#STRUCTURE"][n], log);
				REL_assign_opened(Large, Small, FileContent["#STRUCTURE"][n], log);
			} else {
				NR_assign_closed(Large, FileContent["#STRUCTURE"][n], log);
				NR_assign_opened(Large, FileContent["#STRUCTURE"][n], log);
			}			
		}

    if (BasisModel == "GEOMETRIC") {
      // User defined basis.
      int num_funcs = 0;
      double alpha = 0., beta = 0.;
      vector<double> lambda;
      for (int L = 0; L < FileContent["#BASIS"].size(); L++) {
        stringstream stream(FileContent["#BASIS"][L]);
        stream >> num_funcs >> alpha >> beta;

        lambda.clear();
        lambda.resize(num_funcs);
        lambda[0] = alpha;
        for (int n = 1; n < num_funcs; n++) lambda[n] = lambda[n-1]*beta;
        for (int k = 0; k < Large.size(); k++) {
          if (Large[k].l != L) continue;

          Large[k].lambda = lambda;
          if (Hamiltonian == "DHF") Small[k].lambda = lambda;
        }
        
      }
    } else {
      // Geometric basis set.
      int num_lines = 1;
      int curr_line = 1;
      int L = 0;
      double val = 0;
      stringstream stream(FileContent["#BASIS"][0]);
      stream >> num_lines;
      for (int n = 1; n < FileContent["#BASIS"].size(); n++) {
        stringstream stream(FileContent["#BASIS"][n]);
			  if (curr_line <= num_lines) {
          curr_line++;
          stream >> val;
          for (int k = 0; k < Large.size(); k++) {
            if (Large[k].l != L) continue;

            Large[k].lambda.push_back(val);
            if (Hamiltonian == "DHF") Small[k].lambda.push_back(val);
          }
        } else {
          curr_line = 1;
          stream >> num_lines;
          L++; // Implicit assumption about ordering User defined basis set.
        }

		  }
    }

	}
}


int CardIn::REL_assign_closed(vector<Gspinor> & Large, vector<Gspinor> & Small, string line, ofstream & log)
{ 
	// Given a non-relativistic value of angular momentum L
	// Initialize RELATIVISTIC orbiral.

	if (line.compare(0, 1, "[")) {
		log << "No closed shell notations in the input file." << endl;
		return 1;
	}

	map<string, int> shells = {{"[He]" , 1}, {"[Ne]", 2}, {"[Ar]", 3}, {"[Kr]", 4}, {"[Xe]", 5}, {"[Rn]", 6}};
	
	try {
		int N_max = shells.at(line.substr(0, 4));
		int Occ = 0;
		
		// g orbitals are not occupies in noble gasses.
		for (int L = 0; L < N_max; L++) {
			Occ = 2*L + 2;

			Large.push_back(Gspinor());
			Large.back().l = L;
			Large.back().j = L + 0.5;
			Large.back().k = -L-1;
			Small.push_back(Gspinor());
			Small.back().l = L;
			Small.back().j = L + 0.5;
			Small.back().k = -L-1;
			
			for (int N = L + 1; N <= N_max; N++) {
				Large.back().occup.push_back(Occ);
				Small.back().occup.push_back(Occ);
			}

			if (L > 0) {
				Occ = 2*L;
				Large.push_back(Gspinor());
				Large.back().l = L;
				Large.back().j = L - 0.5;	
				Large.back().k = L;	
				Small.push_back(Gspinor());
				Small.back().l = L;
				Small.back().j = L - 0.5;
				Small.back().k = L;
				for (int N = L + 1; N <= N_max; N++) {
					Large.back().occup.push_back(Occ);
					Small.back().occup.push_back(Occ);
				}
			}

			if (L == 1) N_max--; // 3d orbital is occupied only after 4s.
			if (L == 2) N_max--; // 4f orbital is occupied onlly after 6s.
		}
	}	catch (const out_of_range& oor) {
		log << "Unknown symbol for closed shells - [X]." << endl;
		return 1;
	}

	return 0;
}

int CardIn::NR_assign_closed(vector<Gspinor> & Large, string line, ofstream & log)
{

	if (line.compare(0, 1, "[")) {
		log << "No closed shell notations in the input file." << endl;
		return 1;
	}

	map<string, int> shells = {{"[He]" , 1}, {"[Ne]", 2}, {"[Ar]", 3}, {"[Kr]", 4}, {"[Xe]", 5}, {"[Rn]", 6}};

	try {
		int N_max = shells.at(line.substr(0, 4));
		int Occ = 0;

		if (N_max < 1) return 1;
		
		for (int L = 0; L < N_max; L++) {
			int Occ = 4*L + 2;
			Large.push_back(Gspinor());
			Large.back().l = L;
			Large.back().k = 0;
			for (int N = L + 1; N <= N_max; N++) Large.back().occup.push_back(Occ);
			if (L == 1) N_max--; // 3d orbital is occupied only after 4s.
			if (L == 2) N_max--; // 4f orbital is occupied onlly after 6s.
		}
	}	catch (const out_of_range& oor) {
		log << "Unknown symbol for closed shells - [X]." << endl;
		return 1;
	}

  for (int K = 0; K < Large.size(); K++) {
    Large[K].energy.resize(Large[K].occup.size());
    Large[K].C_expansion.resize(Large[K].occup.size(), vector<double>(Large[K].lambda.size(), 0.));
  }

	return 0;
}

int CardIn::REL_assign_opened(vector<Gspinor> & Large, vector<Gspinor> & Small, string line, ofstream & log)
{
	// Read all open shell orbitals.
	vector<int> N;
	vector<int> L;
	vector<int> occ;

	read_opened_line(line, N, L, occ);
	
	// Assign opened shell orbitals.
	double Occ = 0;
	for (int a = 0; a < N.size(); a++) {
		bool found = false;
		for (int i = 0; i < Large.size(); i++) {
			if (Large[i].l == L[a]) {
				found = true;
				Occ = occ[a] * 2.*abs(Large[i].k)/(4*L[a] + 2);
				Large[i].occup.push_back(Occ);
				Small[i].occup.push_back(Occ);
				if (occ[a] != 4*L[a] + 2) Large[i].is_open = true;
				if (occ[a] != 4*L[a] + 2) Small[i].is_open = true;
			}
		}
		if (!found) {
			Large.push_back(Gspinor());
			Large.back().l = L[a];
			Large.back().j = L[a] + 0.5;
			Large.back().k = -L[a]-1;
      Occ = occ[a]*2.*abs(Large.back().k)/(4*L[a] + 2);
			Large.back().occup.push_back(Occ);
			if (occ[a] != 4*L[a] + 2) Large.back().is_open = true;
			Small.push_back(Gspinor());
			Small.back().l = L[a];
			Small.back().j = L[a] + 0.5;
			Small.back().k = -L[a]-1;
			Small.back().occup.push_back(Occ);
			if (occ[a] != 4*L[a] + 2) Small.back().is_open = true;
			if (L[a] > 0) {
				Large.push_back(Gspinor());
				Large.back().l = L[a];
				Large.back().j = L[a] - 0.5;	
				Large.back().k = L[a];
        Occ = occ[a]*2.*abs(Large.back().k)/(4*L[a] + 2);
				Large.back().occup.push_back(Occ);	
				if (occ[a] != 4*L[a] + 2) Large.back().is_open = true;
				Small.push_back(Gspinor());
				Small.back().l = L[a];
				Small.back().j = L[a] - 0.5;
				Small.back().k = L[a];
				Small.back().occup.push_back(Occ);
				if (occ[a] != 4*L[a] + 2) Small.back().is_open = true;
			}						
		}

	}
	
  for (int K = 0; K < Large.size(); K++) {
    Large[K].energy.resize(Large[K].occup.size());
    Small[K].energy.resize(Small[K].occup.size());
    Large[K].C_expansion.resize(Large[K].occup.size(), vector<double>(Large[K].lambda.size(), 0.));
    Small[K].C_expansion.resize(Small[K].occup.size(), vector<double>(Small[K].lambda.size(), 0.));
  }
	return 0;
}

int CardIn::NR_assign_opened(vector<Gspinor> & Large, string line, ofstream & log)
{
	// Read all open shell orbitals.
	vector<int> N;
	vector<int> L;
	vector<int> occ;

	read_opened_line(line, N, L, occ);

	// Assign opened shell orbitals.
	for (int a = 0; a < N.size(); a++) {
		bool found = false;
		for (int i = 0; i < Large.size(); i++) {
			if (Large[i].l == L[a]) {
				found = true;
				Large[i].occup.push_back(occ[a]);
				if (occ[a] != 4*L[a] + 2) Large[i].is_open = true;
			}
		}
		if (!found) {
			Large.push_back(Gspinor());
			Large.back().l = L[a];
			Large.back().k = 0;
			Large.back().occup.push_back(occ[a]);
			Large.back().is_open = true;					
		}

	}
	
	return 0;
}

void CardIn::read_opened_line(string line, vector<int> & N, vector<int> & L, vector<int> & Occ)
{
	string copy_line = line;
	vector<string> orbitals;
	size_t open = line.find('(');
	size_t close = line.find(')');
	while (open != string::npos) {
		orbitals.push_back(copy_line.substr(open+1, close-1));
		copy_line = copy_line.substr(close+1, string::npos);

		open = copy_line.find('(');
		close = copy_line.find(')');
	}
	// Traslate them into list of N, L, and Occ.
	map<string, int> L_list = {{"s", 0}, {"p", 1}, {"d", 2}, {"f", 3}};

	for (int n = 0; n < orbitals.size(); n++) {
		int tmp = 0;
		open = orbitals[n].find(' ');
		// Occupancy.
		copy_line = orbitals[n].substr(open+1, string::npos);
		stringstream stream(copy_line);
		stream >> tmp;
		Occ.push_back(tmp);
		// Angular momentum.
		copy_line = orbitals[n].substr(open-1, 1);
		L.push_back(L_list.at(copy_line));
		// Prime quantum number.
		copy_line = orbitals[n].substr(0, open-1);
		stream << copy_line;
		stream >> tmp;
		N.push_back(tmp);
	}
}

CardIn::~CardIn()
{
}
