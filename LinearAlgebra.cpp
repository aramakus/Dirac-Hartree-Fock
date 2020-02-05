#include "LinearAlgebra.h"
#include "OnePartInt.h"
#include "TwoPartInt.h"
#include <iostream>
#include <fstream>

using namespace Constant;

LinearAlgebra::LinearAlgebra(CardIn & Input) : input(Input)
{
}

void LinearAlgebra::SolveFock(vector<Gspinor> & Large, vector<Gspinor> & Small)
{
	I.resize(Large.size());
	Overlap.resize(Large.size());
	Dim.resize(Large.size());

	for (int i = 0; i < Large.size(); i++) {
		Dim[i] = Large[i].lambda.size();
		I[i].resize(2*Dim[i],2*Dim[i]);
		Overlap[i] = MatrixXd::Zero(2*Dim[i],2*Dim[i]);
	}

	OnePartInt Calc(input, Large, Small);

	for (int K = 0; K < Large.size(); K++) {
		for (int i = 0; i < Dim[K]; i++) {
			for (int j = 0; j < Dim[K]; j++) {
				I[K](i+Dim[K], j) = Calc.Kinetic_SL(i, j, Small[K], Large[K]);
				I[K](j, i+Dim[K]) = I[K](i+Dim[K], j);
				if (j >= i) {
					Overlap[K](i, j) = Calc.S_LL(i, j, Large[K]);
					Overlap[K](j, i) = Overlap[K](i, j);
					Overlap[K](i+Dim[K], j+Dim[K]) = Calc.S_SS(i, j, Small[K]);
					Overlap[K](j+Dim[K], i+Dim[K]) = Overlap[K](i+Dim[K], j+Dim[K]);
					I[K](i, j) = Calc.Nuclear_LL(i, j, Large[K]);
					I[K](j, i) = I[K](i, j);
					I[K](i+Dim[K], j+Dim[K]) = Calc.Nuclear_SS(i, j, Small[K]) - 2*Overlap[K](i+Dim[K], j+Dim[K])/Alpha2;
					I[K](j+Dim[K], i+Dim[K]) = I[K](i+Dim[K], j+Dim[K]);
				}
			}
		}

		GeneralizedSelfAdjointEigenSolver<MatrixXd> Magic(I[K], Overlap[K]);
		EigVecs.push_back(Magic.eigenvectors());
		EigVals.push_back(Magic.eigenvalues());
	}
}

void LinearAlgebra::SolveFock(vector<Gspinor> & Large)
{
	I.resize(Large.size());
	Overlap.resize(Large.size());
	Dim.resize(Large.size());

	for (int i = 0; i < Large.size(); i++) {
		Dim[i] = Large[i].lambda.size();
		I[i].resize(Dim[i],Dim[i]);
		Overlap[i] = MatrixXd::Zero(Dim[i],Dim[i]);
	}

	OnePartInt Calc(input, Large);

	for (int K = 0; K < Large.size(); K++) {
		for (int i = 0; i < Dim[K]; i++) {
			for (int j = i; j < Dim[K]; j++) {
				Overlap[K](i, j) = Calc.S_LL(i, j, Large[K]);
				Overlap[K](j, i) = Overlap[K](i, j);
				I[K](i, j) = Calc.Kinetic_NR(i, j, Large[K]) + Calc.Nuclear_NR(i, j, Large[K]);
				I[K](j, i) = I[K](i, j);
			}
		}

		GeneralizedSelfAdjointEigenSolver<MatrixXd> Magic(I[K], Overlap[K]);
		EigVecs.push_back(Magic.eigenvectors());
		EigVals.push_back(Magic.eigenvalues());
	}
}

void LinearAlgebra::IterateHF(vector<Gspinor> & Large)
{
	// Non-relativistic atomic Hartree-Fock calculation.
	// Solves set of coupled generalized eigenvalue problems:
	//
	//	F(l)C(l) = e(l)S(l)C(l)
	//
	// l is an orbital quantum number. Equations for different l are linked through Hartree-Fock potential,
	// which depends on all l.
	// SolveFock SHOULD BE RUN BEFOREHAND TO GET AN INITIAL GUESS !!!
	printf("=============================================================================\n");
	printf("===                  Non-relativistic Hartree-Fock                        ===\n");
	printf("=============================================================================\n");
	printf("%4s %15s %11s %18s %11s\n", "Iter", "Nucl + Kinetic", "Coulomb", "Total Energy", "Сhange");
	cG.resize(Large.size());
	oG.resize(Large.size());

  Dens.cLL.resize(Large.size());
  Dens.oLL.resize(Large.size());

  Dens_test.resize(Large.size());
  G.resize(Large.size());

	for (int K = 0; K < Large.size(); K++) {
    //Dens.cLL[K] = MatrixXd::Zero(Dim[K],Dim[K]);
    //Dens.oLL[K] = MatrixXd::Zero(Dim[K],Dim[K]);

		cG[K].resize(Dim[K],Dim[K]);
		oG[K].resize(Dim[K],Dim[K]);
    
    G[K].cl=vector<vector<double>>(Dim[K], vector<double>(Dim[K], 0.));
    G[K].op=vector<vector<double>>(Dim[K], vector<double>(Dim[K], 0.));

    Dens_test[K].cLL=vector<vector<double>>(Dim[K], vector<double>(Dim[K], 0.));
    Dens_test[K].oLL=vector<vector<double>>(Dim[K], vector<double>(Dim[K], 0.));
	}
	vector<MatrixXd> Fock(Large.size());
	// Set start up Density.
  oEigVecs = EigVecs;
  oEigVals = EigVals;

	make_cDens(1, Large);
	make_oDens(1, Large);

	double p = 0.6;
	double E0 = 0, E1 = 0, dE = 0; // Zero and first order energies.

	// Primary loop.
	int it = -1; 
	do {
		it++;
		// Old energy.
		dE = -1*(E1 + E0);
		E0 = 0;
		E1 = 0;

		// Get electron repulsion matrices cG and oG.
		for (int K1 = 0; K1 < Large.size(); K1++) {
			for (int K2 = K1; K2 < Large.size(); K2++) {
				if (K1 == K2 && Large[K1].is_open) make_oG(Large, K1);
				else make_cG(Large, K1, K2); 
			}
		}
    // Fill in a remaining half of the matrixes.
		//sym_Matr(cG);
		//sym_Matr(oG);

    sym_Matr(G);

		// Fill lower triangle of Coulomb matrix.
		for (int K = 0; K < Large.size(); K++) {
			Fock[K] = I[K] + cG[K];
			
			GeneralizedSelfAdjointEigenSolver<MatrixXd> Magic(Fock[K], Overlap[K]);
			EigVecs[K] = Magic.eigenvectors();
			EigVals[K] = Magic.eigenvalues();
			
			if (Large[K].is_open) {
				//MatrixXd B = make_cProj(Large, K);
				Fock[K] = Fock[K] + oG[K];
				//Fock[K] = Fock[K] - Overlap[K] * B * Fock[K] - Fock[K] * B * Overlap[K] \
                  + Overlap[K] * B * Fock[K] * B * Overlap[K];
				
				GeneralizedSelfAdjointEigenSolver<MatrixXd> oMagic(Fock[K], Overlap[K]);
        oEigVecs[K] = oMagic.eigenvectors();
				oEigVals[K] = oMagic.eigenvalues();
			}
		}

		// New energy.
    bool is_open = false;
    int dim = 0;
    
    make_cDens(p, Large);
		make_oDens(p, Large);

		for (int K = 0; K < Large.size(); K++) {
      is_open = Large[K].is_open;
      dim = Dim[K];
			for (int i = 0; i < dim; i++) {
				for (int j = i; j < dim; j++) {
          if (i == j) {
            E0 += I[K](i, i)*Dens_test[K].cLL[i][i];
            E1 += G[K].cl[i][i]*Dens_test[K].cLL[i][i];
            if (is_open) E1 += G[K].op[i][i]*Dens_test[K].cLL[i][i];
          } else {
            E0 += 2*I[K](i, j)*Dens_test[K].cLL[i][j];
            E1 += 2*G[K].cl[i][j]*Dens_test[K].cLL[i][j];
            if (is_open) E1 += 2*G[K].op[i][j]*Dens_test[K].cLL[i][j];
          }
				}
			}
      G[K].cl = vector<vector<double>>(dim, vector<double>(dim, 0));
      G[K].op = vector<vector<double>>(dim, vector<double>(dim, 0));
			cG[K] = MatrixXd::Zero(Dim[K],Dim[K]);
			oG[K] = MatrixXd::Zero(Dim[K],Dim[K]);
		}

		E1 *= 0.5;
		dE += E0 + E1;

		//if (fabs(dE/(E0 + E1)) < 0.00001) p = 0.9;
		printf("%3d  %5.8f  %8.8f  %8.8f  %2.5E\n", it, E0, E1, E0+E1, fabs(dE/(E0 + E1)));
	} while (it < max_iter && fabs(dE/(E0+E1)) > Etoller);
	printf("=============================================================================\n");
	printf("===                       Orbital Energies (a.u.)                         ===\n");
	printf("=============================================================================\n");
	for (int K = 0; K < Large.size(); K++) {
		printf("\n");
		for (int i = 0; i < Large[K].occup.size(); i++) {
			if (Large[K].occup[i] == 4*Large[K].l + 2)	printf("  %d%c :  [%-6.8f]\n", i+1+Large[K].l, L_to_symb.at(Large[K].l), EigVals[K][i]);
			else printf("  %d%c :  [%-6.8f]\n", i+1+Large[K].l, L_to_symb.at(Large[K].l), oEigVals[K][i]);
		}
	}
	printf("\n");
}


void LinearAlgebra::IterateHF(vector<Gspinor> & Large, vector<Gspinor> & Small)
{
	// Atomic Dirac-Hartree-Fock calculation.
	// Solves set of coupled generalized eigenvalue problems:
	//
	//	F(kappa)C(kappa) = e(kappa)S(kappa)C(kappa)
	//
	// kappa is a relativistic orbital quantum number. Equations for different kappa
	//  are linked through Hartree-Fock potential, which depends on all kappa.
	// SolveFock SHOULD BE RUN BEFOREHAND TO GET AN INITIAL GUESS !!!
	printf("=============================================================================\n");
	printf("===                        Dirac-Hartree-Fock                             ===\n");
	printf("=============================================================================\n");
	printf("%4s %15s %11s %19s %12s\n", "Iter", "Nucl+Kinetic", "Coulomb", "Total Energy", "Сhange");
	cG.resize(Large.size());
	oG.resize(Large.size());

	Dens.cLL.resize(Large.size());
	Dens.cLS.resize(Large.size());
	Dens.cSS.resize(Large.size());
  Dens.oLL.resize(Large.size());
	Dens.oLS.resize(Large.size());
	Dens.oSS.resize(Large.size());

	for (int K = 0; K < Large.size(); K++) {
		Dens.cLL[K] = MatrixXd::Zero(Dim[K],Dim[K]);
		Dens.cLS[K] = MatrixXd::Zero(Dim[K],Dim[K]);
		Dens.cSS[K] = MatrixXd::Zero(Dim[K],Dim[K]);
    Dens.oLL[K] = MatrixXd::Zero(Dim[K],Dim[K]);
		Dens.oLS[K] = MatrixXd::Zero(Dim[K],Dim[K]);
		Dens.oSS[K] = MatrixXd::Zero(Dim[K],Dim[K]);
		cG[K] = MatrixXd::Zero(2*Dim[K],2*Dim[K]);
		oG[K] = MatrixXd::Zero(2*Dim[K],2*Dim[K]);
	}
	vector<MatrixXd> Fock(Large.size());
	// Set stratup Density. First Dim[K] eigenstates are the Sea states.
  oEigVecs = EigVecs;
  oEigVals = EigVals;

	make_cDens(1, Large, Small);
  make_oDens(1, Large, Small);

	double p = 0.6;
	double E0 = 0, E1 = 0, dE = 0; // Zero and first order energies.
	int it = -1;

	TwoPartInt Get(Large, Small);
	// Primary loop.
	do {
		it++;
		// Old energy.
		dE = -1*(E1 + E0);
		E0 = 0;
		E1 = 0;

		// Get electron repulsion matrices cG and oG.
		for (int K1 = 0; K1 < Large.size(); K1++) {
			for (int K2 = K1; K2 < Large.size(); K2++) {
				if (Large[K1].is_open && (Large[K1].l == Large[K2].l || K1 == K2)) make_oG(Large, Small, K1, K2);
				else make_cG(Large, Small, K1, K2);
			}
		}
    // Fill in lower triangle in Coulomb matrixes.
		REL_sym_Matr(cG);
		REL_sym_Matr(oG);
		
		for (int K = 0; K < Large.size(); K++) {
			Fock[K] = I[K] + cG[K];

			GeneralizedSelfAdjointEigenSolver<MatrixXd> Magic(Fock[K], Overlap[K]);
			EigVecs[K] = Magic.eigenvectors();
			EigVals[K] = Magic.eigenvalues();
      
      if (Large[K].is_open) {
				//MatrixXd B = make_cProj(Large, K);
				Fock[K] = Fock[K] + oG[K];
				//Fock[K] = Fock[K] - Overlap[K] * B * Fock[K] - Fock[K] * B * Overlap[K] \
                  + Overlap[K] * B * Fock[K] * B * Overlap[K];
				
				GeneralizedSelfAdjointEigenSolver<MatrixXd> oMagic(Fock[K], Overlap[K]);
				oEigVals[K] = oMagic.eigenvalues();
        oEigVecs[K] = oMagic.eigenvectors();
			}
		}

		// Set up Density. First Dim[K] eigenstates are the Sea states.
		make_cDens(p, Large, Small);
    make_oDens(p, Large, Small);
		// New Energy.
		bool is_open = false;
		for (int K = 0; K < Large.size(); K++) {
			is_open = Large[K].is_open;
			for (int i = 0; i < Dim[K]; i++) {
				for (int j = i; j < Dim[K]; j++) {
					if (i == j) {
						E0 += I[K](i, i)*Dens.cLL[K](i, i) + I[K](i+Dim[K], i+Dim[K])*Dens.cSS[K](i, i);
						E0 += 2*I[K](i, i+Dim[K])*Dens.cLS[K](i, i);

						E1+= cG[K](i, i)*Dens.cLL[K](i, i);
						E1+= cG[K](i+Dim[K], i+Dim[K])*Dens.cSS[K](i, i);
						E1+= 2*cG[K](i, i+Dim[K])*Dens.cLS[K](i, i);
						
						if (is_open) {
							E1+= oG[K](i, i)*Dens.cLL[K](i, i);
							E1+= oG[K](i+Dim[K], i+Dim[K])*Dens.cSS[K](i, i);
							E1+= 2*oG[K](i, i+Dim[K])*Dens.cLS[K](i, i);
						}
					} else {
						E0 += 2*(I[K](i, j)*Dens.cLL[K](i, j) + I[K](i+Dim[K], j+Dim[K])*Dens.cSS[K](i, j));
						E0 += 2*(I[K](i+Dim[K], j)*Dens.cLS[K](j, i) + I[K](i, j+Dim[K])*Dens.cLS[K](i, j));

						E1+= 2*cG[K](i, j)*Dens.cLL[K](i, j);
						E1+= 2*cG[K](i+Dim[K], j+Dim[K])*Dens.cSS[K](i, j);
						E1+= 2*(cG[K](i, j+Dim[K])*Dens.cLS[K](i, j) + cG[K](j, i+Dim[K])*Dens.cLS[K](j, i));

						if(is_open) {
							E1+= 2*oG[K](i, j)*Dens.cLL[K](i, j);
							E1+= 2*oG[K](i+Dim[K], j+Dim[K])*Dens.cSS[K](i, j);
							E1+= 2*(oG[K](i, j+Dim[K])*Dens.cLS[K](i, j) + oG[K](j, i+Dim[K])*Dens.cLS[K](j, i));							
						}
					}
				}
			}
			cG[K] = MatrixXd::Zero(2*Dim[K],2*Dim[K]);
      oG[K] = MatrixXd::Zero(2*Dim[K],2*Dim[K]);
		}

		E1 *= 0.5;
		dE += E0 + E1;
		printf("%-3d    %-8.8f    %-8.8f    %-8.8f    %-2.5E\n", it, E0, E1, E0+E1, fabs(dE/(E0 + E1)));
	} while (it < max_iter && fabs(dE/(E0+E1)) > Etoller);
	// Output final energies to the screen.
	printf("=============================================================================\n");
	printf("===                       Orbital Energies (a.u.)                         ===\n");
	printf("=============================================================================\n");
	string myStr;
  double energy;
	for (int K = 0; K < Large.size(); K++) {
		printf("\n");
		for (int i = 0; i < Large[K].occup.size(); i++) {
      myStr = SymbAng(Large[K].k);
      if (Large[K].occup[i] == 2*abs(Large[K].k)) energy = EigVals[K][i+Dim[K]];
      else energy = oEigVals[K][i+Dim[K]];
			
			printf("  %d%s :  [%6.8f]\n", i+1+Large[K].l, myStr.c_str(), energy);
		}
	}
	printf("\n");
}

void LinearAlgebra::make_cDens(double p, vector<Gspinor> & Large)
{
	// Construct non-relativistic density matrix.
	double DensLLtmp = 0;
	int ClsOcc = 0;
  int dim = 0;
	for (int K = 0; K < Large.size(); K++) {
		ClsOcc = 4 * Large[K].l + 2;
    dim = Dim[K];
		for (int i = 0; i < dim; i++) {
			for (int j = i; j < dim; j++) {
				DensLLtmp = 0;
				for (int N = 0; N < Large[K].occup.size(); N++) {
					if (Large[K].occup[N] == ClsOcc) DensLLtmp += Large[K].occup[N] * EigVecs[K](i, N) * EigVecs[K](j, N);
					else DensLLtmp += Large[K].occup[N] * oEigVecs[K](i, N) * oEigVecs[K](j, N);
				}
				Dens_test[K].cLL[i][j] = (1-p)*Dens_test[K].cLL[i][j] + p*DensLLtmp;
				Dens_test[K].cLL[j][i] = Dens_test[K].cLL[i][j];
        //Dens.cLL[K](i, j) = (1-p)*Dens.cLL[K](i, j) + p*DensLLtmp;
        //Dens.cLL[K](j, i) = (1-p)*Dens.cLL[K](i, j) + p*DensLLtmp;
			}
		}
	}

}

void LinearAlgebra::make_oDens(double p, vector<Gspinor> & Large)
{
	// Construct non-relativistic open shell density matrix.
	double DensLLtmp = 0;
	// Find open shell in K. Only one shell in K1 can be opened!!!! TODO: add support for multiple opened shells.
	int N = -1;
	double scaled_occ = 0;
	int closed_occ = 0;
  int dim = 0;

	for (int K = 0; K < Large.size(); K++) {
		if (!Large[K].is_open) continue;
		closed_occ = 4 * Large[K].l + 2;
		// Locate open shell.
		N = -1;
		for (int n = 0; n < Large[K].occup.size(); n++) {
			if (Large[K].occup[n] != closed_occ) {
				N = n;
				break;
			}
		}
    dim = Dim[K];
		if (N == -1 || Large[K].occup[N] == 0) {
      Dens_test[K].oLL.clear();
      Dens_test[K].oLL.resize(dim, vector<double>(dim, 0.));
      continue;
    }
		// Subtract Large[K].occup[N] as it is already added to cDens.
		scaled_occ = (Large[K].occup[N] - 1)*closed_occ/(closed_occ - 1) - Large[K].occup[N];
		for (int i = 0; i < Dim[K]; i++) {
			for (int j = i; j < Dim[K]; j++) {
				DensLLtmp = scaled_occ * oEigVecs[K](i, N) * oEigVecs[K](j, N);
				
				Dens_test[K].oLL[i][j] = (1-p)*Dens_test[K].oLL[i][j] + p*DensLLtmp;
				Dens_test[K].oLL[j][i] = Dens_test[K].oLL[i][j];
			}
		}
	}
}

void LinearAlgebra::make_cG(vector<Gspinor> & Large, int K1, int K2)
{
	// Function that construct Coulomb repulsion matrix.
	double eAB = 0;
	double angular = 0;
  int dim1 = Dim[K1], dim2 = Dim[K2];

	// Two electron integral calculator.
	TwoPartInt Get(Large);
  // Direct.
	for (int i1 = 0; i1 < dim1; i1++) {
		for (int j1 = i1; j1 < dim1; j1++) {
			for (int i2 = 0; i2 < dim2; i2++) {
				for (int j2 = i2; j2 < dim2; j2++) {
					if (K1 == K2 && i2 < i1) continue;
					eAB = Get.R_k_LL(0, i1, i2, j1, j2, Large[K1], Large[K2], Large[K1], Large[K2]);

					if (i2 == j2) G[K1].cl[i1][j1] += Dens_test[K2].cLL[i2][j2] * eAB;
					else G[K1].cl[i1][j1] += 2*Dens_test[K2].cLL[i2][j2] * eAB;

					if (i1 == i2 && K1 == K2) continue;
					if (i1 == j1) G[K2].cl[i2][j2] += Dens_test[K1].cLL[i1][j1] * eAB;
					else G[K2].cl[i2][j2] += 2*Dens_test[K1].cLL[i1][j1] * eAB;
				}
			}
		}
	}
	// Exchange.
	for (int i1 = 0; i1 < dim1; i1++) {
		for (int j1 = 0; j1 < dim1; j1++) {
			if (i1 > j1 && K1 == K2) continue;
			for (int i2 = 0; i2 < dim2; i2++) {
				for (int j2 = 0; j2 < dim2; j2++) {
					if (i1 > j1 && i2 > j2) continue;
					eAB = 0;
					for (int L = abs(Large[K1].l - Large[K2].l); L <= Large[K1].l + Large[K2].l; L += 2) {
						angular = Get.Wigner3j(Large[K1].l, Large[K2].l, L);
						angular *= angular;
						eAB -= angular*Get.R_k_LL(L, i1, i2, j2, j1, Large[K1], Large[K2], Large[K2], Large[K1]);
					}
					eAB *= 0.5;
					if (i1 <= j1) G[K1].cl[i1][j1] += Dens_test[K2].cLL[i2][j2] * eAB;
					if (i2 <= j2 && K1 != K2) G[K2].cl[i2][j2] += Dens_test[K1].cLL[i1][j1] * eAB;
				}
			}
		}
	}

}

void LinearAlgebra::make_oG(vector<Gspinor> & Large, int K)
{
	// Function that construct Coulomb repulsion matrix.
	double eAB = 0;
	double angular = 0;
	int dim = Dim[K];

	// Two electron integral calculator.
	TwoPartInt Get(Large);

	// Direct.
	for (int i1 = 0; i1 < dim; i1++) {
		for (int j1 = i1; j1 < dim; j1++) {
			for (int i2 = i1; i2 < dim; i2++) {
				for (int j2 = i2; j2 < dim; j2++) {
					eAB = Get.R_k_LL(0, i1, i2, j1, j2, Large[K], Large[K], Large[K], Large[K]);

					if (i2 == j2) {
						G[K].cl[i1][j1] += Dens_test[K].cLL[i2][j2] * eAB;
						G[K].op[i1][j1] += Dens_test[K].oLL[i2][j2] * eAB;
					} else {
						G[K].cl[i1][j1] += 2*Dens_test[K].cLL[i2][j2] * eAB;
						G[K].op[i1][j1] += 2*Dens_test[K].oLL[i2][j2] * eAB;
					}

					if (i1 == i2) continue;
					if (i1 == j1) {
						G[K].cl[i2][j2] += Dens_test[K].cLL[i1][j1] * eAB;
						G[K].op[i2][j2] += Dens_test[K].oLL[i1][j1] * eAB;
					} else {
						G[K].cl[i2][j2] += 2*Dens_test[K].cLL[i1][j1] * eAB;
						G[K].op[i2][j2] += 2*Dens_test[K].oLL[i1][j1] * eAB;
					}
				}
			}
		}
	}
	// Exchange.
	for (int i1 = 0; i1 < dim; i1++) {
		for (int j1 = i1; j1 < dim; j1++) {
			for (int i2 = 0; i2 < dim; i2++) {
				for (int j2 = 0; j2 < dim; j2++) {
					eAB = 0;
					for (int L = abs(Large[K].l - Large[K].l); L <= Large[K].l + Large[K].l; L += 2) {
						angular = Get.Wigner3j(Large[K].l, Large[K].l, L);
						angular *= angular;
						eAB -= angular*Get.R_k_LL(L, i1, i2, j2, j1, Large[K], Large[K], Large[K], Large[K]);
					}
					eAB *= 0.5;

					G[K].cl[i1][j1] += Dens_test[K].cLL[i2][j2] * eAB;
					G[K].op[i1][j1] += Dens_test[K].oLL[i2][j2] * eAB;
				}
			}
		}
	}

}

void LinearAlgebra::make_cDens(double p, vector<Gspinor> & Large, vector<Gspinor> & Small)
{
  // Construct relativistic density matrix.
	double DensLLtmp = 0;
	double DensLStmp = 0;
	double DensSLtmp = 0;
	double DensSStmp = 0;
  int ClsOcc = 0;
  MatrixXd * pEigVects;
	for (int K = 0; K < Large.size(); K++) {
    ClsOcc = 2 * abs(Large[K].k);

		for (int i = 0; i < Dim[K]; i++) {
			for (int j = i; j < Dim[K]; j++) {
				DensLLtmp = 0;
				DensLStmp = 0;
				DensSLtmp = 0;
				DensSStmp = 0;
        // As with NONREL, close shells see open orbital as a closed one with scaled occupancy.
        for (int N = 0; N < Large[K].occup.size(); N++) {
          if (Large[K].occup[N] == ClsOcc) pEigVects = &EigVecs[K];
          else pEigVects = &oEigVecs[K];
          
          DensLLtmp += Large[K].occup[N] * pEigVects->operator()(i, N+Dim[K]) * pEigVects->operator()(j, N+Dim[K]);
          DensLStmp += Large[K].occup[N] * pEigVects->operator()(i, N+Dim[K]) * pEigVects->operator()(j+Dim[K], N+Dim[K]);
          DensSLtmp += Large[K].occup[N] * pEigVects->operator()(i+Dim[K], N+Dim[K]) * pEigVects->operator()(j, N+Dim[K]);
          DensSStmp += Large[K].occup[N] * pEigVects->operator()(i+Dim[K], N+Dim[K]) * pEigVects->operator()(j+Dim[K], N+Dim[K]);
				}
				Dens.cLL[K](i, j) = (1-p)*Dens.cLL[K](i, j) + p*DensLLtmp;
				Dens.cLL[K](j, i) = Dens.cLL[K](i, j);

				Dens.cSS[K](i, j) = (1-p)*Dens.cSS[K](i, j) + p*DensSStmp;
				Dens.cSS[K](j, i) = Dens.cSS[K](i, j);

				Dens.cLS[K](i, j) = (1-p)*Dens.cLS[K](i, j) + p*DensLStmp;
				Dens.cLS[K](j, i) = (1-p)*Dens.cLS[K](j, i) + p*DensSLtmp;
			}
		}
	}
}

void LinearAlgebra::make_oDens(double p, vector<Gspinor> & Large, vector<Gspinor> & Small)
{
	// Construct relativistic open-shell density matrix.
	double DensLLtmp = 0;
	double DensLStmp = 0;
	double DensSLtmp = 0;
	double DensSStmp = 0;
	// Find open shell in K. Only one shell in K1 can be opened!!!! TODO: add support for multiple opened shells.
	int N = -1;
	double scaled_occ = 0;
	int closed_occ = 0;

	for (int K = 0; K < Large.size(); K++) {
		if (!Large[K].is_open) continue;
		closed_occ = 4 * Large[K].l + 2;
		// Locate open shell.
		N = -1;
		for (int n = 0; n < Large[K].occup.size(); n++) {
			if (Large[K].occup[n] != 2*abs(Large[K].k)) {
				N = n;
				break;
			}
		}
		if (N == -1 || Large[K].occup[N] == 0) {
      Dens.oLL[K] = MatrixXd::Zero(Dim[K], Dim[K]);
      Dens.oSS[K] = MatrixXd::Zero(Dim[K], Dim[K]);
      Dens.oLS[K] = MatrixXd::Zero(Dim[K], Dim[K]);

      continue;
    }
		
    scaled_occ = closed_occ*0.5*Large[K].occup[N]/abs(Large[K].k);

		scaled_occ = 2*abs(Large[K].k)*(scaled_occ - 1)/(closed_occ - 1) - Large[K].occup[N];
		for (int i = 0; i < Dim[K]; i++) {
			for (int j = i; j < Dim[K]; j++) {
				DensLLtmp = scaled_occ * oEigVecs[K](i, N+Dim[K]) * oEigVecs[K](j, N+Dim[K]);
				DensLStmp = scaled_occ * oEigVecs[K](i, N+Dim[K]) * oEigVecs[K](j+Dim[K], N+Dim[K]);
				DensSLtmp = scaled_occ * oEigVecs[K](i+Dim[K], N+Dim[K]) * oEigVecs[K](j, N+Dim[K]);
				DensSStmp = scaled_occ * oEigVecs[K](i+Dim[K], N+Dim[K]) * oEigVecs[K](j+Dim[K], N+Dim[K]);

				Dens.oLL[K](i, j) = (1-p)*Dens.oLL[K](i, j) + p*DensLLtmp;
				Dens.oLL[K](j, i) = Dens.oLL[K](i, j);

				Dens.oSS[K](i, j) = (1-p)*Dens.oSS[K](i, j) + p*DensSStmp;
				Dens.oSS[K](j, i) = Dens.oSS[K](i, j);

				Dens.oLS[K](i, j) = (1-p)*Dens.oLS[K](i, j) + p*DensLStmp;
				Dens.oLS[K](j, i) = (1-p)*Dens.oLS[K](j, i) + p*DensSLtmp;
			}
		}
	}

}


void LinearAlgebra::make_cG(vector<Gspinor> & Large, vector<Gspinor> & Small, int K1, int K2)
{
  double eABll = 0, eABls = 0, eABsl = 0, eABss = 0;
	double angular = 0;
  int dim1 = Dim[K1], dim2 = Dim[K2];

  TwoPartInt Get(Large, Small);

  // Direct.
  for (int i1 = 0; i1 < dim1; i1++) {
    for (int j1 = i1; j1 < dim1; j1++) {
      for (int i2 = 0; i2 < dim2; i2++) {
        for (int j2 = i2; j2 < dim2; j2++) {
          if (K1 == K2 && i2 < i1) continue;
          eABll = Get.R_k_LL(0, i1, i2, j1, j2, Large[K1], Large[K2], Large[K1], Large[K2]);
          eABss = Get.R_k_SS(0, i1, i2, j1, j2, Small[K1], Small[K2], Small[K1], Small[K2]);
          eABls = Get.R_k_LS(0, i1, i2, j1, j2, Large[K1], Small[K2], Large[K1], Small[K2]);
          eABsl = Get.R_k_LS(0, i2, i1, j2, j1, Large[K2], Small[K1], Large[K2], Small[K1]);

          if (i2 == j2) {
            cG[K1](i1, j1) += Dens.cLL[K2](i2, j2) * eABll + Dens.cSS[K2](i2, j2) * eABls;
            cG[K1](i1+dim1, j1+dim1) += Dens.cSS[K2](i2, j2) * eABss + Dens.cLL[K2](i2, j2) * eABsl;
          } else {
            cG[K1](i1, j1) += 2*(Dens.cLL[K2](i2, j2) * eABll + Dens.cSS[K2](i2, j2) * eABls);
            cG[K1](i1+dim1, j1+dim1) += 2 * (Dens.cSS[K2](i2, j2) * eABss + Dens.cLL[K2](i2, j2) * eABsl);
          }

          if (i1 == i2 && K1 == K2) continue;
          if (i1 == j1) {
            cG[K2](i2, j2) += Dens.cLL[K1](i1, j1) * eABll + Dens.cSS[K1](i1,j1) * eABsl;
            cG[K2](i2+dim2, j2+dim2) += Dens.cSS[K1](i1, j1) * eABss + Dens.cLL[K1](i1, j1) * eABls;
          } else {
            cG[K2](i2, j2) += 2 * (Dens.cLL[K1](i1, j1) * eABll + Dens.cSS[K1](i1,j1) * eABsl);
            cG[K2](i2+dim2, j2+dim2) += 2 * (Dens.cSS[K1](i1, j1) * eABss + Dens.cLL[K1](i1, j1) * eABls);
          }
        }
      }
    }
  }
  // Exchange.				
  for (int i1 = 0; i1 < dim1; i1++) {
    for (int j1 = 0; j1 < dim1; j1++) {
      //R(LSLS) & R(SLSL) - less symmetries.
      for (int i2 = 0; i2 < dim2; i2++) {
        for (int j2 = 0; j2 < dim2; j2++) {
          eABls = 0;
          for (int L = (int)fabs(Large[K1].j - Large[K2].j); L <= (int)(Large[K1].j + Large[K2].j); L++) {
            if ( (L + Large[K1].l + Large[K2].l) % 2 != 0 ) continue;
            angular = Get.Wigner3j(Large[K1].j, Large[K2].j, L);
            angular *= angular;
            eABls -= angular*Get.R_k_LS(L, i1, i2, j2, j1, Large[K1], Small[K2], Large[K2], Small[K1]);
          }
          cG[K1](i1, j1+dim1) += Dens.cLS[K2](j2, i2) * eABls;
          if (K1 != K2) cG[K2](j2, i2+dim2) += Dens.cLS[K1](i1, j1) * eABls;
        }
      }
      //R(LLLL) & R(SSSS) - non-relativistic symmetries.
      if (i1 > j1 && K1 == K2) continue;
      for (int i2 = 0; i2 < dim2; i2++) {
        for (int j2 = 0; j2 < dim2; j2++) {
          if (i1 > j1 && i2 > j2) continue;
          eABll = 0;
          eABss = 0;
          for (int L = (int)fabs(Large[K1].j - Large[K2].j); L <= (int)(Large[K1].j + Large[K2].j); L++) {
            if ( (L + Large[K1].l + Large[K2].l) % 2 != 0 ) continue;
            angular = Get.Wigner3j(Large[K1].j, Large[K2].j, L);
            angular *= angular;
            eABll -= angular*Get.R_k_LL(L, i1, i2, j2, j1, Large[K1], Large[K2], Large[K2], Large[K1]);
            eABss -= angular*Get.R_k_SS(L, i1, i2, j2, j1, Small[K1], Small[K2], Small[K2], Small[K1]);
          }
          if (i1 <= j1) {
            cG[K1](i1, j1) += Dens.cLL[K2](i2, j2) * eABll;
            cG[K1](i1+dim1, j1+dim1) += Dens.cSS[K2](i2, j2) * eABss;
          }								
          if (i2 <= j2 && K1 != K2) {
            cG[K2](i2, j2) += Dens.cLL[K1](i1, j1) * eABll;
            cG[K2](i2+dim2, j2+dim2) += Dens.cSS[K1](i1, j1) * eABss;
          }
        }
      }
    }
  }

}

void LinearAlgebra::make_oG(vector<Gspinor> & Large, vector<Gspinor> & Small, int K1, int K2)
{
  double eABll = 0, eABls = 0, eABsl = 0, eABss = 0;
	double angular = 0;
  int dim1 = Dim[K1], dim2 = Dim[K2];

  TwoPartInt Get(Large, Small);

  // Direct.
  for (int i1 = 0; i1 < dim1; i1++) {
    for (int j1 = i1; j1 < dim1; j1++) {
      for (int i2 = 0; i2 < dim2; i2++) {
        for (int j2 = i2; j2 < dim2; j2++) {
          if (K1 == K2 && i2 < i1) continue;
          eABll = Get.R_k_LL(0, i1, i2, j1, j2, Large[K1], Large[K2], Large[K1], Large[K2]);
          eABss = Get.R_k_SS(0, i1, i2, j1, j2, Small[K1], Small[K2], Small[K1], Small[K2]);
          eABls = Get.R_k_LS(0, i1, i2, j1, j2, Large[K1], Small[K2], Large[K1], Small[K2]);
          eABsl = Get.R_k_LS(0, i2, i1, j2, j1, Large[K2], Small[K1], Large[K2], Small[K1]);

          if (i2 == j2) {
            cG[K1](i1, j1) += Dens.cLL[K2](i2, j2) * eABll + Dens.cSS[K2](i2, j2) * eABls;
            cG[K1](i1+dim1, j1+dim1) += Dens.cSS[K2](i2, j2) * eABss + Dens.cLL[K2](i2, j2) * eABsl;
						oG[K1](i1, j1) += Dens.oLL[K2](i2, j2) * eABll + Dens.oSS[K2](i2, j2) * eABls;
            oG[K1](i1+dim1, j1+dim1) += Dens.oSS[K2](i2, j2) * eABss + Dens.oLL[K2](i2, j2) * eABsl;
          } else {
            cG[K1](i1, j1) += 2*(Dens.cLL[K2](i2, j2) * eABll + Dens.cSS[K2](i2, j2) * eABls);
            cG[K1](i1+dim1, j1+dim1) += 2 * (Dens.cSS[K2](i2, j2) * eABss + Dens.cLL[K2](i2, j2) * eABsl);
						oG[K1](i1, j1) += 2*(Dens.oLL[K2](i2, j2) * eABll + Dens.oSS[K2](i2, j2) * eABls);
            oG[K1](i1+dim1, j1+dim1) += 2 * (Dens.oSS[K2](i2, j2) * eABss + Dens.oLL[K2](i2, j2) * eABsl);
          }

          if (i1 == i2 && K1 == K2) continue;
          if (i1 == j1) {
            cG[K2](i2, j2) += Dens.cLL[K1](i1, j1) * eABll + Dens.cSS[K1](i1,j1) * eABsl;
            cG[K2](i2+dim2, j2+dim2) += Dens.cSS[K1](i1, j1) * eABss + Dens.cLL[K1](i1, j1) * eABls;
						oG[K2](i2, j2) += Dens.oLL[K1](i1, j1) * eABll + Dens.oSS[K1](i1,j1) * eABsl;
            oG[K2](i2+dim2, j2+dim2) += Dens.oSS[K1](i1, j1) * eABss + Dens.oLL[K1](i1, j1) * eABls;
          } else {
            cG[K2](i2, j2) += 2 * (Dens.cLL[K1](i1, j1) * eABll + Dens.cSS[K1](i1,j1) * eABsl);
            cG[K2](i2+dim2, j2+dim2) += 2 * (Dens.cSS[K1](i1, j1) * eABss + Dens.cLL[K1](i1, j1) * eABls);
						oG[K2](i2, j2) += 2 * (Dens.oLL[K1](i1, j1) * eABll + Dens.oSS[K1](i1,j1) * eABsl);
            oG[K2](i2+dim2, j2+dim2) += 2 * (Dens.oSS[K1](i1, j1) * eABss + Dens.oLL[K1](i1, j1) * eABls);
          }
        }
      }
    }
  }
  // Exchange.				
  for (int i1 = 0; i1 < dim1; i1++) {
    for (int j1 = 0; j1 < dim1; j1++) {
      //R(LSLS) & R(SLSL) - less symmetries.
      for (int i2 = 0; i2 < dim2; i2++) {
        for (int j2 = 0; j2 < dim2; j2++) {
          eABls = 0;
          for (int L = (int)fabs(Large[K1].j - Large[K2].j); L <= (int)(Large[K1].j + Large[K2].j); L++) {
            if ( (L + Large[K1].l + Large[K2].l) % 2 != 0 ) continue;
            angular = Get.Wigner3j(Large[K1].j, Large[K2].j, L);
            angular *= angular;
            eABls -= angular*Get.R_k_LS(L, i1, i2, j2, j1, Large[K1], Small[K2], Large[K2], Small[K1]);
          }
          cG[K1](i1, j1+dim1) += Dens.cLS[K2](j2, i2) * eABls;
					oG[K1](i1, j1+dim1) += Dens.oLS[K2](j2, i2) * eABls;
          if (K1 != K2) {
						cG[K2](j2, i2+dim2) += Dens.cLS[K1](i1, j1) * eABls;
						oG[K2](j2, i2+dim2) += Dens.oLS[K1](i1, j1) * eABls;
					}
        }
      }
      //R(LLLL) & R(SSSS) - non-relativistic symmetries.
      if (i1 > j1 && K1 == K2) continue;
      for (int i2 = 0; i2 < dim2; i2++) {
        for (int j2 = 0; j2 < dim2; j2++) {
          if (i1 > j1 && i2 > j2) continue;
          eABll = 0;
          eABss = 0;
          for (int L = (int)fabs(Large[K1].j - Large[K2].j); L <= (int)(Large[K1].j + Large[K2].j); L++) {
            if ( (L + Large[K1].l + Large[K2].l) % 2 != 0 ) continue;
            angular = Get.Wigner3j(Large[K1].j, Large[K2].j, L);
            angular *= angular;
            eABll -= angular*Get.R_k_LL(L, i1, i2, j2, j1, Large[K1], Large[K2], Large[K2], Large[K1]);
            eABss -= angular*Get.R_k_SS(L, i1, i2, j2, j1, Small[K1], Small[K2], Small[K2], Small[K1]);
          }
          if (i1 <= j1) {
            cG[K1](i1, j1) += Dens.cLL[K2](i2, j2) * eABll;
            cG[K1](i1+dim1, j1+dim1) += Dens.cSS[K2](i2, j2) * eABss;
						oG[K1](i1, j1) += Dens.oLL[K2](i2, j2) * eABll;
            oG[K1](i1+dim1, j1+dim1) += Dens.oSS[K2](i2, j2) * eABss;
          }								
          if (i2 <= j2 && K1 != K2) {
            cG[K2](i2, j2) += Dens.cLL[K1](i1, j1) * eABll;
            cG[K2](i2+dim2, j2+dim2) += Dens.cSS[K1](i1, j1) * eABss;
						oG[K2](i2, j2) += Dens.oLL[K1](i1, j1) * eABll;
            oG[K2](i2+dim2, j2+dim2) += Dens.oSS[K1](i1, j1) * eABss;
          }
        }
      }
    }
  }

}

// Standard atomic notations for Rel orbitals.
string LinearAlgebra::SymbAng(int kappa)
{
	int l = kappa;
	if (kappa < 0) l = -1 - kappa;

	string Result(1, L_to_symb.at(l));
	int jx2 = 2*abs(kappa)-1;

	Result += "(" + to_string(jx2) + "/" + to_string(2)+ ")";

	return Result;
}

void LinearAlgebra::sym_Matr(vector<MatrixXd> & M)
{
	// Most of the matrixes in this code are symmetric, so only
	// M(i, j) , i <= j are evaluated. Function fills the remaining elements.
	int rows = 0;
	int cols = 0;
	for (int K = 0; K < M.size(); K++) {
		rows = M[K].rows();
		cols = M[K].cols();
		for (int i = 0; i < rows; i++) {
			for (int j = i+1; j < cols; j++) {
				M[K](j, i) = M[K](i, j);
			}
		}
	}

}

void LinearAlgebra::sym_Matr(vector<Fast_Coulomb> & M)
{
	// Most of the matrixes in this code are symmetric, so only
	// M(i, j) , i <= j are evaluated. Function fills the remaining elements.
	int dim = 0;
	for (int K = 0; K < M.size(); K++) {
		dim = Dim[K];
		for (int i = 0; i < dim; i++) {
			for (int j = i; j < dim; j++) {
				cG[K](i, j) = M[K].cl[i][j];
        cG[K](j, i) = M[K].cl[i][j];
        M[K].cl[j][i] = M[K].cl[i][j];
        oG[K](i, j) = M[K].op[i][j];
        oG[K](j, i) = M[K].op[i][j];
        M[K].op[j][i] = M[K].op[i][j];
			}
		}
	}

}

void LinearAlgebra::REL_sym_Matr(vector<MatrixXd> & M)
{
	// Does same as relativistic, but keeps track of LL, SS, and LS blocks.
	int dim2 = 0, dim = 0;
	for (int K = 0; K < M.size(); K++) {
		dim2 = M[K].rows();
    dim = dim2/2;
		for (int i = 0; i < dim; i++) {
			for (int j = i+1; j < dim; j++) {
				M[K](j, i) = M[K](i, j);
        M[K](j+dim, i+dim) = M[K](i+dim, j+dim);
        M[K](i+dim, j) = M[K](j, i+dim);
        M[K](j+dim, i) = M[K](i, j+dim);
			}
		}
	}

}

MatrixXd LinearAlgebra::make_cProj(vector<Gspinor> & Large, int K)
{
	MatrixXd Result = MatrixXd::Zero(Dim[K],Dim[K]);
	double DensLLtmp = 0;
	int N = -1;
	int ClsOcc = 4 * Large[K].l + 2;
	if (!Large[K].is_open) return Result;
	// Locate open shell.
	for (int n = 0; n < Large[K].occup.size(); n++) {
		if (Large[K].occup[n] != ClsOcc) {
			N = n;
			break;
		}
	}

	for (int i = 0; i < Dim[K]; i++) {
		for (int j = i; j < Dim[K]; j++) {
			DensLLtmp = 0;
			DensLLtmp = Large[K].occup[N] * EigVecs[K](i, N) * EigVecs[K](j, N);
			Result(i, j) = DensLLtmp;
			Result(j, i) = DensLLtmp;
		}
	}

	return Result;
}


LinearAlgebra::~LinearAlgebra()
{
}