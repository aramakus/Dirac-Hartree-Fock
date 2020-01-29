#include "OnePartInt.h"

using namespace Constant;

OnePartInt::OnePartInt(CardIn & Input, vector<Gspinor> & Large, vector<Gspinor> & Small) : input(Input)
{
	int dim = 0;
 
 	zeta = sqrt(1.5)*100000*Bohr / (0.836*pow(input.Amass, 1. / 3.) + 0.57);
	// Normalize basis functions.
	for (int K = 0; K < Large.size(); K++) {
		dim = Large[K].lambda.size();
		Large[K].Norm.resize(dim);
		Small[K].Norm.resize(dim);
		for (int i = 0; i < dim; i++) {
			Large[K].Norm[i] = sqrt(2 * pow(2 * Large[K].lambda[i], Large[K].l + 1.5)
										/ Gamma(Large[K].l + 1.5));
			Small[K].Norm[i] = sqrt(2 * pow(2 * Small[K].lambda[i], Small[K].l + 0.5)
										/ Gamma(Large[K].l + 2.5));
		}
	}

}

OnePartInt::OnePartInt(CardIn & Input, vector<Gspinor> & Large) : input(Input)
{
	int dim = 0;
 
 	zeta = sqrt(1.5)*100000*Bohr / (0.836*pow(input.Amass, 1. / 3.) + 0.57);
	// Normalize basis functions.
	for (int K = 0; K < Large.size(); K++) {
		dim = Large[K].lambda.size();
		Large[K].Norm.resize(dim);
		for (int i = 0; i < dim; i++) {
			Large[K].Norm[i] = sqrt(2 * pow(2 * Large[K].lambda[i], Large[K].l + 1.5)
										/ Gamma(Large[K].l + 1.5));
		}
	}
}

double OnePartInt::S_LL(int a, int b, Gspinor & M)
{
	double a_lmbd = M.lambda[a];
	double b_lmbd = M.lambda[b];

	return pow(2*sqrt(a_lmbd*b_lmbd) / (a_lmbd + b_lmbd), M.l + 1.5);
}

double OnePartInt::S_SS(int a, int b, Gspinor & M)
{
	double a_lmbd = M.lambda[a];
	double b_lmbd = M.lambda[b];
	
	return pow(2*sqrt(a_lmbd*b_lmbd) / (a_lmbd + b_lmbd), M.l + 2.5);
}

double OnePartInt::Kinetic_SL(int a, int b, Gspinor & M_S, Gspinor & M_L)
{
	return M_L.Norm[b]/M_S.Norm[b]*S_SS(a,b,M_L)/Alpha;
}

double OnePartInt::Nuclear_LL(int a, int b, Gspinor & M)
{
	double nrm_ab = M.Norm[a] * M.Norm[b];
	double exp_ab = M.lambda[a] + M.lambda[b];

	return -1 * nrm_ab * input.Z * OddRINT(2 * M.l + 1, exp_ab, zeta);
}

double OnePartInt::Nuclear_SS(int a, int b, Gspinor & M)
{
	double nrm_ab = M.Norm[a] * M.Norm[b];
	double exp_ab = M.lambda[a] + M.lambda[b];
	int l = M.l;

	if (M.k < 0) {
		return -4 * nrm_ab * input.Z * M.lambda[a] * M.lambda[b] * OddRINT(2*l + 3, exp_ab, zeta);
	}		
	else {
		return -1 * nrm_ab * input.Z * ((2*l + 1)*(2*l + 1)*OddRINT(2*l - 1, exp_ab, zeta) -
			    2 * exp_ab * (2*l + 1)*OddRINT(2*l + 1, exp_ab, zeta) + 4 * M.lambda[a] * M.lambda[b] * OddRINT(2*l + 3, exp_ab, zeta));
	}
}

double OnePartInt::Kinetic_NR(int a, int b, Gspinor & M)
{
	double nrm_ab = M.Norm[a] * M.Norm[b];
	double exp_ab = M.lambda[a] + M.lambda[b];
	int l = M.l;
	
	return 0.5 * nrm_ab * (2 * l + 3) * M.lambda[a] * M.lambda[b] * Gamma(l+1.5) / pow(exp_ab, l + 2.5);
}


double OnePartInt::Nuclear_NR(int a, int b, Gspinor & M)
{
	return -0.5 * input.Z * M.Norm[a] * M.Norm[b] * Gamma(M.l) 
	    / pow(M.lambda[a] + M.lambda[b], M.l + 1);
}

OnePartInt::~OnePartInt()
{
}
