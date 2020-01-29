#include "AuxFunc.h"
#include <algorithm>

AuxFunc::AuxFunc()
{
	// initializa integer and half-integrer Gamma functions
	// Gamma[0] = Gamma(0)
	// Gamma[1] = Gamma(1/2)
	// Gamma[2n] = Gamma(n)
	// Gamma[2n+1] = Gamma(n+1/2)
	gamma.resize(40);
	gamma[0] = 1;
	gamma[1] = sqrt(M_PI);
	gamma[2] = 1;
	for (int n = 1; n < gamma.size() / 2; n++) {
		if (n != 1) gamma[2 * n] = gamma[2 * n - 2] * n;
		gamma[2 * n + 1] = gamma[2 * n - 1] * (2 * n - 1) * 0.5;
	}
}

double AuxFunc::Gamma(double x)
{
	int index = (int)(2*x);
	if (index < gamma.size()) return gamma[index];
	else return 0;
}

double AuxFunc::OddRINT(int n, double lambda, double zeta)
{
	// Odd R integral:
	// int r^n exp(-lambda r^2) erf(-zeta r)
	double x = lambda / zeta / zeta;
	double poly = 1;// Polinomial part
	double T0 = lambda + zeta * zeta;
	double ratio = zeta * zeta / T0;
	double Result = 0.5*zeta / sqrt(T0) / lambda;

	if (n % 2 != 1) return 0; // This formula only works for odd n
	if (n > 21 || n < 1) return 0; // n can take values form 1 to 21 inclusive

	for (int i = 1; i <= 0.5*(n - 1); i++) {
		Result *= 0.5 * ratio / lambda;
	}
	switch (n) {
		case 3:
			poly = 2 + 3 * x;
			break;
		case 5:
			poly = 8 + x*(20 + x * 15);
			break;
		case 7:
			poly = 16 + x*(56 + x*(70 + x * 35));
			Result *= 3;
			break;
		case 9:
			poly = 128 + x*(576 + x*(1008 + x*(840 + x * 315)));
			Result *= 3;
			break;
		case 11:
			poly = 256 + x*(1408 + x*(3168 + x*(3696 + x*(2310 + x * 693))));
			Result *= 15;
			break;
		case 13:
			poly = 1024 + x*(6656 + x*(18304 + x*(27456 + x*(24024 + x*(12012 + x * 3003)))));
			Result *= 45;
			break;
		case 15:
			poly = 2048 + x*(15360 + x*(49920 + x*(91520 + x*(102960 + x*(72072 + x*(30030 + x * 6435))))));
			Result *= 315;
			break;
		case 17:
			poly = 32768 + x*(278528 + x*(1044480 + x*(2263040 +
				x*(3111680 + x*(2800512 + x*(1633632 + x*(583440 + x * 109395)))))));
			Result *= 315;
			break;
		case 19:
			poly = 65536 + x*(6222592 + x*(2646016 + x*(6615040 +
				x*(10749440 + x*(11824384 + x*(8868288 + x*(4434144 + x*(1385670 + x * 230945))))))));
			Result *= 2835;
			break;
		case 21:
			poly = 262144 + x*(2752512 + x*(13074432 + x*(37044224 +
				x*(69457920 + x*(90295296 + x*(82770688 + x*(53209728 +
				x*(23279256 + x*(6466460 + x * 969969)))))))));
			Result *= 14175;
			break;
	}


	Result *= poly;
	return Result;
}

// Natural number power x^n.
double AuxFunc::Power(double x, int n)
{
	double Result = 0;

	if (n == 0) Result = 1.;
	if (n == 1) Result = x;
	else Result = x * this->Power(x, n - 1);

	return Result;
}


double AuxFunc::BetaFunc(int p, int q, double z)
{
	double X = 1.;
	double C_nk = 1.;
	double Result = 1. / (p + 0.5);
	
	if (p < 1 || q < 1) return 0.;
	for (int i = 1; i < q; i++) {
		X *= (-z);
		C_nk *= (q - i);
		C_nk /= i;
		Result += C_nk * X / (p + 0.5 + i);
	}
	Result *= Power(z, p);
	Result *= sqrt(z);//z^(p+0.5)
	
	return Result;
}

// This function evaluates log(Num!/Denom!).
// t is used to calculate factorial summations in Wigner3j funstion.
// Allows to avoid getting huge factorial products.
double AuxFunc::LogFF(float Num, float Denom)
{
	double Result = 0.;

	if (Num > Denom) {
		for (int i = (int)Denom + 1; i <= (int)Num; i++) {
			Result += log((double)i);
		}
	}
	else if (Num < Denom) {
		for (int i = (int)Num + 1; i <= (int)Denom; i++) {
			Result -= log((double)i);
		}
	}
	else Result = 0.;

	return Result;
}

double AuxFunc::Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
{
	double Result = 0.;
	double fracpart, intpart, J = j1 + j2 + j3;

	//some particular cases can be calculated in simpler and faster way.
	if (m1 + m2 + m3 != 0) return 0.;
	if (j1 < 0. || j2 < 0. || j3 < 0.) return 0.;
	fracpart = modf(J, &intpart);// J should be integer.
	if (fracpart != 0) return 0;
	if (m1 == 0 && m2 == 0) {
		// /j1 j2 j3\ - Edmonds, "Angular momentum in qantum mechanics", p 50
		// \0  0  0 /
		fracpart = modf(0.5*J, &intpart);
		if (fracpart != 0) return 0;
		else {
			vector<double> factor(3);
			factor[0] = J / 2. - j1;
			factor[1] = J / 2. - j2;
			factor[2] = J / 2. - j3;

			sort(factor.begin(), factor.end());// Smallest first, largest last.

			for (int i = 0; i < factor.size(); i++) {
				fracpart = modf(factor[i], &intpart);// J should be integer
				if (fracpart != 0) return 0;
			}

			Result += LogFF(2 * factor[0], factor[2]);//taking the largest factor under the square root in denominator (factor[2]!)*(factor[2]!)
			Result += LogFF(2 * factor[1], factor[2]);
			Result += LogFF(2 * factor[2], J + 1);
			Result *= 0.5;// end of square root 

			Result += LogFF(J / 2, factor[1]);
			Result += LogFF(1, factor[0]);//smallest factor goes unpaired with numerator

			fracpart = modf(0.25*J, &intpart);//fracpart = 0 if J/2 is even 
			if (fracpart == 0) Result = exp(Result);
			else Result = -exp(Result);
		}
	}
	else if ((m1 == 0 || m2 == 0 || m3 == 0) && (m1 == 0.5 || m2 == 0.5 || m3 == 0.5)) {
		// /ja   jb   J\ - Brink and Satcher, "Angular momentum", p 138
		// \0.5 -0.5  0 /
		double ja, jb;
		int columns_permutation = 1;//account for sign flip

		if (m1 == 0) {
			J = j1;
			if (m2 == 0.5) {
				ja = j2;
				jb = j3;
			}
			else {
				ja = j3;
				jb = j2;
				columns_permutation = -1;
			}
		}
		else if (m1 == 0.5) {
			ja = j1;
			if (m2 == 0) {
				J = j2;
				jb = j3;
				columns_permutation = -1;
			}
			else {
				J = j3;
				jb = j2;
			}
		}
		else {
			jb = j1;
			if (m2 == 0) {
				J = j2;
				ja = j3;
			}
			else {
				J = j3;
				ja = j2;
				columns_permutation = -1;
			}
		}

		double K;
		fracpart = modf(0.5*(ja + jb + J), &intpart);
		if (fracpart == 0) K = J;
		else K = J + 1;

		vector<double> factor_qsrt(3);//sort factorials under square root
		factor_qsrt[0] = ja + jb - J;
		factor_qsrt[1] = ja + J - jb;
		factor_qsrt[2] = jb + J - ja;
		sort(factor_qsrt.begin(), factor_qsrt.end());

		vector<double> factor(3);//sort factorials under square root
		factor[0] = (ja + jb - K) / 2;
		factor[1] = (ja + K - jb - 1) / 2;
		factor[2] = (jb + J - ja) / 2;
		sort(factor.begin(), factor.end());

		Result += LogFF(factor_qsrt[0], factor[2]);
		Result += LogFF(factor_qsrt[1], factor[2]);
		Result += LogFF(factor_qsrt[2], (J + ja + jb + 1));
		Result *= 0.5;

		Result += LogFF(1, factor[0]);
		Result += LogFF(0.5*(K + ja + jb), factor[1]);

		fracpart = modf(0.25*(K + ja + jb) - 0.5, &intpart);
		if (fracpart != 0) columns_permutation *= -1;
		Result = columns_permutation*exp(Result)*2. / sqrt((2 * ja + 1)*(2 * jb + 1));
	}
	else {
		vector<double> Regge(9);
		vector<double> ReggeInt(9);
		Regge[0] = j2 + j3 - j1;//R(11)
		Regge[1] = j1 - j2 + j3;//R(12)
		Regge[2] = j1 + j2 - j3;//R(13)
		Regge[3] = j1 + m1;//R(21)
		Regge[4] = j2 + m2;//R(22)
		Regge[5] = j3 + m3;//R(23)
		Regge[6] = j1 - m1;//R(31)
		Regge[7] = j2 - m2;//R(32)
		Regge[8] = j3 - m3;//R(33)

		// if any Regge symbol is <0 or not integer, 3j symbol = 0
		for (int i = 0; i < Regge.size(); i++) {
			if (Regge[i] < 0.) return 0.;
			fracpart = modf(Regge[i], &intpart);
			if (Regge[i] != intpart) return 0.;
			ReggeInt[i] = (int)intpart;
		}

		//using 8.3.29 of Varshalovich.
		Result = 0.;
	}

	return Result;
}

AuxFunc::~AuxFunc()
{
}