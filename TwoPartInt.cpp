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
#include "TwoPartInt.h"

TwoPartInt::TwoPartInt(vector<Gspinor> & Large)
{
	// 3j symbol /l_a l_b k\ , max(l_a) = max(l_b) = 6, max(k) = 12
	//           \ 0  0   0/
    int lmax = 6;
    int kmax = 12;
    Stored3j.clear();
    Stored3j.resize(lmax+1);
    for (auto& s: Stored3j) {
        s.resize(lmax+1);
    }
    for (int i = 0; i <= lmax; i++) {
        for (int j = 0; j <= lmax; j++) {
            for (int k = abs(i-j); k <= i+j; k+=2 ) {
                Stored3j[i][j].push_back(AuxFunc::Wigner3j(i, j, k, 0, 0, 0));
            }
        }
    }
}

TwoPartInt::TwoPartInt(vector<Gspinor> & Large, vector<Gspinor> & Small)
{
	// 3j symbol /j_a  j_b k\ , max(j_a) = max(j_b) = 13/2, max(k) = 13
	//           \0.5 -0.5 0/
    int jmax = 6;
    int kmax = 12;
    Stored3j.clear();
    Stored3j.resize(jmax+1);
    for (auto& s: Stored3j) {
        s.resize(jmax+1);
    }
    float ja = 0.5;
    float jb = 0.5;
    for (int i = 0; i <= jmax; i++) {
        ja = i + 0.5;
        for (int j = 0; j <= jmax; j++) {
            jb = j + 0.5;
            for (int k = abs(i-j); k <= i+j+1; k++) {
                Stored3j[i][j].push_back(AuxFunc::Wigner3j(ja, jb, k, 0.5, -0.5, 0));
            }
        }
    }
}

double TwoPartInt::Wigner3j(float j1, float j2, int j3) 
{
    if ( (int)(j1 + j2 + j3) % 2 != 0 && (int)(2*j1) % 2 == 0) return 0;
    if ( j1 + j2 < j3 || j1 + j3 < j2 || j2 + j3 < j1) return 0;
    if ( fabs(j1 - j2) > j3 || fabs(j1 - j3) > j2 || fabs(j2-j3) > j1) return 0;
    int i = (int)j1;
    int j = (int)j2;
    int k = (j3 - (int)fabs(j1-j2));
    if ((int)(2*j1) % 2 == 0) k /= 2;

    return Stored3j[i][j][k]; 
}

// The Good: T = L, T' = L. Only Good in non-relativistic calculations :-)
double TwoPartInt::R_k_LL(int k, int a, int b, int c, int d, 
                          Gspinor & A, Gspinor & B, Gspinor & C, Gspinor & D)
{
    int nAC = A.l + C.l + 2;
    int nBD = B.l + D.l + 2;
    if (nAC + k % 2 == 1 || nBD + k % 2 == 1) return 0.;
    double lamAC = A.lambda[a] + C.lambda[c];
    double lamBD = B.lambda[b] + D.lambda[d];
   
    double Result = A.Norm[a]*B.Norm[b]*C.Norm[c]*D.Norm[d];
    Result *= RadInt(k, nAC, nBD, lamAC, lamBD);

    return Result;
}

// The Bad: T = L, T' = S.
double TwoPartInt::R_k_LS(int k, int a, int b, int c, int d,
                          Gspinor & A, Gspinor & B, Gspinor & C, Gspinor & D)
{
    if (A.l + C.l + k % 2 == 1 || B.l + D.l + k % 2 == 1) return 0.;

    double lamAC = A.lambda[a] + C.lambda[c];
    double lamBD = B.lambda[b] + D.lambda[d];
    double Result = 0.;

    int nAC = A.l + C.l + 2;
    int nBD = B.l + D.l + 2;

    if (B.k > 0 || D.k > 0) {
        Result = RadInt(k, nAC, nBD, lamAC, lamBD);
        Result *= -2*(B.lambda[b]*(D.k + D.l + 1) + 
                      D.lambda[d]*(B.k + B.l + 1));
        if (B.k > 0 && D.k > 0) {
            nBD -= 2;
            Result += (2*D.l + 1)*(2*B.l + 1)*RadInt(k, nAC, nBD, lamAC, lamBD);
        }
    }
    nBD = B.l + D.l + 4;
    Result += 4*B.lambda[b]*D.lambda[d]*RadInt(k, nAC, nBD, lamAC, lamBD);
    Result *= A.Norm[a]*B.Norm[b]*C.Norm[c]*D.Norm[d];

    return Result;
}

// The Ugly: T = S, T' = S.
double TwoPartInt::R_k_SS(int k, int a, int b, int c, int d,
                          Gspinor & A, Gspinor & B, Gspinor & C, Gspinor & D)
{
    if (A.l + C.l + k % 2 == 1 || B.l + D.l + k % 2 == 1) return 0.;

    double lamAC = A.lambda[a] + C.lambda[c];
    double lamBD = B.lambda[b] + D.lambda[d];
    double Result = 0.;

    int nAC = A.l + C.l + 4;
    int nBD = B.l + D.l + 4;

    Result = 4*B.lambda[b]*D.lambda[d]*RadInt(k, nAC, nBD, lamAC, lamBD);
    double Tmp = 0;
    if (B.k > 0 || D.k > 0) {
        nBD = B.l + D.l + 2;
        Tmp = RadInt(k, nAC, nBD, lamAC, lamBD);
        Tmp *= -2*(B.lambda[b]*(D.k + D.l + 1) + 
                   D.lambda[d]*(B.k + B.l + 1));
        if (B.k > 0 && D.k > 0) {
            nBD = B.l + D.l;
            Tmp += (D.k + D.l + 1)*(B.k + B.l + 1)*RadInt(k, nAC, nBD, lamAC, lamBD);
        }
    }
    Result += Tmp;
    Result *= 4*A.lambda[a]*C.lambda[c];  

    if (A.k > 0 || C.k > 0) {
        nAC = A.l + C.l + 2;
        Tmp = 0.;
        if (B.k > 0 || D.k > 0) {
            nBD = B.l + D.l + 2;
            Tmp = RadInt(k, nAC, nBD, lamAC, lamBD);
            Tmp *= -2*(B.lambda[b]*(D.k + D.l + 1) + 
                    D.lambda[d]*(B.k + B.l + 1));
            if (B.k > 0 && D.k > 0) {
                nBD = B.l + D.l;
                Tmp += (D.k + D.l + 1)*(B.k + B.l + 1)*RadInt(k, nAC, nBD, lamAC, lamBD);
            }
        }
        nBD = B.l + D.l + 4;
        Tmp += 4*B.lambda[b]*D.lambda[d]*RadInt(k, nAC, nBD, lamAC, lamBD);
        Tmp *= -2*(A.lambda[a]*(C.k + C.l + 1) + 
                   C.lambda[c]*(A.k + A.l + 1));
        Result += Tmp;

        if (A.k > 0 && C.k > 0) {
            nAC = A.l + C.l;
            Tmp = 0.;
            if (B.k > 0 || D.k > 0) {
                nBD = B.l + D.l + 2;
                Tmp = RadInt(k, nAC, nBD, lamAC, lamBD);
                Tmp *= -2*(B.lambda[b]*(D.k + D.l + 1) + 
                        D.lambda[d]*(B.k + B.l + 1));
                if (B.k > 0 && D.k > 0) {
                    nBD = B.l + D.l;
                    Tmp += (D.k + D.l + 1)*(B.k + B.l + 1)*RadInt(k, nAC, nBD, lamAC, lamBD);
                }
            }
            nBD = B.l + D.l + 4;
            Tmp += 4*B.lambda[b]*D.lambda[d]*RadInt(k, nAC, nBD, lamAC, lamBD);
            Tmp *= (A.k + A.l + 1)*(C.k + C.l + 1);
            Result += Tmp;
        }       
    }

    Result *= A.Norm[a]*B.Norm[b]*C.Norm[c]*D.Norm[d];
    return Result;
}


/* Radial integral:
   int dr1 dr2 r1^nAC exp(-lambdaAC*r1^2) (r<)^k/(r>)^{k+1} r1^nBD exp(-lambdaBD*r2^2),
   H.M. Quiney
   Handbook of Molecular Physics and Quantum Chemistry,
   Volume 2, Part 6, Chapter 31, pp 696–716.
*/
double TwoPartInt::RadInt(int k, int nAC, int nBD,
                          double lamAC, double lamBD)
{
    double z1 = lamAC/(lamAC+lamBD);
    double z2 = lamBD/(lamAC+lamBD);
    int p1 = (nAC + k)/2;// + 1/2
    int q1 = (nBD - k)/2;
    int p2 = (nAC - k)/2;
    int q2 = (nBD + k)/2;// + 1/2

    double Result = BetaFunc(p1, q1, z1) / Power(lamAC, p1) / Power(lamBD, q1) / sqrt(lamAC) +
                    BetaFunc(q2, p2, z2) / Power(lamAC, p2) / Power(lamBD, q2) / sqrt(lamBD);
    Result *= 0.25*Gamma(0.5*(nAC + nBD + 1));

    return Result;
}

TwoPartInt::~TwoPartInt()
{
}