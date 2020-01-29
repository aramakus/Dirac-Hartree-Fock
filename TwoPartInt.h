#pragma once

#include "AuxFunc.h"
#include "CardIn.h"

using namespace std;

class TwoPartInt : AuxFunc
{
    // Calculates two electron atomic integrals.
public:
    // Pre-calculates non-relativistic 3j symbols /la k lb\
    //                                            \ 0 0 0 /
    TwoPartInt(vector<Gspinor> & Large);
    // Pre-calculates relativistic 3j symbols /ja  k  jb \
    //                                        \0.5 0 -0.5/
    TwoPartInt(vector<Gspinor> & Large, vector<Gspinor> & Small);

    // Radial Coulomb integral of the order 'k':
    // R_k_TT' = int dr1 dr2 M*(T,a,r1)M(T,c,r1) (r<)^k/(r>)^{k+1} M*(T',b,r2)M(T',d,r2)
    // 1) T = L, T' = L
    double R_k_LL(int k, int a, int b, int c, int d,
                  Gspinor & A, Gspinor & B, Gspinor & C, Gspinor & D);
    // 2) T = L, T' = S
    double R_k_LS(int k, int a, int b, int c, int d,
                  Gspinor & A, Gspinor & B, Gspinor & C, Gspinor & D);
    // 3) T = S, T' = S
    double R_k_SS(int k, int a, int b, int c, int d,
                  Gspinor & A, Gspinor & B, Gspinor & C, Gspinor & D);
	// Wigner3j for particular m1,m2,m3.
	// 3j symbol /l_a l_b k\ , max(l_a) = max(l_b) = 6, max(k) = 12
	//           \ 0  0   0/
	// 3j symbol /j_a  j_b k\ , max(j_a) = max(j_b) = 13/2, max(k) = 13
	//           \0.5 -0.5 0/
    double Wigner3j(float j1, float j2, int j3);

    ~TwoPartInt();
private:
    double RadInt(int k, int nAC, int nBD, double lamAC, double lamBD);
    vector<vector<vector<double>>> Stored3j;
};