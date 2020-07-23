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

#include "AuxFunc.h"
#include "CardIn.h"

using namespace std;

class TwoPartInt : AuxFunc
{
    // Class to calculate two electron atomic integrals.
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