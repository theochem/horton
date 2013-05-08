// Horton is a Density Functional Theory program.
// Copyright (C) 2011-2013 Toon Verstraelen <Toon.Verstraelen@UGent.be>
//
// This file is part of Horton.
//
// Horton is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// Horton is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#include <stdexcept>


long fill_cartesian_polynomials(double* output, long lmax) {
    // Shell l=0
    if (lmax <= 0) return -1;

    // Shell l=1
    if (lmax <= 1) return 0;

    // Shell l=2
    output[3] = output[0]*output[0];
    output[4] = output[0]*output[1];
    output[5] = output[0]*output[2];
    output[6] = output[1]*output[1];
    output[7] = output[1]*output[2];
    output[8] = output[2]*output[2];
    if (lmax <= 2) return 3;

    // Shell l=3
    output[9] = output[0]*output[3];
    output[10] = output[0]*output[4];
    output[11] = output[0]*output[5];
    output[12] = output[0]*output[6];
    output[13] = output[0]*output[7];
    output[14] = output[0]*output[8];
    output[15] = output[1]*output[6];
    output[16] = output[1]*output[7];
    output[17] = output[1]*output[8];
    output[18] = output[2]*output[8];
    if (lmax <= 3) return 9;

    // Shell l=4
    output[19] = output[0]*output[9];
    output[20] = output[0]*output[10];
    output[21] = output[0]*output[11];
    output[22] = output[0]*output[12];
    output[23] = output[0]*output[13];
    output[24] = output[0]*output[14];
    output[25] = output[0]*output[15];
    output[26] = output[0]*output[16];
    output[27] = output[0]*output[17];
    output[28] = output[0]*output[18];
    output[29] = output[1]*output[15];
    output[30] = output[1]*output[16];
    output[31] = output[1]*output[17];
    output[32] = output[1]*output[18];
    output[33] = output[2]*output[18];
    if (lmax <= 4) return 19;

    // Shell l=5
    output[34] = output[0]*output[19];
    output[35] = output[0]*output[20];
    output[36] = output[0]*output[21];
    output[37] = output[0]*output[22];
    output[38] = output[0]*output[23];
    output[39] = output[0]*output[24];
    output[40] = output[0]*output[25];
    output[41] = output[0]*output[26];
    output[42] = output[0]*output[27];
    output[43] = output[0]*output[28];
    output[44] = output[0]*output[29];
    output[45] = output[0]*output[30];
    output[46] = output[0]*output[31];
    output[47] = output[0]*output[32];
    output[48] = output[0]*output[33];
    output[49] = output[1]*output[29];
    output[50] = output[1]*output[30];
    output[51] = output[1]*output[31];
    output[52] = output[1]*output[32];
    output[53] = output[1]*output[33];
    output[54] = output[2]*output[33];
    if (lmax <= 5) return 34;

    // Shell l=6
    output[55] = output[0]*output[34];
    output[56] = output[0]*output[35];
    output[57] = output[0]*output[36];
    output[58] = output[0]*output[37];
    output[59] = output[0]*output[38];
    output[60] = output[0]*output[39];
    output[61] = output[0]*output[40];
    output[62] = output[0]*output[41];
    output[63] = output[0]*output[42];
    output[64] = output[0]*output[43];
    output[65] = output[0]*output[44];
    output[66] = output[0]*output[45];
    output[67] = output[0]*output[46];
    output[68] = output[0]*output[47];
    output[69] = output[0]*output[48];
    output[70] = output[0]*output[49];
    output[71] = output[0]*output[50];
    output[72] = output[0]*output[51];
    output[73] = output[0]*output[52];
    output[74] = output[0]*output[53];
    output[75] = output[0]*output[54];
    output[76] = output[1]*output[49];
    output[77] = output[1]*output[50];
    output[78] = output[1]*output[51];
    output[79] = output[1]*output[52];
    output[80] = output[1]*output[53];
    output[81] = output[1]*output[54];
    output[82] = output[2]*output[54];
    if (lmax <= 6) return 55;

    // Shell l=7
    output[83] = output[0]*output[55];
    output[84] = output[0]*output[56];
    output[85] = output[0]*output[57];
    output[86] = output[0]*output[58];
    output[87] = output[0]*output[59];
    output[88] = output[0]*output[60];
    output[89] = output[0]*output[61];
    output[90] = output[0]*output[62];
    output[91] = output[0]*output[63];
    output[92] = output[0]*output[64];
    output[93] = output[0]*output[65];
    output[94] = output[0]*output[66];
    output[95] = output[0]*output[67];
    output[96] = output[0]*output[68];
    output[97] = output[0]*output[69];
    output[98] = output[0]*output[70];
    output[99] = output[0]*output[71];
    output[100] = output[0]*output[72];
    output[101] = output[0]*output[73];
    output[102] = output[0]*output[74];
    output[103] = output[0]*output[75];
    output[104] = output[0]*output[76];
    output[105] = output[0]*output[77];
    output[106] = output[0]*output[78];
    output[107] = output[0]*output[79];
    output[108] = output[0]*output[80];
    output[109] = output[0]*output[81];
    output[110] = output[0]*output[82];
    output[111] = output[1]*output[76];
    output[112] = output[1]*output[77];
    output[113] = output[1]*output[78];
    output[114] = output[1]*output[79];
    output[115] = output[1]*output[80];
    output[116] = output[1]*output[81];
    output[117] = output[1]*output[82];
    output[118] = output[2]*output[82];
    if (lmax <= 7) return 83;

    throw std::domain_error("Encountered lmax > 7.");
}
