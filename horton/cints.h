/*******************************************************************************

  The initial version of this file was taken from the PyQuante quantum chemistry
  program suite, version 1.6.4. PyQuante is written by Richard P. Muller.

  Copyright (c) 2004, Richard P. Muller. All Rights Reserved.

  PyQuante version 1.2 and later is covered by the modified BSD license:

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  - Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

  - Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  - Neither the name of Dr. Muller nor the names of its contributors may be used
    to endorse or promote products derived from this software without specific
    prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  You may contact the author of PyQuante via email at rmuller@sandia.gov.

********************************************************************************

  This file was later modified by Toon Verstraelen for its use in Horton.

********************************************************************************

  The following piece of text avoids that the header of the file gets squashed:
  no_update_headers

********************************************************************************

  The equations herein are based upon: 'Gaussian Expansion Methods for Molecular
  Orbitals.' H. Taketa, S. Huzinaga, and K. O-ohata.; J. Phys. Soc. Japan, 21,
  2313, 1966. [THO paper]

*******************************************************************************/


#ifndef HORTON_CINTS_H
#define HORTON_CINTS_H


double fB(int i, int l1, int l2, double px, double ax, double bx,
      int r, double g);
double Bfunc(int i, int r, double g);

double coulomb_repulsion(double xa, double ya, double za, double norma,
                int la, int ma, int na, double alphaa,
                double xb, double yb, double zb, double normb,
                int lb, int mb, int nb, double alphab,
                double xc, double yc, double zc, double normc,
                int lc, int mc, int nc, double alphac,
                double xd, double yd, double zd, double normd,
                int ld, int md, int nd, double alphad);

double *B_array(int l1, int l2, int l3, int l4, double p, double a,
        double b, double q, double c, double d,
        double g1, double g2, double delta);

double B_term(int i1, int i2, int r1, int r2, int u, int l1, int l2,
             int l3, int l4, double Px, double Ax, double Bx,
             double Qx, double Cx, double Dx, double gamma1,
             double gamma2, double delta);
double kinetic(double alpha1, int l1, int m1, int n1,
              double xa, double ya, double za,
              double alpha2, int l2, int m2, int n2,
              double xb, double yb, double zb);
double overlap(double alpha1, int l1, int m1, int n1,
              double xa, double ya, double za,
              double alpha2, int l2, int m2, int n2,
              double xb, double yb, double zb);
double overlap_1D(int l1, int l2, double PAx,
             double PBx, double gamma);
double nuclear_attraction(double x1, double y1, double z1, double norm1,
                 int l1, int m1, int n1, double alpha1,
                 double x2, double y2, double z2, double norm2,
                 int l2, int m2, int n2, double alpha2,
                 double x3, double y3, double z3);
double A_term(int i, int r, int u, int l1, int l2,
             double PAx, double PBx, double CPx, double gamma);
double *A_array(int l1, int l2, double PA, double PB,
               double CP, double g);

int fact(int n);
int fact2(int n);
double dist2(double x1, double y1, double z1,
            double x2, double y2, double z2);
double dist(double x1, double y1, double z1,
           double x2, double y2, double z2);
double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb);
int binomial(int a, int b);

double Fgamma(double m, double x);
double gamm_inc(double a, double x);

int fact_ratio2(int a, int b);

double product_center_1D(double alphaa, double xa,
             double alphab, double xb);

double three_center_1D(double xi, int ai, double alphai,
                  double xj, int aj, double alphaj,
                  double xk, int ak, double alphak);

/* Routines from Numerical Recipes */
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);


#endif
