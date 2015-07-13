// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2016 The HORTON Development Team
//
// This file is part of HORTON.
//
// HORTON is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HORTON is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

// UPDATELIBDOCTITLE: Evaluation of integrals of Gaussian basis functions

#ifndef HORTON_GBASIS_INTS_H
#define HORTON_GBASIS_INTS_H

#include "libint2.h"
#include "horton/gbasis/calc.h"
#include "horton/gbasis/iter_pow.h"


class GB2Integral : public GBCalculator {
    protected:
        long shell_type0, shell_type1;
        const double *r0, *r1;
        IterPow2 i2p;
    public:
        GB2Integral(long max_shell_type);
        void reset(long shell_type0, long shell_type1, const double* r0, const double* r1);
        virtual void add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1) = 0;
        void cart_to_pure();
        const long get_shell_type0() const {return shell_type0;}
        const long get_shell_type1() const {return shell_type1;}
};

/** @brief
 Compute the overlap integrals in a Gaussian orbital basis.
 */
class GB2OverlapIntegral: public GB2Integral {
    public:
        GB2OverlapIntegral(long max_shell_type) : GB2Integral(max_shell_type) {};
        virtual void add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1);
};

/** @brief
 Compute the kinetic integrals in a Gaussian orbital basis.
 */
class GB2KineticIntegral: public GB2Integral {
    public:
        GB2KineticIntegral(long max_shell_type) : GB2Integral(max_shell_type) {};
        virtual void add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1);
};

/** @brief
 Compute the nuclear attraction integrals in a Gaussian orbital basis.
 */
class GB2NuclearAttractionIntegral: public GB2Integral {
    private:
        double* charges;
        double* centers;
        long ncharge;

        double* work_g0;
        double* work_g1;
        double* work_g2;
        double* work_boys;
    public:
        GB2NuclearAttractionIntegral(long max_shell_type, double* charges, double* centers, long ncharge);
        ~GB2NuclearAttractionIntegral();
        virtual void add(double coeff, double alpha0, double alpha1, const double* scales0, const double* scales1);
};

/** @brief
        Compute the (multipole) moment integrals in a Gaussian orbital basis.
        < gto_a | (x - C_x)^l (y - C_y)^m (z - C_z)^n | gto_b >.
 */
class GB2MomentIntegral: public GB2Integral {
    private:
        long* xyz;          //!< Powers for x, y and z of the multipole moment.
        double* center;     //!< The origin w.r.t. to which the multipole moment is computed.

    public:
        /** @brief
                Initialize Moment integral calculator

            @param max_shell_type
                The highest angular momentum index suported

            @param xyz
                The powers of x,y,z in the integrals (l, m, n).

            @param center
                The center [C_x, C_y, C_z] around which the moment integrals arecomputed
        */
        GB2MomentIntegral(long max_shell_type, long* xyz, double* center);

        /** @brief
                Add integrals for a pair of primite shells to the current contraction.

            @param coeff
                The contraction coefficient for the current primitive.

            @param alpha0
                The exponent of the primitive shell 0.

            @param alpha1
                The exponent of the primitive shell 1.

            @param scales0
                The normalization constants for the basis functions in primitive shell 0.

            @param scales1
                The normalization constants for the basis functions in primitive shell 1.
          */
        virtual void add(double coeff, double alpha0, double alpha1,
                         const double* scales0, const double* scales1);
};


class GB4Integral : public GBCalculator {
    protected:
        long shell_type0, shell_type1, shell_type2, shell_type3;
        const double *r0, *r1, *r2, *r3;
    public:
        GB4Integral(long max_shell_type);
        virtual void reset(long shell_type0, long shell_type1, long shell_type2, long shell_type3, const double* r0, const double* r1, const double* r2, const double* r3);
        virtual void add(double coeff, double alpha0, double alpha1, double alpha2, double alpha3, const double* scales0, const double* scales1, const double* scales2, const double* scales3) = 0;
        void cart_to_pure();

        const long get_shell_type0() const {return shell_type0;}
        const long get_shell_type1() const {return shell_type1;}
        const long get_shell_type2() const {return shell_type2;}
        const long get_shell_type3() const {return shell_type3;}
};


typedef struct {
    unsigned int am;
    const double* r;
    double alpha;
} libint_arg_t;

class GB4LibInt : public GB4Integral {
    private:
        Libint_eri_t erieval;
        libint_arg_t libint_args[4];
        long order[4];
        double ab[3], cd[3];
        double ab2, cd2;
    public:
        GB4LibInt(long max_shell_type);
        ~GB4LibInt();
        virtual void reset(long shell_type0, long shell_type1, long shell_type2, long shell_type3, const double* r0, const double* r1, const double* r2, const double* r3);
        virtual void add(double coeff, double alpha0, double alpha1, double alpha2, double alpha3, const double* scales0, const double* scales1, const double* scales2, const double* scales3);

        virtual void laplace_of_potential(double rho, double t, long mmax, double* output) = 0;
};


class GB4ElectronRepulsionIntegralLibInt : public GB4LibInt {
    public:
        GB4ElectronRepulsionIntegralLibInt(long max_shell_type) : GB4LibInt(max_shell_type) {};
        virtual void laplace_of_potential(double rho, double t, long mmax, double* output);
};


#endif
