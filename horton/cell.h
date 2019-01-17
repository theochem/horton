// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2019 The HORTON Development Team
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

// UPDATELIBDOCTITLE: Unit cell code to specify periodic boundary conditions

#ifndef HORTON_CELL_H
#define HORTON_CELL_H

/** @brief
        3D/2D/1D periodic boundary conditions and derived quantities.

    Upon construction, an object of this class acts as a read-only
    representation of the periodic boundary conditions. Reciprocal cell vectors,
    vector lengths and spacings between planes are computed immediately. All
    sorts of manipulations of fractional/Cartesian coordinates are supported.

    Note that this implementation is specific for 3D systems, eventhough lower-
    dimensional periodic boundary conditions are supported. In case of 1D or 2D
    PBC, the cell vectors are internally extended with orthogonal basis vectors
    to guarantee an invertible transformation between cartesian and fractional
    coordinates.
 */
class Cell {
    private:
        double rvecs[9], gvecs[9];
        double rlengths[3], glengths[3];
        double rspacings[3], gspacings[3];
        double volume;
        int nvec;
    public:
        /** @brief
                Construct a Cell object.

            @param _rvecs
                A pointer to 3*nvec doubles that represent the real-space
                vectors in row-major ordering.

            @param _nvec
                The number of cell vectors. This corresponds to the
                dimensionality of the periodic boundary conditions.
        */
        Cell(double* _rvecs, int _nvec);

        /** @brief
                Apply the minimum image convention to a real-space relative
                vector.

            @param delta
                A pointer to 3 doubles with the relative vector. It will be
                modified inplace.

            This is an approximate implementation that sometimes fails
            in very skewed cells. This is a common weakness in most
            implementations. For more details see:
            http://scicomp.stackexchange.com/questions/3107/minimum-image-convention-for-triclinic-unit-cell
        */
        void mic(double* delta) const;

        /** @brief
                Convert a real-space vector to fractional coordinates.

            @param cart
                A pointer to 3 doubles containing the input real-space vector.

            @param frac
                A pointer to 3 doubles in which the output is written.
         */
        void to_frac(double* cart, double* frac) const;

        /** @brief
                Convert a fractional coordinates vector to real-space
                coordinates. This can also be interpreted as making a linear
                combination of real-space cell vectors.

            @param frac
                A pointer to 3 doubles containing the input fractional
                coordinates.

            @param cart
                A pointer to 3 doubles to which the output is written
         */
        void to_cart(double* frac, double* cart) const;

        /** @brief
                Construct a linear combination of reciprocal cell vectors.

            @param coeffs
                A pointer to 3 doubles containing the coefficients for the
                linear combination.

            @param gvec
                A pointer to 3 doubles to which the output is written
         */
        void g_lincomb(double* coeffs, double* gvec) const;

        /** @brief
                Compute a dot-product of a real-space vector with each cell
                vector.

            @param cart
                A pointer to 3 doubles containing the input fractional
                coordinates.

            @param dots
                A pointer to 3 doubles to which the output is written.
         */
        void dot_rvecs(double* cart, double* dots) const;

        /** @brief
                Add a linear combination of cell vectors to delta.

            @param delta
                A pointer to 3 doubles for the real-space vector to which the
                linear combination is added inplace.

            @param coeffs
                A pointer to 3 doubles with the coefficients of the linear
                combination.
         */
        void add_rvec(double* delta, long* coeffs) const;


        /** Returns the dimensionality of the periodic boundary conditions. */
        int get_nvec() const {return nvec;};
        /** Return the volume (or area or length) of the cell. */
        double get_volume() const {return volume;};
        /** Return the spacing between the i-th real-space crystal plane */
        double get_rspacing(int i) const;
        /** Return the spacing between the i-th reciprocal crystal plane */
        double get_gspacing(int i) const;
        /** Return the length of the i-th real-space cell vector */
        double get_rlength(int i) const;
        /** Return the length of the i-th reciprocal cell vector */
        double get_glength(int i) const;

        /* The following block of methods is solely useful for the Python
           wrapper. */

        /** Write a copy of rvecs to _rvecs (3*nvec doubles) */
        void copy_rvecs(double* _rvecs) const;
        /** Write a copy of gvecs to _gvecs (3*nvec doubles) */
        void copy_gvecs(double* _gvecs) const;
        /** Write a copy of rlengths to _rlengths (nvec doubles) */
        void copy_rlengths(double* _rlengths) const;
        /** Write a copy of glengths to _glengths (nvec doubles) */
        void copy_glengths(double* _glengths) const;
        /** Write a copy of rspacings to _rspacings (nvec doubles) */
        void copy_rspacings(double* _rspacings) const;
        /** Write a copy of gspacings to _gspacings (nvec doubles) */
        void copy_gspacings(double* _gspacings) const;

        /** @brief
                Get ranges of periodic images without a cutoff radius.

            @param center
                A pointer to three doubles that specify the center of the cutoff
                sphere in real-space.

            @param rcut
                The cutoff radius.

            @param ranges_begin
                A pointer to nvec longs to which the begin of each range of
                periodic images along a periodic boundary condition is written.

            @param ranges_end
                A pointer to nvec longs to which the end of each range of
                periodic images along a periodic boundary condition is written.
                Then end values are non-inclusive as in Python ranges.

            This function effectively defines a supercell that is guaranteed to
            enclose the cutoff sphere.
         */
        void set_ranges_rcut(double* center, double rcut, long* ranges_begin,
            long* ranges_end) const;
};

/**
    @brief
        A standardized modulo operation that works across all compilers.

    @return
        i % shape if pbc is non-zero. -1 is returned otherwise.
 */
long smart_wrap(long i, long shape, long pbc);

#endif
