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

// UPDATELIBDOCTITLE: Evaluate functions on grids and derive Fock matrices from potentials on grids.

/**
    @file fns.h
    @brief Evaluate functions on grids and derive Fock matrices from potentials on grids.

    Functions on grid points are calculated in three steps, which are controlled by
    methods in the GBasis class in gbasis.h. Also code is shared with the Fock build from
    potentials on grids.

    GB1GridFn functions require only a single loop to compute properties of basis
    functions. See following methods in gbasis.h:

        GOBasis::compute_grid1_exp
        GOBasis::compute_grid1_dm
        GOBasis::compute_grid1_fock
        GOBasis::compute_grid_point1 (called by the previous three)

    GB2GridFn functions requires double loop to compute properties of basis pairs of basis
    functions. See following methods in gbasis.h:

        GOBasis::compute_grid2_dm
        GOBasis::compute_grid_point2 (called by compute_grid2_dm)

    The methods `reset`, `add`, `cart_to_pure` and `compute_point_from_*` or
    `compute_fock_from_pot`, are called from functions in gbasis.cpp (GOBasis). For a
    given grid point...

    1) The method `reset` is called with basic information of the contraction and the
       position of the grid point.

    2) The method `add` is called to evaluate properties (function value and/or
       derivatives) of primitives on the grid point. Code in gbasis.cpp takes care of
       contractions and calls `add` with the proper prefactor `coeff`. The `add` method
       adds results for one primitive in `work_cart`.

    3) `cart_to_pure` to transform the results for a contraction to pure functions when
       needed.

    4) After looping over all the (pairs of) contractions (steps 1 to 3), one of the
       following two happens, still just for one grid point

    4a) The compute_point_from_* method is called to convert the properties of the basis
        functions into the (density) function(s) of interest, making use of either
        expansion coefficients of the orbitals or first-order density matrix coefficients.

    4b) The compute_fock_from_pot method is called to convert a potential on a grid into
        a Fock operator.

    The above work flow is repeated for every grid point.


    The Cartesian polynomials in the Gaussian primitives (and/or their derivatives) are
    computed only once for a given contraction when calling the `reset` method. This is
    done by calling `fill_cartesian_polynomials`, which computes all mononomials in
    alphabetical order up to a given order. It returns an offset, which is the position of
    the first mononomial of the highest order. This offset is used to quickly find the
    desired mononomial for evaluating each primitive. The offset is stored as a class
    attribute and used on the `add` method:

    - There is only only offset attribute in many cases, i.e. when only mononomials of one
      order are needed.

    - When mononomials of different orders are needed, their offsets are computed in the
      `reset` method: offset_l1, offset_h1, offset_l2 and offset_h2. The suffixes h and l
      refer to higher and lower.

    - Keep in mind that the derivative (toward x) of a primitive (with x^k) includes two
      terms (one with a mononomial x^{k-1} and one with x^(k+1}). So in the case of GGA
      and even more so for MGGA, many mononomials are needed and several relevant offsets
      must be stored.

    Given the position of a mononomial, related mononomials (increasing or decreasing one
    or two powers) are "easily" found. Given a mononomial:

    - mono = x**nx * y**ny * z**nz: offset + i
    - mono * x: offset_h1 + i
    - mono / x: offset_l1 + i
    - mono * y: offset_h1 + i + 1 + ny + nz
    - mono / y: offset_l1 + i - ny - nz
    - mono * z: offset_h1 + i + 2 + ny + nz
    - mono / z: offset_l1 + i - 1 - ny - nz
    - mono * (x*x): offset_h2 + i
    - mono / (x*x): offset_l2 + i
    - mono * (x*y): offset_h2 + i + 1 + ny + nz
    - mono / (x*y): offset_l2 + i - ny - nz
    - mono * (x*z): offset_h2 + i + 2 + ny + nz
    - mono / (x*z): offset_l2 + i - 1 - ny - nz
    - mono * (y*y): offset_h2 + i + 3 + 2*ny + 2*nz
    - mono / (y*y): offset_l2 + i + 1 - 2*ny - 2*nz
    - mono * (y*z): offset_h2 + i + 4 + 2*ny + 2*nz
    - mono / (y*z): offset_l2 + i - 2*ny - 2*nz
    - mono * (z*z): offset_h2 + i + 5 + 2*ny + 2*nz
    - mono / (z*z): offset_l2 + i - 1 - 2*ny - 2*nz

    These rules are used on the `add` methods.
  */

#ifndef HORTON_GBASIS_FNS_H_
#define HORTON_GBASIS_FNS_H_

#include "horton/gbasis/calc.h"
#include "horton/gbasis/common.h"
#include "horton/gbasis/iter_pow.h"

/** @brief
      Base class for grid calculators that require only a single loop over all basis
      functions.
  */
class GB1GridFn : public GBCalculator  {
 public:
  /** @brief
        Construct a GB1GridFn object.

      @param max_shell_type
        The maximum shell type in the basis set.

      @param dim_work
        A multiplier for the size of the work array, e.g. when multiple results need to
        be stored, such as the orbitals and its gradient.

      @param dim_output
        The number of results for each grid point, e.g. 3 for a density gradient.
    */
  GB1GridFn(long max_shell_type, long dim_work, long dim_output);

  /** @brief
        Reset calculator for a new contraction.

      @param shell_type0
        Shell type of contraction.

      @param r0
        Center of the contraction. (size=3)

      @param point
        Cartesian coordinates of the grid point. (size=3)
   */
  virtual void reset(long shell_type0, const double* r0, const double* point);

  //! Convert results in work array from Cartesian to pure functions where needed.
  void cart_to_pure();

  //! Shell type of the contraction.
  const long get_shell_type0() const {return shell_type0;}

  //! Multiplier for the size of the work array
  long get_dim_work() {return dim_work;}

  //! Number of results per grid point.
  long get_dim_output() {return dim_output;}

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.

      @param coeff
        Contraction coefficient of current primitive shells.

      @param alpha0
        Exponent of the primitive shell.

      @param scales0
        Normalization prefactors for primitives in the shell.
    */
  virtual void add(double coeff, double alpha0, const double* scales0) = 0;

 protected:
  const long dim_work;    //!< Multiplier for the size of the work array.
  const long dim_output;  //!< Number of results per grid point.
  long shell_type0;       //!< Shell type of the current contraction.
  const double *r0;       //!< Center of the current contraction.
  const double *point;    //!< Grid point at which the fn is evaluated
  IterPow1 i1p;           //!< Iterator over Cartesian powers for given ang mom.
};


/** @brief
      Base class for GB1 grid calculators that use the expansion coefficients of the
      orbitals.
  */
class GB1ExpGridFn : public GB1GridFn  {
 public:
  /** @brief
        Construct a GB1ExpGridFn object.

      @param max_shell_type
        The maximum shell type in the basis set.

      @param nfn
        The number of orbitals (occupied and virtual).

      @param dim_work
        A multiplier for the size of the work array, e.g. when multiple results need to
        be stored, such as an orbital and its gradient.

      @param dim_output
        The number of results, i.e. elements in output and pot arguments, for each grid
        point, e.g. 3 for a density gradient.
    */
  GB1ExpGridFn(long max_shell_type, long nfn, long dim_work, long dim_output)
      : GB1GridFn(max_shell_type, dim_work, dim_output), nfn(nfn) {}

  /** @brief
        Compute (final) results for a given grid point.

      @param work_basis
        Properties of basis functions computed for the current grid point. (Work done by
        add method.) (size=nbasis*dim_work)

      @param coeffs
        The orbital expansion coefficients. (size=nbasis*nfn)

      @param nbasis
        The number of basis functions.

      @param output
        The output array for the current grid point. (size=dim_output)
    */
  virtual void compute_point_from_exp(double* work_basis, double* coeffs, long nbasis,
                                      double* output) = 0;

 protected:
  long nfn;  //!< Number of orbitals (occupied and virtual).
};


/** @brief
      Evaluates a selection of orbitals on a grid.

      Content of work_basis (at one grid point):
        [0] Basis function value.
      Content of the argument 'output' (at one grid point):
        [0-norb] Values of the orbitals.
  */
class GB1ExpGridOrbitalFn : public GB1ExpGridFn  {
 public:
  /** @brief
        Construct a GB1ExpGridOrbitalFn object.

      @param max_shell_type
        The maximum shell type in the basis set.

      @param nfn
        The number of orbitals (occupied and virtual).

      @param iorbs
        An array with orbitals to be computed.

      @param norb
        The number of elements in iorbs.
    */
  GB1ExpGridOrbitalFn(long max_shell_type, long nfn, long* iorbs, long norb)
      : GB1ExpGridFn(max_shell_type, nfn, 1, norb), poly_work{0.0}, offset(0),
        iorbs(iorbs), norb(norb) {}

  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);

  //! Compute (final) results for a given grid point. (See base class for details.)
  virtual void compute_point_from_exp(double* work_basis, double* coeffs,
                                      long nbasis, double* output);

 protected:
  double poly_work[MAX_NCART_CUMUL];  //!< Work array with Cartesian polynomials.
  long offset;  //!< Offset for the polynomials for the density.
  long* iorbs;  //!< Array of indices of orbitals to be evaluated on grid.
  long norb;    //!< The number of elements in iorbs.
};


/** @brief
      Evaluates the gradient of a selection of orbitals on a grid.

      Content of work_basis (at one grid point):
        [0] Basis function derivative toward x.
        [1] Basis function derivative toward y.
        [2] Basis function derivative toward z.
      Content of the argument 'output' (at one grid point):
        [i*norb + 0] derivative of orbital i toward x.
        [i*norb + 1] derivative of orbital i toward y.
        [i*norb + 2] derivative of orbital i toward z.
 */
class GB1ExpGridOrbGradientFn : public GB1ExpGridFn  {
 public:
  /** @brief
        Construct a GB1ExpGridOrbGradientFn object.

      @param max_shell_type
        The maximum shell type in the basis set.

      @param nfn
        The number of orbitals (occupied and virtual).

      @param iorbs
        An array with orbitals to be computed.

      @param norb
        The number of elements in iorbs.
    */
  explicit GB1ExpGridOrbGradientFn(long max_shell_type, long nfn, long* iorbs, long norb)
      : GB1ExpGridFn(max_shell_type, nfn, 3, norb*3), poly_work{0.0}, offset(0),
        offset_l1(0), offset_h1(0), iorbs(iorbs), norb(norb) {}
  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);


  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_exp(double* work_basis, double* coeffs,
                                      long nbasis, double* output);

 protected:
  long* iorbs;     //!< Array of indices of orbitals to be evaluated on grid.
  long norb;       //!< The number of elements in iorbs.
  double poly_work[MAX_NCART_CUMUL_D];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
};


/** @brief
      Base class for GB1 grid calculators that use the first-order density matrix
      coefficients.
  */
class GB1DMGridFn : public GB1GridFn  {
 public:
  //! Construct a GB1GridFn object. (See base class for details.)
  GB1DMGridFn(long max_shell_type, long dim_work, long dim_output)
      : GB1GridFn(max_shell_type, dim_work, dim_output) {}

  /** @brief
        Compute the final result on one grid point.

      @param work_basis
        Properties of basis functions computed for the current grid point. (Work done by
        add method.) (size=nbasis*dim_work)

      @param dm
        The coefficients of the first-order density matrix. (size=nbasis*nbasis)

      @param nbasis
        The number of basis functions.

      @param output
        The output array for the current grid point. (size=dim_output)

      @param epsilon
        A cutoff value used to discard small contributions.

      @param dmmaxrow
        The maximum value of the density matrix on each row. (size=nbasis)
    */
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow) = 0;

  /** @brief
        Add contribution to Fock matrix from one grid point.

      The chain rule is used to transform grid potential data (in one point, see reset
      method) into a Fock matrix contribution.

      @param pot
        The value of the potential at the grid point. This may be multiple values, e.g in
        case of GGA, this contains four elements: the functional derivative of the energy
        w.r.t. the density and the components of the gradient. (size=dim_output)

      @param work_basis
        Properties of the orbital basis in the current grid point, typically the value of
        the basis function and optionally first or second derivatives toward x, y and z,
        all evaluated in `point` (see reset method). (size=nbasis*dim_work)

      @param nbasis
        The number of basis functions.

      @param fock
        The Fock matrix to which the result will be added. (size=nbasis*nbasis)
    */
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock) = 0;
};


/** @brief
      Compute just the electron density on a grid.

    Content of work_basis (at one grid point):
      [0] Basis function value.
    Content of the argument 'output' and the energy derivative in 'pot' (at one grid point):
      [0] The electron density
  */
class GB1DMGridDensityFn : public GB1DMGridFn  {
 public:
  /** @brief
        Construct a GB1DMGridDensityFn object.

      @param max_shell_type
        The maximum shell type in the basis set.
    */
  explicit GB1DMGridDensityFn(long max_shell_type)
      : GB1DMGridFn(max_shell_type, 1, 1), poly_work{0.0}, offset(0) {}

  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);

  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow);

  //! Add contribution to Fock matrix for one grid point. (See base class for details.)
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock);

 private:
  double poly_work[MAX_NCART_CUMUL];  //!< Work array with Cartesian polynomials.
  long offset;  //!< Offset for the polynomials for the density
};


/** @brief
      Compute gradient of the electron density on a grid.

    Content of work_basis (at one grid point):
      [0] Basis function value.
      [1] Basis function derivative toward x.
      [2] Basis function derivative toward y.
      [3] Basis function derivative toward z.
    Content of the argument 'output' and the energy derivative in 'pot' (at one grid point):
      [0] Density derivative toward x.
      [1] Density derivative toward y.
      [2] Density derivative toward z.
  */
class GB1DMGridGradientFn : public GB1DMGridFn  {
 public:
  /** @brief
        Construct a GB1DMGridGradientFn object.

      @param max_shell_type
        The maximum shell type in the basis set.
    */
  explicit GB1DMGridGradientFn(long max_shell_type)
      : GB1DMGridFn(max_shell_type, 4, 3), poly_work{0.0}, offset(0),
        offset_l1(0), offset_h1(0) {}

  /** @brief
        Construct a GB1DMGridGradientFn object.

      @param max_shell_type
        The maximum shell type in the basis set.

      @param dim_output
        The number of results for each grid point. This may be different from 3 for
        derived classes.
   */
  explicit GB1DMGridGradientFn(long max_shell_type, long dim_output)
      : GB1DMGridFn(max_shell_type, 4, dim_output), poly_work{0.0}, offset(0),
        offset_l1(0), offset_h1(0) {}

  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);

  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow);

  //! Add contribution to Fock matrix for one grid point. (See base class for details.)
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock);

 protected:
  double poly_work[MAX_NCART_CUMUL_D];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
};


/** @brief
      Compute density and gradient on a grid.

    Content of work_basis (at one grid point):
      [0] Basis function value.
      [1] Basis function derivative toward x.
      [2] Basis function derivative toward y.
      [3] Basis function derivative toward z.
    Content of the argument 'output' and the energy derivative in 'pot' (at one grid point):
      [0] Density.
      [1] Density derivative toward x.
      [2] Density derivative toward y.
      [3] Density derivative toward z.
  */
class GB1DMGridGGAFn : public GB1DMGridGradientFn  {
 public:
  /** @brief
        Construct a GB1DMGridGGAFn object.

      @param max_shell_type
        The maximum shell type in the basis set.
    */
  explicit GB1DMGridGGAFn(long max_shell_type): GB1DMGridGradientFn(max_shell_type, 4) {}

  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow);

  //! Add contribution to Fock matrix for one grid point. (See base class for details.)
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock);
};


/** @brief
      Compute kinetic energy density on a grid.

    Content of work_basis (at one grid point):
      [0] Basis function derivative toward x.
      [1] Basis function derivative toward y.
      [2] Basis function derivative toward z.
    Content of the argument 'output' and the energy derivative in 'pot' (at one grid point):
      [0] Kinetic energy density.
  */
class GB1DMGridKineticFn : public GB1DMGridFn  {
 public:
  /** @brief
        Construct a GB1DMGridKineticFn object.

      @param max_shell_type
        The maximum shell type in the basis set.
    */
  explicit GB1DMGridKineticFn(long max_shell_type)
      : GB1DMGridFn(max_shell_type, 3, 1), poly_work{0.0}, offset(0), offset_l1(0),
        offset_h1(0) {}

  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);

  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow);

  //! Add contribution to Fock matrix for one grid point. (See base class for details.)
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock);

 private:
  double poly_work[MAX_NCART_CUMUL_D];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density.
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
};


/** @brief
      Compute density Hessian on a grid: xx, xy, xz, yy, yz, zz.

    Content of work_basis (at one grid point):
      [0] Basis function value.
      [1] Basis function derivative toward x.
      [2] Basis function derivative toward y.
      [3] Basis function derivative toward z.
      [4] Basis function derivative toward xx.
      [5] Basis function derivative toward xy.
      [6] Basis function derivative toward xz.
      [7] Basis function derivative toward yy.
      [8] Basis function derivative toward yz.
      [9] Basis function derivative toward zz.
    Content of the argument 'output' and the energy derivative in 'pot' (at one grid point):
      [0] Density derivative toward xx.
      [1] Density derivative toward xy.
      [2] Density derivative toward xz.
      [3] Density derivative toward yy.
      [4] Density derivative toward yz.
      [5] Density derivative toward zz.
  */
class GB1DMGridHessianFn : public GB1DMGridFn  {
 public:
  /** @brief
        Construct a GB1DMGridHessianFn object.

      @param max_shell_type
        The maximum shell type in the basis set.
    */
  explicit GB1DMGridHessianFn(long max_shell_type)
      : GB1DMGridFn(max_shell_type, 10, 6), poly_work{0.0}, offset(0),
        offset_l1(0), offset_h1(0), offset_l2(0), offset_h2(0) {}

  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);

  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow);

  //! Add contribution to Fock matrix for one grid point. (See base class for details.)
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock);

 private:
  double poly_work[MAX_NCART_CUMUL_DD];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density.
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
  long offset_l2;  //!< Lower offset for the polynomials for the hessian.
  long offset_h2;  //!< Higher offset for the polynomials for the hessian.
};


/** @brief
      Compute MGGA properties on a grid: density, gradient, laplacian and kinetic energy
      density.

    Content of work_basis (at one grid point):
      [0] Basis function value.
      [1] Basis function derivative toward x.
      [2] Basis function derivative toward y.
      [3] Basis function derivative toward z.
      [4] Basis function Laplacian.
    Content of the argument 'output' and the energy derivative in 'pot' (at one grid point):
      [0] Density.
      [1] Density derivative toward x.
      [2] Density derivative toward y.
      [3] Density derivative toward z.
      [4] Laplacian of the density.
      [5] Kinetic energy density.
  */
class GB1DMGridMGGAFn : public GB1DMGridFn  {
 public:
  /** @brief
        Construct a GB1DMGridMGGAFn object.

      @param max_shell_type
        The maximum shell type in the basis set.
    */
  explicit GB1DMGridMGGAFn(long max_shell_type)
      : GB1DMGridFn(max_shell_type, 5, 6), poly_work{0.0}, offset(0), offset_l1(0),
        offset_h1(0), offset_l2(0), offset_h2(0) {}

  //! Reset calculator for a new contraction. (See base class for details.)
  virtual void reset(long _shell_type0, const double* _r0, const double* _point);

  /** @brief
        Add contributions to work array for current grid point and given primitive shell.
        (See base class for more details.)
    */
  virtual void add(double coeff, double alpha0, const double* scales0);

  //! Compute the final result on one grid point. (See base class for details.)
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow);

  //! Add contribution to Fock matrix for one grid point. (See base class for details.)
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* fock);

 private:
  double poly_work[MAX_NCART_CUMUL_DD];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density.
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
  long offset_l2;  //!< Lower offset for the polynomials for the hessian.
  long offset_h2;  //!< Higher offset for the polynomials for the hessian.
};


/** @brief
      Base class for grid calculators that require a double loop over all basis functions.
  */
class GB2DMGridFn : public GBCalculator  {
 public:
  //! Construct a GB2DMGridFn object. (See base class for details.)
  explicit GB2DMGridFn(long max_shell_type);

  /** @brief
        Reset calculator for a new pair of contractions.

      @param shell_type0
        Shell type of first contraction.

      @param shell_type1
        Shell type of second contraction.

      @param r0
        Center of the first contraction. (size=3)

      @param r1
        Center of the second contraction. (size=3)

      @param point
        Cartesian coordinates of the grid point. (size=3)
   */
  void reset(long shell_type0, long shell_type1, const double* r0, const double* r1,
             const double* point);

  //! Convert results in work array from Cartesian to pure functions where needed.
  void cart_to_pure();

  //! Shell type of the first contraction.
  const long get_shell_type0() const {return shell_type0;}

  //! Shell type of the second contraction.
  const long get_shell_type1() const {return shell_type1;}

  /** @brief
        Add contributions to work array for current grid point and given pair of primitive
        shells.

      @param coeff
        Product of contraction coefficients of current primitive shells.

      @param alpha0
        Exponent of the first shell.

      @param alpha1
        Exponent of the second shell.

      @param scales0
        Normalization prefactors for primitives in first shell.

      @param scales1
        Normalization prefactors for primitives in second shell.
    */
  virtual void add(double coeff, double alpha0, double alpha1, const double* scales0,
                   const double* scales1) = 0;

 protected:
  long shell_type0;     //!< Shell type of the first contraction.
  long shell_type1;     //!< Shell type of the second contraction.
  const double *r0;     //!< Center of first basis contraction.
  const double *r1;     //!< Center of second basis contraction.
  const double *point;  //!< Grid point at which the fn is evaluated.
  IterPow2 i2p;         //!< Double loop iterator over Cartesian powers of given ang mom.
};


/** @brief
      Calculator for Hartree potential on a set of grid points.

    This requires a double loop over all basis functions because the contributions from
    the products of basis functions cannot completely be factorized.

    Content of work_basis (at one grid point):
      [0] The Hartree potential due to the product of the two basis functions.
 */
class GB2DMGridHartreeFn : public GB2DMGridFn  {
 public:
  //! Construct a GB2DMGridHartreeFn object. (See base class for details.)
  explicit GB2DMGridHartreeFn(long max_shell_type);
  ~GB2DMGridHartreeFn();

  /** @brief
        Add contributions to work array for current grid point and given pair of primitive
        shells. (See base class for details.)
    */
  virtual void add(double coeff, double alpha0, double alpha1, const double* scales0,
                   const double* scales1);

 private:
  double* work_g0;    //!< Work array for results from NAI helper code (x).
  double* work_g1;    //!< Work array for results from NAI helper code (y).
  double* work_g2;    //!< Work array for results from NAI helper code (z).
  double* work_boys;  //!< Work array for Boys function values for different m values.
};


#endif  // HORTON_GBASIS_FNS_H_
