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

// UPDATELIBDOCTITLE: Evaluation of functions expanded in a Gaussian basis

#ifndef HORTON_GBASIS_FNS_H_
#define HORTON_GBASIS_FNS_H_

#include "horton/gbasis/calc.h"
#include "horton/gbasis/common.h"
#include "horton/gbasis/iter_pow.h"


/** @brief
      Base class for grid calculators that require only a single loop over all basis
      functions.

    All these functions split the work into two distinct steps. These steps are called
    from functions in gbasis.cpp (GOBasis). For a given grid point...

    1) The add method is called several times with all Gaussian primitives to evaluate
       properties of these primitives on the grid points. After this stage, the code in
       gbasis.cpp will take care of combining everything into proper contractions. It
       also calls the cart_to_pure method to get transform the results to pure functions
       when needed.

    2) (optional)
       The compute_point_from_* method is called to convert the properties of the basis
       functions into the function of interest, making use of either expansion
       coefficients of the orbitals or first-order density matrix coefficients.

    3) (optional)
       The compute_fock_from_pot method is called to convert a potential on a grid into
       a Fock operator.
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
        The number of results for each grid point, e.g. 3 for a density gradient.
    */
  GB1ExpGridFn(long max_shell_type, long nfn, long dim_work, long dim_output)
      : GB1GridFn(max_shell_type, dim_work, dim_output), nfn(nfn) {}

  /** @brief
        Compute (final) results for a given grid point.

      @param work_basis
        Properties of basis functions computed for the current grid point. (Work done by
        add method.)

      @param coeffs
        The orbital expansion coefficients.

      @param nbasis
        The number of basis functions.

      @param output
        The output array for the current grid point.
    */
  virtual void compute_point_from_exp(double* work_basis, double* coeffs, long nbasis,
                                      double* output) = 0;

 protected:
  long nfn;  //!< Number of orbitals (occupied and virtual).
};


/** @brief
      Evaluates a selection of orbitals on a grid.
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
        add method.)

      @param dm
        The coefficients of the first-order density matrix.

      @param nbasis
        The number of basis functions.

      @param output
        The output array for the current grid point.

      @param epsilon
        A cutoff value used to discard small contributions.

      @param dmmaxrow
        The maximum value of the density matrix on each row.
    */
  virtual void compute_point_from_dm(double* work_basis, double* dm, long nbasis,
                                     double* output, double epsilon, double* dmmaxrow) = 0;

  /** @brief
        Add contribution to Fock matrix for one grid point.

      @param pot
        The value of the potential on the grid point. (This may be multiple values, e.g in
        case of GGA.)

      @param work_basis
        Properties of the orbital basis in the current grid point.

      @param nbasis
        The number of basis functions.

      @param output
        The Fock matrix to which the result will be added.
    */
  virtual void compute_fock_from_pot(double* pot, double* work_basis, long nbasis,
                                     double* output) = 0;
};


//! Compute just the electron density on a grid.
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
                                     double* output);

 private:
  double poly_work[MAX_NCART_CUMUL];  //!< Work array with Cartesian polynomials.
  long offset;  //!< Offset for the polynomials for the density
};


//! Compute gradient of the electron density on a grid.
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
                                     double* output);

 protected:
  double poly_work[MAX_NCART_CUMUL_D];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
};


//! Compute density and gradient on a grid.
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
                                     double* output);
};


//! Compute kinetic energy density on a grid.
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
                                     double* output);

 private:
  double poly_work[MAX_NCART_CUMUL_D];  //!< Work array with Cartesian polynomials.
  long offset;     //!< Offset for the polynomials for the density.
  long offset_l1;  //!< Lower offset for the polynomials for the gradient.
  long offset_h1;  //!< Higher offset for the polynomials for the gradient.
};


/** @brief
      Compute density Hessian on a grid: xx, xy, xz, yy, yz, zz.
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
                                     double* output);

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
                                     double* output);

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
