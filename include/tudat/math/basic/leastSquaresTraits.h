/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *      Madsen, K., Nielsen, H., and Tingleff, O., Methods for Non-Linear Least Squares Problems, 2nd ed.,
 *          Technical University of Denmark, Faculty of Informatics and Mathematical Modelling, April 2004.
 */

#ifndef TUDAT_LEASTSQUARESTRAITS_H
#define TUDAT_LEASTSQUARESTRAITS_H

#include <Eigen/SparseCore>
#include <optional>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

namespace tudat
{

namespace linear_algebra
{

struct StorageClass {};
struct Dense: public StorageClass {};
struct Sparse: public StorageClass {};

struct SolverClass {};
struct SVD: public SolverClass {};
struct QR: public SolverClass {};

/*! Traits class generalizing over different possible matrix representations.
 * Mainly used to control whether operations are performed on Dense or Sparse
 * matrices, but can also be used to switch to single precision.
 */
template <typename Real, typename Storage>
struct MatrixTraits {};

/*! Convenience traits class to get `MatrixTraits` from an Eigen matrix type.
 */
template <typename E>
struct from_eigen {};

template <typename T>
struct from_eigen<Eigen::MatrixX<T>> {
   using traits = MatrixTraits<T, Dense>;
   using value_type = typename traits::value_type;
   using dense_vector_type = Eigen::VectorX<T>;
};

template <typename T>
struct from_eigen<Eigen::SparseMatrix<T>> {
    using traits = MatrixTraits<T, Sparse>;
    using value_type = typename traits::value_type;
    using dense_vector_type = Eigen::VectorX<T>;
};

/*! Traits class generalizing over different solvers implemented in Eigen.
 */
template <typename Matrix, typename Solver>
struct SolverTraits {};

template <typename Real>
struct MatrixTraits<Real, Dense>
{
    using value_type = Real;
    using matrix_type = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

    /*! Function to set an element in a dense matrix.
     *
     * Usually it is more efficient to use a generating expression, instead of
     * setting individual elements.
     *
     * \param A matrix
     * \param i row
     * \param j column
     */
    inline static void set_element(matrix_type &A, unsigned i, unsigned j, double value) {
        A(i, j) = value;
    }

    inline static void set_zero(matrix_type &A, unsigned i, unsigned j) {
        A(i, j) = value_type(0);
    }

    /*! Extends a square matrix with a second one in a symmetric way. The matrix `A` is
     * grown to fit `B` in the extended rows, and `B^T` in the extended columns.
     */
    matrix_type extend_symmetric( matrix_type &A, matrix_type const &B ) {
        unsigned a_n     = A->cols( );
        unsigned b_nrows = B->rows( );
        unsigned b_ncols = B->cols( );

        if ( A->cols( ) != A.rows( ) )
        {
            std::ostringstream msg;
            msg << "Internal error: expected square matrix, got size " << A->rows() << " x " << A->cols();
            throw std::runtime_error( msg.str( ) );
        }

        if ( b_ncols != a_n )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible with partials" );
        }

        A.conservativeResize( a_n + b_nrows, a_n + b_nrows );
        A.block( a_n,   0, b_nrows, b_ncols ) = B;
        A.block(   0, a_n, b_ncols, b_nrows ) = B.transpose( );
        A.block( a_n, a_n, b_nrows, b_nrows ).setZero( );
    }

    template <typename SparseRange>  // requires std::range<R>
    inline static void set_from_triplets(matrix_type &A, SparseRange const &rng) {
        for (auto &[i, j, v] : rng) {
            A(i, j) = v;
        }
    }

    /*! Function to generate a matrix from a function.
     *
     * This creates a matrix and only fills those elements for which the given
     * function returns a value.
     *
     * \param rows total number of rows.
     * \param cols total number of columns.
     * \param f The given function should take row and column arguments and
     * return a `std::optional<value_type>`.
     * \tparam Fn Type of the function argument.
     * \returns The generated matrix.
     */
    template <typename Fn>
    static matrix_type generate(unsigned rows, unsigned cols, Fn f) {
        matrix_type result{rows, cols};
        for (unsigned j = 0; j < rows; ++j) {
            for (unsigned i = 0; i < cols; ++i) {
                if (std::optional<value_type> v = f(j, i)) {
                    result(j, i) = *v;
                }
            }
        }
        return result;
    }

    /*! Fill a matrix using a range.
     *
     * Fills a matrix using a range that iterates 3-tuples of row, column, value.
     * This is the preferred way to fill a sparse matrix, since a range can hint at its
     * size.
     *
     * \param rows Number of rows in the matrix.
     * \param cols Number of columns in the matrix.
     * \param r Range over [row, col, value] triples.
     * \returns The generated matrix.
     */
    template <typename SparseRange>
    static matrix_type from_sparse_range(unsigned rows, unsigned cols, SparseRange const &r) {
        matrix_type result{rows, cols};  // default initializes elements to 0
        for (auto &[i, j, v] : r) {
            result(i, j) = v;
        }
        return result;
    }

    //! Function that compresses a sparse matrix. This is needed for sparse
    //! matrices, empty for dense matrices.
    inline static void make_compressed(matrix_type &A) {}
};

template <typename Real>
struct MatrixTraits<Real, Sparse>
{
    using value_type = Real;
    using matrix_type = Eigen::SparseMatrix<Real>;

    /*! Function to set an element in a sparse matrix.
     *
     * Usually it is more efficient to use a generating expression, instead of
     * setting individual elements.
     *
     * \param A matrix
     * \param i row
     * \param j column
     */
    inline static void set_element(matrix_type &A, unsigned i, unsigned j, double value) {
        A.insert(i, j) = value;
    }

    /*! Function to set an element to zero.
     *
     * In sparse matrices this is a noop.
     */
    inline static void set_zero(matrix_type &A, unsigned i, unsigned j) {}

    template <typename SparseRange>  // requires std::range<R>
    inline static void set_from_triplets(matrix_type &A, SparseRange const &rng) {
        A.setFromTriplets(rng.begin(), rng.end());
    }

    //! Function that compresses a sparse matrix. Empty for dense matrices.
    inline static void make_compressed(matrix_type &A) {
        A.makeCompressed();
    }

    /*! Fill a matrix using a range.
     *
     * Fills a matrix using a range that iterates 3-tuples of row, column, value.
     * This is the preferred way to fill a sparse matrix, since a range can hint at its
     * size.
     *
     * \param rows Number of rows in the matrix.
     * \param cols Number of columns in the matrix.
     * \param r Range over [row, col, value] triples.
     * \returns The generated matrix.
     */
    template <typename SparseRange>
    static matrix_type from_sparse_range(unsigned rows, unsigned cols, SparseRange const &r) {
        matrix_type result{rows, cols};
        result.setFromTriplets(r.begin(), r.end());
        return result;
    }
};

}  // namespace linear_algebra

}  // namespace tudat

#endif
