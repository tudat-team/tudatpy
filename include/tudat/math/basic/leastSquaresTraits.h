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

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>

namespace tudat
{

namespace linear_algebra
{

struct Dense {};
struct Sparse {};

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

template <typename Real>
struct MatrixTraits<Real, Dense>
{
    using storage_type = Dense;
    using value_type = Real;
    using matrix_type = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;

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
        matrix_type result = matrix_type::Zero(rows, cols);
        for (const auto& entry : r) {
            result(entry.row(), entry.col()) = entry.value();
        }
        return result;
    }

    //! Function that compresses a sparse matrix. This is needed for sparse
    //! matrices, empty for dense matrices.
    inline static void make_compressed(matrix_type &A) {}

    static Eigen::VectorX<Real> normalize_columns(matrix_type& A) {
        Eigen::VectorX<Real> normalizationTerms = Eigen::VectorX<Real>::Ones(A.cols());
        for (int column = 0; column < A.cols(); column++) {
            Eigen::VectorX<Real> currentVector = A.block(0, column, A.rows(), 1);
            const Real minimum = currentVector.minCoeff();
            const Real maximum = currentVector.maxCoeff();
            if (std::fabs(minimum) > maximum) {
                normalizationTerms(column) = minimum;
            }
            else if (maximum != Real(0)) {
                normalizationTerms(column) = maximum;
            }

            A.block(0, column, A.rows(), 1) = currentVector / normalizationTerms(column);
        }
        return normalizationTerms;
    }

    static void multiply_columns(matrix_type& A, const Eigen::VectorX<Real>& factors) {
        for (int column = 0; column < A.cols(); column++) {
            A.block(0, column, A.rows(), 1) *= factors(column);
        }
    }

    static void multiply_rows_by_sqrt(matrix_type& A, const Eigen::VectorX<Real>& factors) {
        for (int row = 0; row < A.rows(); row++) {
            A.block(row, 0, 1, A.cols()) *= std::sqrt(factors(row));
        }
    }
};

template <typename Real>
struct MatrixTraits<Real, Sparse>
{
    using storage_type = Sparse;
    using value_type = Real;
    using matrix_type = Eigen::SparseMatrix<Real>;

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

    static Eigen::VectorX<Real> normalize_columns(matrix_type& A) {
        Eigen::VectorX<Real> normalizationTerms = Eigen::VectorX<Real>::Ones(A.cols());
        for (int column = 0; column < A.outerSize(); column++) {
            Real minimum = Real(0);
            Real maximum = Real(0);
            for (typename matrix_type::InnerIterator iterator(A, column); iterator; ++iterator) {
                minimum = std::min(minimum, iterator.value());
                maximum = std::max(maximum, iterator.value());
            }

            if (std::fabs(minimum) > maximum) {
                normalizationTerms(column) = minimum;
            }
            else if (maximum != Real(0)) {
                normalizationTerms(column) = maximum;
            }

            for (typename matrix_type::InnerIterator iterator(A, column); iterator; ++iterator) {
                iterator.valueRef() /= normalizationTerms(column);
            }
        }
        A.makeCompressed();
        return normalizationTerms;
    }

    static void multiply_columns(matrix_type& A, const Eigen::VectorX<Real>& factors) {
        for (int column = 0; column < A.outerSize(); column++) {
            for (typename matrix_type::InnerIterator iterator(A, column); iterator; ++iterator) {
                iterator.valueRef() *= factors(column);
            }
        }
        A.makeCompressed();
    }

    static void multiply_rows_by_sqrt(matrix_type& A, const Eigen::VectorX<Real>& factors) {
        for (int column = 0; column < A.outerSize(); column++) {
            for (typename matrix_type::InnerIterator iterator(A, column); iterator; ++iterator) {
                iterator.valueRef() *= std::sqrt(factors(iterator.row()));
            }
        }
        A.makeCompressed();
    }
};

}  // namespace linear_algebra

}  // namespace tudat

#endif
