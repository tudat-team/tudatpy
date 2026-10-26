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

#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SVD>

namespace tudat
{

namespace linear_algebra
{

struct StorageClass {
};
struct Dense : public StorageClass {
};
struct Sparse : public StorageClass {
};

struct SolverClass {
};
struct SVD : public SolverClass {
};
struct QR : public SolverClass {
};

/*! Traits class generalizing over different possible matrix representations.
 * Mainly used to control whether operations are performed on Dense or Sparse
 * matrices, but can also be used to switch to single precision.
 */
template< typename Real, typename Storage >
struct MatrixTraits {
};

/*! Traits class generalizing over different solvers implemented in Eigen.
 */
template< typename Matrix, typename Solver >
struct SolverTraits {
};

/*! Convenience traits class to get `MatrixTraits` and `SolverTraits` defaults from an Eigen matrix type.
 */
template< typename E >
struct from_eigen {
};

template< typename T >
struct from_eigen< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > > {
    using traits = MatrixTraits< T, Dense >;
    using value_type = T;
    using dense_vector_type = Eigen::Matrix< T, Eigen::Dynamic, 1 >;
    using solver = SVD;
};

template< typename T >
struct from_eigen< Eigen::SparseMatrix< T > > {
    using traits = MatrixTraits< T, Sparse >;
    using value_type = T;
    using dense_vector_type = Eigen::Matrix< T, Eigen::Dynamic, 1 >;
    using solver = QR;
};

template< typename Real >
struct MatrixTraits< Real, Dense > {
    using value_type = Real;
    using matrix_type = Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic >;

    /*! Function to set an element in a dense matrix.
     *
     * Usually it is more efficient to use a generating expression, instead of
     * setting individual elements.
     *
     * \param A matrix
     * \param i row
     * \param j column
     */
    inline static void set_element( matrix_type& A, unsigned i, unsigned j, Real value )
    {
        A( i, j ) = value;
    }

    inline static void set_zero( matrix_type& A, unsigned i, unsigned j )
    {
        A( i, j ) = value_type( 0 );
    }

    //! Function to create a matrix of zeros of the given size.
    inline static matrix_type zeros( unsigned rows, unsigned cols )
    {
        return matrix_type::Zero( rows, cols );
    }

    //! Function to convert a dense matrix to this storage. Passes the matrix through without copying.
    inline static const matrix_type& from_dense( const matrix_type& A )
    {
        return A;
    }

    //! Function to convert a sparse matrix to this storage.
    inline static matrix_type from_sparse( const Eigen::SparseMatrix< Real >& A )
    {
        return matrix_type( A );
    }

    //! Function to convert a matrix in this storage to a dense matrix.
    inline static matrix_type to_dense( matrix_type A )
    {
        return A;
    }

    /*! Extends a square matrix with a second one in a symmetric way. The matrix `A` is
     * grown to fit `B` in the extended rows, and `B^T` in the extended columns.
     */
    inline static void extend_symmetric( matrix_type& A, matrix_type const& B )
    {
        const Eigen::Index n = A.cols( );
        const Eigen::Index c = B.rows( );

        if( A.rows( ) != n )
        {
            std::ostringstream msg;
            msg << "Internal error: expected square matrix, got size " << A.rows( ) << " x " << A.cols( );
            throw std::runtime_error( msg.str( ) );
        }

        if( B.cols( ) != n )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible with partials" );
        }

        A.conservativeResize( n + c, n + c );
        A.block( n, 0, c, n ) = B;
        A.block( 0, n, n, c ) = B.transpose( );
        A.block( n, n, c, c ).setZero( );
    }

    //! Adds `B` to the top-left block of `A` (of the same size as `B`).
    inline static void add_to_top_left( matrix_type& A, matrix_type const& B )
    {
        A.topLeftCorner( B.rows( ), B.cols( ) ) += B;
    }

    template< typename SparseRange >  // requires std::range<R>
    inline static void set_from_triplets( matrix_type& A, SparseRange const& rng )
    {
        for( auto const& [ i, j, v ] : rng )
        {
            A( i, j ) = v;
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
    template< typename Fn >
    static matrix_type generate( unsigned rows, unsigned cols, Fn f )
    {
        matrix_type result = matrix_type::Zero( rows, cols );
        for( unsigned i = 0; i < rows; ++i )
        {
            for( unsigned j = 0; j < cols; ++j )
            {
                if( std::optional< value_type > v = f( i, j ) )
                {
                    result( i, j ) = *v;
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
    template< typename SparseRange >
    static matrix_type from_sparse_range( unsigned rows, unsigned cols, SparseRange const& r )
    {
        matrix_type result = matrix_type::Zero( rows, cols );
        for( auto const& [ i, j, v ] : r )
        {
            result( i, j ) = v;
        }
        return result;
    }

    //! Function that compresses a sparse matrix. This is needed for sparse
    //! matrices, empty for dense matrices.
    inline static void make_compressed( matrix_type& A ) {}
};

template< typename Real >
struct MatrixTraits< Real, Sparse > {
    using value_type = Real;
    using matrix_type = Eigen::SparseMatrix< Real >;

    /*! Function to set an element in a sparse matrix.
     *
     * Usually it is more efficient to use a generating expression, instead of
     * setting individual elements.
     *
     * \param A matrix
     * \param i row
     * \param j column
     */
    inline static void set_element( matrix_type& A, unsigned i, unsigned j, Real value )
    {
        A.insert( i, j ) = value;
    }

    /*! Function to set an element to zero.
     *
     * In sparse matrices this is a noop.
     */
    inline static void set_zero( matrix_type& A, unsigned i, unsigned j ) {}

    //! Function to create an empty (all-zero) sparse matrix of the given size.
    inline static matrix_type zeros( unsigned rows, unsigned cols )
    {
        return matrix_type( rows, cols );
    }

    //! Function to convert a dense matrix to this storage. Only exact zeros are dropped.
    inline static matrix_type from_dense( const Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic >& A )
    {
        return A.sparseView( );
    }

    //! Function to convert a sparse matrix to this storage. Passes the matrix through without copying.
    inline static const matrix_type& from_sparse( const matrix_type& A )
    {
        return A;
    }

    //! Function to convert a matrix in this storage to a dense matrix.
    inline static Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic > to_dense( const matrix_type& A )
    {
        return Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic >( A );
    }

    /*! Extends a square sparse matrix with a second one in a symmetric way. Same semantics as the dense
     * version: `A` is grown to fit `B` in the extended rows, and `B^T` in the extended columns. Rebuilt
     * from triplets, since sparse matrices cannot be resized in place without losing their entries.
     */
    inline static void extend_symmetric( matrix_type& A, matrix_type const& B )
    {
        const Eigen::Index n = A.cols( );
        const Eigen::Index c = B.rows( );

        if( A.rows( ) != n )
        {
            std::ostringstream msg;
            msg << "Internal error: expected square matrix, got size " << A.rows( ) << " x " << A.cols( );
            throw std::runtime_error( msg.str( ) );
        }

        if( B.cols( ) != n )
        {
            throw std::runtime_error( "Error when performing constrained least-squares, constraints are incompatible with partials" );
        }

        std::vector< Eigen::Triplet< Real > > triplets;
        triplets.reserve( A.nonZeros( ) + 2 * B.nonZeros( ) );
        for( int k = 0; k < A.outerSize( ); ++k )
        {
            for( typename matrix_type::InnerIterator it( A, k ); it; ++it )
            {
                triplets.emplace_back( it.row( ), it.col( ), it.value( ) );
            }
        }
        for( int k = 0; k < B.outerSize( ); ++k )
        {
            for( typename matrix_type::InnerIterator it( B, k ); it; ++it )
            {
                triplets.emplace_back( n + it.row( ), it.col( ), it.value( ) );
                triplets.emplace_back( it.col( ), n + it.row( ), it.value( ) );
            }
        }

        A.resize( n + c, n + c );
        A.setFromTriplets( triplets.begin( ), triplets.end( ) );
    }

    //! Adds `B` to the top-left block of `A`; `B` is padded with zeros to the size of `A` if needed.
    inline static void add_to_top_left( matrix_type& A, matrix_type const& B )
    {
        if( A.rows( ) == B.rows( ) && A.cols( ) == B.cols( ) )
        {
            A += B;
        }
        else
        {
            matrix_type paddedB = B;
            paddedB.conservativeResize( A.rows( ), A.cols( ) );
            A += paddedB;
        }
    }

    template< typename SparseRange >  // requires std::range<R>
    inline static void set_from_triplets( matrix_type& A, SparseRange const& rng )
    {
        A.setFromTriplets( rng.begin( ), rng.end( ) );
    }

    //! Function that compresses a sparse matrix. Empty for dense matrices.
    inline static void make_compressed( matrix_type& A )
    {
        A.makeCompressed( );
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
    template< typename SparseRange >
    static matrix_type from_sparse_range( unsigned rows, unsigned cols, SparseRange const& r )
    {
        matrix_type result( rows, cols );
        result.setFromTriplets( r.begin( ), r.end( ) );
        return result;
    }
};

/*! Solver traits for dense matrices, using a singular value decomposition.
 *
 * Emits a warning on stdout when the condition number exceeds `limitConditionNumberForWarning`.
 * A NaN limit disables the check.
 */
template< typename Real >
struct SolverTraits< Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic >, SVD > {
    using matrix_type = Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic >;
    using vector_type = Eigen::Matrix< Real, Eigen::Dynamic, 1 >;

    static vector_type solve( const matrix_type& A, const vector_type& b, Real limitConditionNumberForWarning )
    {
        auto svd = A.template jacobiSvd< Eigen::ComputeThinU | Eigen::ComputeThinV >( );
        if( limitConditionNumberForWarning == limitConditionNumberForWarning )
        {
            const vector_type& singularValues = svd.singularValues( );
            Real conditionNumber = singularValues( 0 ) / singularValues( singularValues.rows( ) - 1 );
            if( conditionNumber > limitConditionNumberForWarning )
            {
                std::cerr << "Warning when performing least squares, condition number is " << conditionNumber << std::endl;
            }
        }
        return svd.solve( b );
    }
};

/*! Solver traits for sparse matrices, using a sparse QR decomposition.
 *
 * Eigen's sparse QR does not provide a condition number estimate, so the limit is ignored.
 */
template< typename Real >
struct SolverTraits< Eigen::SparseMatrix< Real >, QR > {
    using matrix_type = Eigen::SparseMatrix< Real >;
    using vector_type = Eigen::Matrix< Real, Eigen::Dynamic, 1 >;

    static vector_type solve( const matrix_type& A, const vector_type& b, Real /* limitConditionNumberForWarning */ )
    {
        Eigen::SparseQR< matrix_type, Eigen::COLAMDOrdering< int > > qr;
        qr.compute( A );
        if( qr.info( ) != Eigen::Success )
        {
            throw std::runtime_error( "Error when solving sparse least-squares system: sparse QR decomposition failed" );
        }
        return qr.solve( b );
    }
};

}  // namespace linear_algebra

}  // namespace tudat

#endif
