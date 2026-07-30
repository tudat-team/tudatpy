#ifndef TUDATPY_QUAD_PRECISION_TYPE_CASTERS_H
#define TUDATPY_QUAD_PRECISION_TYPE_CASTERS_H

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <tudat/config.hpp>

#if TUDAT_HIGH_PRECISION_STATE_SCALAR_IS_CPP_BIN_FLOAT_QUAD

namespace pybind11
{
namespace detail
{

//! Convert the configured quad scalar to and from Python's binary64 float.
template<>
struct type_caster< tudat::HighPrecisionStateScalar > {
    using QuadScalar = tudat::HighPrecisionStateScalar;

    PYBIND11_TYPE_CASTER( QuadScalar, const_name( "float" ) );

    bool load( handle source, bool convert )
    {
        if( !source || ( !convert && !PyFloat_Check( source.ptr( ) ) ) )
        {
            return false;
        }

        const double doubleValue = PyFloat_AsDouble( source.ptr( ) );
        if( PyErr_Occurred( ) )
        {
            PyErr_Clear( );
            return false;
        }

        value = QuadScalar( doubleValue );
        return true;
    }

    static handle cast( const QuadScalar& source, return_value_policy, handle )
    {
        return PyFloat_FromDouble( static_cast< double >( source ) );
    }
};

//! Convert dense Eigen matrices containing quad scalars through float64 NumPy arrays.
template< int Rows, int Columns, int Options, int MaximumRows, int MaximumColumns >
struct type_caster<
        Eigen::Matrix< tudat::HighPrecisionStateScalar, Rows, Columns, Options, MaximumRows, MaximumColumns >,
        enable_if_t< is_eigen_dense_plain<
                Eigen::Matrix< tudat::HighPrecisionStateScalar, Rows, Columns, Options, MaximumRows, MaximumColumns > >::value > > {
    using QuadScalar = tudat::HighPrecisionStateScalar;
    using QuadMatrix = Eigen::Matrix< QuadScalar, Rows, Columns, Options, MaximumRows, MaximumColumns >;
    using DoubleMatrix = Eigen::Matrix< double, Rows, Columns, Options, MaximumRows, MaximumColumns >;
    using DoubleCaster = make_caster< DoubleMatrix >;

    PYBIND11_TYPE_CASTER( QuadMatrix, DoubleCaster::name );

    bool load( handle source, bool convert )
    {
        DoubleCaster doubleCaster;
        if( !doubleCaster.load( source, convert ) )
        {
            return false;
        }

        const DoubleMatrix& doubleValue = doubleCaster;
        value = doubleValue.template cast< QuadScalar >( );
        return true;
    }

    static handle cast( const QuadMatrix& source, return_value_policy policy, handle parent )
    {
        DoubleMatrix doubleValue = source.template cast< double >( );
        return DoubleCaster::cast( std::move( doubleValue ), policy, parent );
    }

    static handle cast( QuadMatrix&& source, return_value_policy policy, handle parent )
    {
        DoubleMatrix doubleValue = source.template cast< double >( );
        return DoubleCaster::cast( std::move( doubleValue ), policy, parent );
    }

    static handle cast( const QuadMatrix* source, return_value_policy policy, handle parent )
    {
        return source == nullptr ? none( ).release( ) : cast( *source, policy, parent );
    }

    static handle cast( QuadMatrix* source, return_value_policy policy, handle parent )
    {
        return source == nullptr ? none( ).release( ) : cast( *source, policy, parent );
    }
};

}  // namespace detail
}  // namespace pybind11

#endif

#endif  // TUDATPY_QUAD_PRECISION_TYPE_CASTERS_H
