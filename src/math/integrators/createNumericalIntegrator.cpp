#include "tudat/math/integrators/createNumericalIntegrator.h"
namespace tudat
{

namespace numerical_integrators
{

std::vector< std::tuple< int, int, int, int > > getStandardCartesianStatesElementsToCheck(
    const int numberOfRows, const int numberOfColumns )
{
    int stateColumn = numberOfColumns - 1;

    if( numberOfRows % 6 != 0 )
    {
        throw std::runtime_error( "Error when getting standard Cartesian element blocks for step-size control; propagated state has incompatible number of rows: " +
                                  std::to_string( numberOfRows ) );
    }

    std::vector< std::tuple< int, int, int, int > > blocks;
    for( int i = 0; i < numberOfRows / 3; i++ )
    {
        blocks.push_back( { i * 3, stateColumn, 3, 1 } );
    }
    return blocks;    if( numberOfColumns != 1 )
    {
        throw std::runtime_error( "Error when getting standard rotational state element blocks for step-size control; propagated state has more than 1 column." );
    }
}

std::vector< std::tuple< int, int, int, int > > getStandardRotationalStatesElementsToCheck(
    const int numberOfRows, const int numberOfColumns )
{
    int stateColumn = numberOfColumns - 1;


    if( numberOfRows % 7 != 0 )
    {
        throw std::runtime_error( "Error when getting standard rotational state element blocks for step-size control; propagated state has incompatible number of rows: " +
                                  std::to_string( numberOfRows ) );
    }

    std::vector< std::tuple< int, int, int, int > > blocks;
    for( int i = 0; i < numberOfRows / 7; i++ )
    {
        blocks.push_back( { i * 7, stateColumn, 4, 1 } );
        blocks.push_back( { i * 7 + 4, stateColumn, 3, 1 } );
    }
    return blocks;
}

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::VectorXd,
//Eigen::VectorXd, double > > createIntegrator< double, Eigen::VectorXd, double >(
//        std::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) > stateDerivativeFunction,
//        const Eigen::VectorXd initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 >,
//Eigen::Matrix< long double, Eigen::Dynamic, 1 >, double > > createIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, double >(
//        std::function< Eigen::Matrix< long double, Eigen::Dynamic, 1 >(
//            const double, const Eigen::Matrix< long double, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
//        const Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::VectorXd,
//Eigen::VectorXd, long double > > createIntegrator< Time, Eigen::VectorXd, long double >(
//        std::function< Eigen::VectorXd( const Time, const Eigen::VectorXd& ) > stateDerivativeFunction,
//        const Eigen::VectorXd initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >,
//Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double > > createIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double >(
//        std::function< Eigen::Matrix< long double, Eigen::Dynamic, 1 >(
//            const Time, const Eigen::Matrix< long double, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
//        const Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );




//template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd,
//Eigen::MatrixXd, double > > createIntegrator< double, Eigen::MatrixXd, double >(
//        std::function< Eigen::MatrixXd( const double, const Eigen::MatrixXd& ) > stateDerivativeFunction,
//        const Eigen::MatrixXd initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >,
//Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, double > > createIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, double >(
//        std::function< Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >(
//            const double, const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction,
//        const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::MatrixXd,
//Eigen::MatrixXd, long double > > createIntegrator< Time, Eigen::MatrixXd, long double >(
//        std::function< Eigen::MatrixXd( const Time, const Eigen::MatrixXd& ) > stateDerivativeFunction,
//        const Eigen::MatrixXd initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );

//template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >,
//Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, long double > > createIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, long double >(
//        std::function< Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >(
//            const Time, const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction,
//        const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );


} // namespace numerical_integrators

} // namespace tudat

