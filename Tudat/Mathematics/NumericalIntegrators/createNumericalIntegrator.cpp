#include "Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h"
namespace tudat
{

namespace numerical_integrators
{

template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::VectorXd,
Eigen::VectorXd, double > > createIntegrator< double, Eigen::VectorXd, double >(
        std::function< Eigen::VectorXd( const double, const Eigen::VectorXd& ) > stateDerivativeFunction,
        const Eigen::VectorXd initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 >,
Eigen::Matrix< long double, Eigen::Dynamic, 1 >, double > > createIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, double >(
        std::function< Eigen::Matrix< long double, Eigen::Dynamic, 1 >(
            const double, const Eigen::Matrix< long double, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
        const Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::VectorXd,
Eigen::VectorXd, long double > > createIntegrator< Time, Eigen::VectorXd, long double >(
        std::function< Eigen::VectorXd( const Time, const Eigen::VectorXd& ) > stateDerivativeFunction,
        const Eigen::VectorXd initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );

template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >,
Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double > > createIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, 1 >, long double >(
        std::function< Eigen::Matrix< long double, Eigen::Dynamic, 1 >(
            const Time, const Eigen::Matrix< long double, Eigen::Dynamic, 1 >& ) > stateDerivativeFunction,
        const Eigen::Matrix< long double, Eigen::Dynamic, 1 > initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );




template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::MatrixXd,
Eigen::MatrixXd, double > > createIntegrator< double, Eigen::MatrixXd, double >(
        std::function< Eigen::MatrixXd( const double, const Eigen::MatrixXd& ) > stateDerivativeFunction,
        const Eigen::MatrixXd initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

template std::shared_ptr< numerical_integrators::NumericalIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >,
Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, double > > createIntegrator< double, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, double >(
        std::function< Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >(
            const double, const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction,
        const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > initialState, std::shared_ptr< IntegratorSettings< double > > integratorSettings );

template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::MatrixXd,
Eigen::MatrixXd, long double > > createIntegrator< Time, Eigen::MatrixXd, long double >(
        std::function< Eigen::MatrixXd( const Time, const Eigen::MatrixXd& ) > stateDerivativeFunction,
        const Eigen::MatrixXd initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );

template std::shared_ptr< numerical_integrators::NumericalIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >,
Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, long double > > createIntegrator< Time, Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >, long double >(
        std::function< Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >(
            const Time, const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic >& ) > stateDerivativeFunction,
        const Eigen::Matrix< long double, Eigen::Dynamic, Eigen::Dynamic > initialState, std::shared_ptr< IntegratorSettings< Time > > integratorSettings );


} // namespace numerical_integrators

} // namespace tudat

