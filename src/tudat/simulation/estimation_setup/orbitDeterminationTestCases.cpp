#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"

namespace tudat
{
namespace unit_tests
{

Eigen::VectorXd getDefaultInitialParameterPerturbation( )
{
    Eigen::VectorXd parameterPerturbations = Eigen::VectorXd( 7 );
    for( int i = 0; i < 3; i++ )
    {
        parameterPerturbations( i ) = 1.0E3;
        parameterPerturbations( i + 3 ) = 1.0E-2;
    }
    parameterPerturbations( 6 ) = 5.0E6;

    return parameterPerturbations;
}

// template std::pair< std::shared_ptr< EstimationOutput< double > >, Eigen::VectorXd > executePlanetaryParameterEstimation< double, double
// >(
//         const int observableType,
//         Eigen::VectorXd parameterPerturbation,
//         Eigen::MatrixXd inverseAPrioriCovariance,
//         const double weight );

// template Eigen::VectorXd executeEarthOrbiterParameterEstimation< double, double >(
//         std::pair< std::shared_ptr< EstimationOutput< double > > , std::shared_ptr< EstimationInput< double, double > > >& podData,
//         const double startTime,
//         const int numberOfDaysOfData,
//         const int numberOfIterations,
//         const bool useFullParameterSet );

// template std::pair< Eigen::VectorXd, bool > executeEarthOrbiterBiasEstimation< double, double >(
//         const bool estimateRangeBiases,
//         const bool estimateTwoWayBiases,
//         const bool useSingleBiasModel,
//         const bool estimateAbsoluteBiases,
//         const bool omitRangeData,
//         const bool useMultiArcBiases );

}  // namespace unit_tests

}  // namespace tudat
