#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicSineCoefficients.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"

namespace tudat
{

namespace estimatable_parameters
{

SphericalHarmonicsSineCoefficients::SphericalHarmonicsSineCoefficients(
        const boost::function< Eigen::MatrixXd( ) > getSineCoefficients,
        const boost::function< void( Eigen::MatrixXd ) > setSineCoefficients,
        const int minimumDegree, const int minimumOrder,
        const int maximumDegree, const int maximumOrder,
        const std::string& associatedBody ):
    EstimatableParameter< Eigen::VectorXd >( spherical_harmonics_sine_coefficient_block, associatedBody ),
    getSineCoefficients_( getSineCoefficients ), setSineCoefficients_( setSineCoefficients ),
    minimumDegree_( minimumDegree ), minimumOrder_( minimumOrder ),
    maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
{
    if( maximumOrder_ > maximumDegree_ )
    {
        maximumOrder_ = maximumDegree_;
    }
    if( minimumOrder_ == 0 )
    {
        minimumOrder_ = 1;
    }

    blockIndices_ = convertShDegreeAndOrderLimitsToBlockIndices(
                minimumDegree_, minimumOrder_, maximumDegree_, maximumOrder_, parameterSize_ );

}


Eigen::VectorXd SphericalHarmonicsSineCoefficients::getParameterValue( )
{
    Eigen::VectorXd parameterVector = Eigen::VectorXd::Zero( parameterSize_ );
    int currentIndex = 0;

    Eigen::MatrixXd coefficientBlock = getSineCoefficients_( );

    for(  std::map< int, std::pair< int, int > >::iterator blockIterator = blockIndices_.begin( );
          blockIterator != blockIndices_.end( ); blockIterator++ )
    {
        parameterVector.segment( currentIndex, blockIterator->second.second ) = coefficientBlock.block(
                    blockIterator->first, blockIterator->second.first, 1, blockIterator->second.second ).transpose( );
        currentIndex += blockIterator->second.second;
    }
    return parameterVector;

}

void SphericalHarmonicsSineCoefficients::setParameterValue( Eigen::VectorXd parameterValue )
{
    int currentIndex = 0;

    Eigen::MatrixXd coefficients = getSineCoefficients_( );

    for(  std::map< int, std::pair< int, int > >::iterator blockIterator = blockIndices_.begin( );
          blockIterator != blockIndices_.end( ); blockIterator++ )
    {
        coefficients.block( blockIterator->first, blockIterator->second.first, 1, blockIterator->second.second ) =
                parameterValue.segment( currentIndex, blockIterator->second.second ).transpose( );
        currentIndex += blockIterator->second.second;
    }
    setSineCoefficients_( coefficients );
}


}

}
