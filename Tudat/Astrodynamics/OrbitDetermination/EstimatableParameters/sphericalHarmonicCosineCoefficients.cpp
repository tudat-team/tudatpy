
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/sphericalHarmonicCosineCoefficients.h"

namespace tudat
{

namespace estimatable_parameters
{

std::map< int, std::pair< int, int > > convertShDegreeAndOrderLimitsToBlockIndices(
        const int minimumDegree, const int minimumOrder, const int maximumDegree, const int maximumOrder,
        int& totalSetSize )
{
    if( maximumOrder > maximumDegree )
    {
        std::cerr<<"Warning when converting sh coefficients to block indices;l maximum order is larger than maximum degree"<<std::endl;
    }

    std::map< int, std::pair< int, int > > blockIndices;

    std::pair< int, int > currentPair;
    int currentOrderLength;
    totalSetSize = 0;
    for( int i = minimumDegree; i <= maximumDegree; i++ )
    {
        currentPair.first = minimumOrder;
        currentOrderLength = 0;
        for( int j = minimumOrder; ( j <= i && j <= maximumOrder ); j++ )
        {
            totalSetSize++;
            currentOrderLength++;
        }
        currentPair.second = currentOrderLength;
        blockIndices[ i ] = currentPair;
    }
    return blockIndices;
}

std::map< int, std::pair< int, int > > convertShDegreeAndOrderLimitsToBlockIndices(
        const int minimumDegree, const int minimumOrder, const int maximumDegree, const int maximumOrder )
{
    int totalSetSize;
    return convertShDegreeAndOrderLimitsToBlockIndices( minimumDegree, minimumOrder, maximumDegree, maximumOrder, totalSetSize );
}


SphericalHarmonicsCosineCoefficients::SphericalHarmonicsCosineCoefficients(
        boost::function< Eigen::MatrixXd( ) > getCosineCoefficients,
        boost::function< void( Eigen::MatrixXd ) > setCosineCoefficients,
        const int minimumDegree, const int minimumOrder,
        const int maximumDegree, const int maximumOrder,
        std::string& associatedBody ):
    EstimatableParameter< Eigen::VectorXd >( spherical_harmonics_cosine_coefficient_block, associatedBody ),
    getCosineCoefficients_( getCosineCoefficients ), setCosineCoefficients_( setCosineCoefficients ),
    minimumDegree_( minimumDegree ), minimumOrder_( minimumOrder ),
    maximumDegree_( maximumDegree ), maximumOrder_( maximumOrder )
{
    if( maximumOrder_ > maximumDegree_ )
    {
        maximumOrder_ = maximumDegree_;
    }

    blockIndices_ = convertShDegreeAndOrderLimitsToBlockIndices(
                minimumDegree_, minimumOrder_, maximumDegree_, maximumOrder_, parameterSize_ );

}


Eigen::VectorXd SphericalHarmonicsCosineCoefficients::getParameterValue( )
{
    Eigen::VectorXd parameterVector = Eigen::VectorXd::Zero( parameterSize_ );
    int currentIndex = 0;

    Eigen::MatrixXd coefficientBlock =
            getCosineCoefficients_( );

    for(  std::map< int, std::pair< int, int > >::iterator blockIterator = blockIndices_.begin( );
          blockIterator != blockIndices_.end( ); blockIterator++ )
    {
        parameterVector.segment( currentIndex, blockIterator->second.second ) = coefficientBlock.block(
                    blockIterator->first, blockIterator->second.first, 1, blockIterator->second.second ).transpose( );
        currentIndex += blockIterator->second.second;
    }
    return parameterVector;

}

void SphericalHarmonicsCosineCoefficients::setParameterValue( const Eigen::VectorXd parameterValue )
{
    Eigen::MatrixXd coefficients = getCosineCoefficients_( );
    int currentIndex = 0;

    for(  std::map< int, std::pair< int, int > >::iterator blockIterator = blockIndices_.begin( );
          blockIterator != blockIndices_.end( ); blockIterator++ )
    {
        coefficients.block( blockIterator->first, blockIterator->second.first, 1, blockIterator->second.second ) =
                parameterValue.segment( currentIndex, blockIterator->second.second ).transpose( );
        currentIndex += blockIterator->second.second;
    }
    setCosineCoefficients_( coefficients );
}

}

}


