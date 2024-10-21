#include "tudat/astro/orbit_determination/observation_partials/clockParameterPartials.h"


namespace tudat
{

namespace observation_partials
{


Eigen::Matrix< double, 1, Eigen::Dynamic > TimingPartialWrtGlobalPolynomialCorrectionCoefficients::getPartialOfClockErrorWrtParameter(
        const double time )
{
    double timeIntoCurrentArc = timingSystem_->getTimeIntoCurrentArc( time );
    Eigen::Matrix< double, 1, Eigen::Dynamic > partialVector = Eigen::Matrix< double, 1, Eigen::Dynamic >::Zero( 1, polynomialPowers_.size( ) );

    for( unsigned int i = 0; i < polynomialPowers_.size( ); i++ )
    {
        partialVector( 0, i ) = std::pow( timeIntoCurrentArc, polynomialPowers_[ i ] );
    }

    return partialVector;
}


Eigen::Matrix< double, 1, Eigen::Dynamic > TimingPartialWrtArcwisePolynomialCorrectionCoefficients::getPartialOfClockErrorWrtParameter(
        const double time )
{
    std::pair< Time, int > timeIntoCurrentArcAndArcIndex = timingSystem_->getTimeIntoCurrentArcAndArcIndex( time );

    Eigen::Matrix< double, 1, Eigen::Dynamic > partialVector = Eigen::Matrix< double, 1, Eigen::Dynamic >::Zero(
                1, numberOfArcs_ * polynomialPowers_.size( ) );

    arcFindIterator_ = std::find( arcIndices_.begin( ), arcIndices_.end( ), timeIntoCurrentArcAndArcIndex.second );
    if( arcFindIterator_ != arcIndices_.end( ) )
    {
        int indexInArcList = std::distance( arcIndices_.begin( ), arcFindIterator_ );
        int startIndex = indexInArcList * polynomialPowers_.size( );

        for( unsigned int i = startIndex; i < startIndex + polynomialPowers_.size( ); i++ )
        {
            partialVector( 0, i ) = std::pow( static_cast< double >( timeIntoCurrentArcAndArcIndex.first ), polynomialPowers_[ i - startIndex ] );
        }

    }
    return partialVector;
}


}

}
