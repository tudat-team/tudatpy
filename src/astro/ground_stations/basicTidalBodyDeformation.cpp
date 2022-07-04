#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Calculate displacement due to degree 2 tide, band-independet love and shida numbers.
Eigen::Vector3d calculateDegreeTwoBasicTidalDisplacement( const double gravitationalParameterRatio,
                                                          const Eigen::Vector3d& stationPositionUnitVector,
                                                          const Eigen::Vector3d& relativeBodyState,
                                                          const double bodyEquatorialRadius,
                                                          const double degreeTwoLoveNumber,
                                                          const double degreeTwoShidaNumber )
{
    // Calculate distance to body causing deformation.
    double bodyDistance = relativeBodyState.norm( );

    //Calculate unit vector in direction of body causing deformation.
    Eigen::Vector3d relativeBodyStateUnitVector = relativeBodyState.normalized( );

    // Pre-calculate inner product of unit vectors for efficiency.
    double positionsInnerProduct = relativeBodyStateUnitVector.dot( stationPositionUnitVector );

    // Calculate site displacement according to Eq. 7.5 of IERS Conventions 2010
    Eigen::Vector3d displacement = degreeTwoLoveNumber * stationPositionUnitVector *
            ( ( 3.0 * positionsInnerProduct * positionsInnerProduct - 1.0 ) / 2.0 ) +
            3.0 * degreeTwoShidaNumber * positionsInnerProduct *
            ( relativeBodyStateUnitVector - positionsInnerProduct * stationPositionUnitVector );

    // Scale result according to aforementioned Eq. 7.5
    double distanceRatio = bodyEquatorialRadius / bodyDistance;
    double distanceRatioSquared = distanceRatio * distanceRatio;

    return displacement * distanceRatioSquared * distanceRatioSquared * bodyDistance * gravitationalParameterRatio;
}


//! Calculate displacement due to degree 3 tide, band-independet love and shida numbers.
Eigen::Vector3d calculateDegreeThreeBasicTidalDisplacement( const double gravitationalParameterRatio,
                                                            const Eigen::Vector3d& stationPositionUnitVector,
                                                            const Eigen::Vector3d& relativeBodyState,
                                                            const double bodyEquatorialRadius,
                                                            const double degreeThreeLoveNumber,
                                                            const double degreeThreeShidaNumber )
{
    // Calculate distance to body causing deformation.
    double bodyDistance = relativeBodyState.norm( );

    //Calculate unit vector in direction of body causing deformation.
    Eigen::Vector3d relativeBodyStateUnitVector = relativeBodyState.normalized( );

    // Pre-calculate inner product of unit vectors, and their square, for efficiency.
    double positionsInnerProduct = relativeBodyStateUnitVector.dot( stationPositionUnitVector );
    double positionsInnerProductSquared = positionsInnerProduct * positionsInnerProduct;

    // Calculate site displacement according to Eq. 7.6 of IERS Conventions 2010
    Eigen::Vector3d displacement = degreeThreeLoveNumber * stationPositionUnitVector *
            ( 2.5 * positionsInnerProductSquared * positionsInnerProduct - 1.5 * positionsInnerProduct ) + degreeThreeShidaNumber * ( 7.5 * positionsInnerProductSquared - 1.5 ) *
            ( relativeBodyStateUnitVector - positionsInnerProduct * stationPositionUnitVector );

    // Scale result according to aforementioned Eq. 7.6
    double distanceRatio = bodyEquatorialRadius / bodyDistance;
    double distanceRatioSquared = distanceRatio * distanceRatio;
    return displacement * distanceRatioSquared * distanceRatioSquared * bodyEquatorialRadius * gravitationalParameterRatio;

}


BasicTidalBodyDeformation::BasicTidalBodyDeformation(
        const std::function< Eigen::Vector6d( const double ) > deformedBodyStateFunction,
        const std::vector< std::function< Eigen::Vector6d( const double ) > >& deformingBodyStateFunctions,
        const std::function< Eigen::Quaterniond( const double ) > deformedBodyRotationFunction,
        const std::function< double( ) > gravitionalParameterOfDeformedBody,
        const std::vector< std::function< double( ) > >& gravitionalParametersOfDeformingBodies,
        const double deformedBodyEquatorialRadius,
        const std::map< int, std::pair< double, double > >& displacementLoveNumbers ):
    deformedBodyStateFunction_( deformedBodyStateFunction ),
    deformingBodyStateFunctions_( deformingBodyStateFunctions ),
    deformedBodyRotationFunction_( deformedBodyRotationFunction ),
    gravitionalParameterOfDeformedBody_( gravitionalParameterOfDeformedBody ),
    gravitionalParametersOfDeformingBodies_( gravitionalParametersOfDeformingBodies ),
    deformedBodyEquatorialRadius_( deformedBodyEquatorialRadius ),
    displacementLoveNumbers_( displacementLoveNumbers ),
    currentTime_( TUDAT_NAN )
{
    if( deformingBodyStateFunctions_.size( ) != gravitionalParametersOfDeformingBodies_.size( ) )
    {
        throw std::runtime_error( "Error in basic tidal body deformation model; different numnber of state and mass functions" );
    }

    numberOfBodies_ = deformingBodyStateFunctions_.size( );

    relativeBodyStates_.resize( numberOfBodies_);
    gravitationalParameterRatios_.resize( numberOfBodies_ );
}

Eigen::Vector3d BasicTidalBodyDeformation::calculateDisplacement(
        const double time,
        const Eigen::Vector3d& bodyFixedPosition )
{
    updateBodyProperties( time );
    Eigen::Vector3d siteDisplacement = Eigen::Vector3d::Zero( );
    Eigen::Vector3d nominalSiteUnitVector = bodyFixedPosition.normalized( );

    for( int i = 0; i < numberOfBodies_; i++ )
    {
        for( auto it : displacementLoveNumbers_ )
        {
            switch( it.first )
            {
            case 2:
                siteDisplacement += calculateDegreeTwoBasicTidalDisplacement( gravitationalParameterRatios_[ i ],
                                                                              nominalSiteUnitVector,
                                                                              relativeBodyStates_[ i ],
                                                                              deformedBodyEquatorialRadius_,
                                                                              it.second.first,
                                                                              it.second.second );
                break;
            case 3:
                siteDisplacement += calculateDegreeThreeBasicTidalDisplacement( gravitationalParameterRatios_[ i ],
                                                                                nominalSiteUnitVector,
                                                                                relativeBodyStates_[ i ],
                                                                                deformedBodyEquatorialRadius_,
                                                                                it.second.first,
                                                                                it.second.second );
                break;

            default:
                throw std::runtime_error( "Error, basic tidal displacements of degrees other than 2 or 3 not implemented" );

                break;
            }
        }
    }
    return siteDisplacement;
}

void BasicTidalBodyDeformation::updateBodyProperties( const double time )
{
    if( !( currentTime_ == time ) )
    {
        Eigen::Vector3d deformedBodyPosition = deformedBodyStateFunction_( time ).segment( 0, 3 );
        Eigen::Quaterniond currentDeformedBodyRotation = deformedBodyRotationFunction_( time );
        double currentDeformedBodyGravitationalParameter = gravitionalParameterOfDeformedBody_( );


        for( int i = 0; i < numberOfBodies_; i++ )
        {
            gravitationalParameterRatios_[ i ] = gravitionalParametersOfDeformingBodies_[ i ]( ) /
                    currentDeformedBodyGravitationalParameter;
            relativeBodyStates_[ i ] = ( currentDeformedBodyRotation * (
                                             deformingBodyStateFunctions_[ i ]( time ).segment( 0, 3 ) -
                                         deformedBodyPosition ) );
        }
        currentTime_ = time;
    }
}

} // namespace basic_astrodynamics

} // namespace tudat
