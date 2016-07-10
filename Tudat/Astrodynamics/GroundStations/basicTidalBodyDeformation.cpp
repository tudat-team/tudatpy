#include <iostream>

#include "Tudat/Astrodynamics/GroundStations/basicTidalBodyDeformation.h"


namespace tudat
{

namespace site_displacements
{

//! Calculate displacement due to degree 2 tide, band-independet love and shida numbers.
Eigen::Vector3d calculateDegreeTwoInPhaseDisplacement( const double gravitationalParameterRatio,
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
    Eigen::Vector3d displacement = degreeTwoLoveNumber * stationPositionUnitVector * ( ( 3.0 *
                                                                                         positionsInnerProduct * positionsInnerProduct - 1.0 ) / 2.0 ) +
            3.0 * degreeTwoShidaNumber * positionsInnerProduct *
            ( relativeBodyStateUnitVector - positionsInnerProduct * stationPositionUnitVector );

    // Scale result according to aforementioned Eq. 7.5
    double distanceRatio = bodyEquatorialRadius / bodyDistance;
    double distanceRatioSquared = distanceRatio * distanceRatio;

    return displacement * distanceRatioSquared * distanceRatioSquared * bodyDistance * gravitationalParameterRatio;
}


//! Calculate displacement due to degree 3 tide, band-independet love and shida numbers.
Eigen::Vector3d calculateDegreeThreeInPhaseDisplacement( const double gravitationalParameterRatio,
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
    Eigen::Vector3d displacement = degreeThreeLoveNumber * stationPositionUnitVector * ( 2.5 * positionsInnerProductSquared * positionsInnerProduct -
                                                                                         1.5 * positionsInnerProduct ) + degreeThreeShidaNumber * ( 7.5 * positionsInnerProductSquared - 1.5 ) *
            ( relativeBodyStateUnitVector - positionsInnerProduct * stationPositionUnitVector );

    // Scale result according to aforementioned Eq. 7.6
    double distanceRatio = bodyEquatorialRadius / bodyDistance;
    double distanceRatioSquared = distanceRatio * distanceRatio;
    return displacement * distanceRatioSquared * distanceRatioSquared * bodyEquatorialRadius * gravitationalParameterRatio;

}

BasicTidalBodyDeformation::BasicTidalBodyDeformation(
        const boost::function< basic_mathematics::Vector6d( const double ) > deformedBodyEphemeris,
        const std::vector< boost::function< basic_mathematics::Vector6d( const double ) > >& deformingBodyEphemerides,
        const boost::function< Eigen::Quaterniond( const double ) > deformedBodyRotation,
        const std::vector< int >& maximumTidalDegrees,
        const boost::function< double( ) > gravitionalParameterOfDeformedBody,
        const std::vector< boost::function< double( ) > >& gravitionalParametersOfDeformingBodies,
        const boost::function< double( ) > deformedBodyEquatorialRadius,
        const std::vector< double >& nominalDisplacementLoveNumbers,
        const std::vector< double >& nominalDisplacementShidaNumbers ):
    deformedBodyEphemeris_( deformedBodyEphemeris ),
    deformingBodyEphemerides_( deformingBodyEphemerides ),
    deformedBodyRotation_( deformedBodyRotation ),
    maximumTidalDegrees_( maximumTidalDegrees ),
    gravitionalParameterOfDeformedBody_( gravitionalParameterOfDeformedBody ),
    gravitionalParametersOfDeformingBodies_( gravitionalParametersOfDeformingBodies ),
    deformedBodyEquatorialRadius_( deformedBodyEquatorialRadius ),
    nominalDisplacementLoveNumbers_( nominalDisplacementLoveNumbers ),
    nominalDisplacementShidaNumbers_( nominalDisplacementShidaNumbers )
{
    numberOfBodies_ = deformingBodyEphemerides.size( );
    assert( maximumTidalDegrees_.size( ) == numberOfBodies_ );
    assert( gravitionalParametersOfDeformingBodies_.size( ) == numberOfBodies_ );

    int maximumTideDegree = 0;
    for( unsigned int i = 0; i < maximumTidalDegrees_.size( ); i++ )
    {
        if( maximumTidalDegrees_[ i ] > maximumTideDegree )
        {
            maximumTideDegree = maximumTidalDegrees_[ i ];
        }
        if( maximumTidalDegrees_[ i ] < 2 )
        {
            //std::cerr<<"Warning, tide of degree lower than 2 requested when making body deformation model."<<std::endl;
        }
    }
    if( ( maximumTideDegree - 1 ) > static_cast< int >( nominalDisplacementLoveNumbers_.size( ) ) )
    {
        //std::cerr<<"Warning, Love numbers not given to sufficient order when making body deformation model"<<std::endl;
    }
    if( ( maximumTideDegree - 1 ) > static_cast< int >( nominalDisplacementShidaNumbers_.size( ) ) )
    {
        //std::cerr<<"Warning, Shida numbers not given to sufficient order when making body deformation model"<<std::endl;
    }
}

Eigen::Vector3d BasicTidalBodyDeformation::calculateSiteDisplacement(
        const double ephemerisTime, const Eigen::Vector3d& nominalSiteUnitVector )
{
    Eigen::Vector3d siteDisplacement = Eigen::Vector3d::Zero( );
    double currentGravitationalParameterRatio;
    Eigen::Vector3d currentRelativeBodyState;
    Eigen::Quaterniond currentDeformedBodyRotation = deformedBodyRotation_( ephemerisTime );// to planet-centered,fixe
    Eigen::Vector3d deformedBodyPosition = deformedBodyEphemeris_( ephemerisTime ).segment( 0, 3 );

    if( deformedBodyPosition.x( ) != deformedBodyPosition.x( ) )
    {
        std::cout<<"Def: "<<deformedBodyPosition.transpose( )<<std::endl;
    }
    for( int i = 0; i < numberOfBodies_; i++ )
    {
        currentGravitationalParameterRatio = gravitionalParametersOfDeformingBodies_[ i ]( ) /
                gravitionalParameterOfDeformedBody_( );
        currentRelativeBodyState = ( currentDeformedBodyRotation * (
                                         deformingBodyEphemerides_[ i ]( ephemerisTime ).segment( 0, 3 ) -
                                         deformedBodyPosition ) );
        if( currentRelativeBodyState.x( ) != currentRelativeBodyState.x( ) )
        {
            std::cout<<"Def2: "<<i<<" "<<deformedBodyPosition.transpose( )<<std::endl;
            std::cout<<"Def2: "<<i<<" "<<deformingBodyEphemerides_[ i ]( ephemerisTime ).segment( 0, 3 ).transpose( )<<std::endl;
            std::cout<<"Def2: "<<i<<" "<<Eigen::Matrix3d( currentDeformedBodyRotation )<<std::endl;
        }
        //std::cout<<"Bod :"<<i<<" "<<currentRelativeBodyState.transpose( )<<std::endl;

        for( int j = 2; j <= maximumTidalDegrees_[ i ]; j++ )
        {
            //std::cout<<"Deg. :"<<maximumTidalDegrees_[ i ]<<" "<<j<<std::endl;
            switch( j )
            {
            case 2:
                siteDisplacement += calculateDegreeTwoInPhaseDisplacement( currentGravitationalParameterRatio,
                                                                           nominalSiteUnitVector,
                                                                           currentRelativeBodyState,
                                                                           deformedBodyEquatorialRadius_( ),
                                                                           nominalDisplacementLoveNumbers_[ 0 ],
                                                                           nominalDisplacementShidaNumbers_[ 0 ] );
                break;
            case 3:
                siteDisplacement += calculateDegreeThreeInPhaseDisplacement( currentGravitationalParameterRatio,
                                                                             nominalSiteUnitVector,
                                                                             currentRelativeBodyState,
                                                                             deformedBodyEquatorialRadius_( ),
                                                                             nominalDisplacementLoveNumbers_[ 1 ],
                                                                             nominalDisplacementShidaNumbers_[ 1 ] );
                break;

            default:
                std::cerr<<"Warning, tidal displacements of degree > 3 ( or < 2) not yet implemented"<<std::endl;

                break;
            }
        }
    }
    return siteDisplacement;
}

}

}
