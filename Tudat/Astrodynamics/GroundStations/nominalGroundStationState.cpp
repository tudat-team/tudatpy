#include <boost/assign/list_of.hpp>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/oblateSpheroidBodyShapeModel.h"
#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"

namespace tudat
{

Eigen::Vector3d dummyPositionVariationFunction( const double time )
{
    return Eigen::Vector3d::Zero( );
}

boost::function< Eigen::Vector3d( const double ) >
DirectGroundStationPositionVariationSettings::createInterpolatedPositionVariationFunction( )
{
    using namespace tudat::interpolators;
    timeStep_ = 600.0;

    std::map< double, Eigen::Vector3d > positionHistory;

    std::vector< Eigen::Vector3d > currentMeanValuesToSubtract;
    currentMeanValuesToSubtract.reserve( positionVariations_.size( ) );
    for( unsigned int i = 0; i < positionVariations_.size( ); i++ )
    {
        if( subtractMeanValues_[ i ] )
        {
            if( meanVariationsToSubtract_.count( i ) > 0 )
            {
                currentMeanValuesToSubtract[ i ] = meanVariationsToSubtract_[ i ]( );
            }
            else
            {
                std::cerr<<"Error when updating ground station variations, expectec mean value function"<<std::endl;
            }
        }
    }

    double time = initialTime_;
    while( time <= finalTime_ )
    {
        positionHistory[ time ] = Eigen::Vector3d::Zero( );

        for( unsigned int i = 0; i < positionVariations_.size( ); i++ )
        {
            positionHistory[ time ] += positionVariations_[ i ]( time );
            if( subtractMeanValues_[ i ] )
            {
                positionHistory[ time ] -= currentMeanValuesToSubtract[ i ];
            }
        }
        time += timeStep_;
    }

    boost::shared_ptr< CubicSplineInterpolator< double, Eigen::Vector3d > > positionInterpolator =
            boost::make_shared< CubicSplineInterpolator< double, Eigen::Vector3d > >( positionHistory );

    return boost::bind( static_cast< Eigen::Vector3d( CubicSplineInterpolator< double, Eigen::Vector3d >::* )( const double )>
                        ( &OneDimensionalInterpolator< double, Eigen::Vector3d >::interpolate ), positionInterpolator, _1 );
}


NominalGroundStationState::NominalGroundStationState(
        const Eigen::Vector3d stationCartesianPosition,
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface,
        const double referenceJulianYear, const bool setTransformation ):
    bodySurface_( bodySurface ),
    positionVariationsFunction_( boost::lambda::constant( Eigen::Vector3d::Zero( ) ) ),
    referenceJulianYear_( referenceJulianYear )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    cartesianPosition_ = stationCartesianPosition;
    sphericalPosition_ = convertCartesianToSpherical( cartesianPosition_ );

    if( setTransformation )
    {
        setTransformationAndUnitVectors( );
    }
}


void NominalGroundStationState::updatePositionVariations( )
{
    if( directPositionVariationSettings_ != NULL )
    {
        positionVariationsFunction_ = directPositionVariationSettings_->createInterpolatedPositionVariationFunction( );
    }
}

Eigen::Vector3d NominalGroundStationState::getCartesianPositionInTime(
        const double secondsSinceEpoch,
        const double inputReferenceEpoch )
{
    //NOTE.: Calculation introduces very minor error julian reference year, not calendar yer is used (sub-0.1 mm for most stations error)
    double yearsSinceReference = ( ( basic_astrodynamics::JULIAN_DAY_ON_J2000 - inputReferenceEpoch ) /
                                   physical_constants::JULIAN_YEAR_IN_DAYS - referenceJulianYear_ ) +
            secondsSinceEpoch / ( physical_constants::JULIAN_YEAR_IN_DAYS * physical_constants::JULIAN_DAY );

    return cartesianPosition_ + siteEccentricity_( secondsSinceEpoch ) +
            positionVariationsFunction_( secondsSinceEpoch );
}

void NominalGroundStationState::resetGroundStationPositionAtEpoch( const Eigen::Vector3d cartesianPosition )
{
    using namespace coordinate_conversions;
    using mathematical_constants::PI;

    cartesianPosition_ = cartesianPosition;
    sphericalPosition_ = convertCartesianToSpherical( cartesianPosition_ );

    setTransformationAndUnitVectors( );
}

void NominalGroundStationState::setTransformationAndUnitVectors( )
{
    geocentricUnitVectors_ = getGeocentricLocalUnitVectors( getLatitude( ), getLongitude( ) );

    bodyFixedToTopocentricFrameRotation_ = getRotationQuaternionFromBodyFixedToTopocentricFrame(
                bodySurface_, getLatitude( ), getLongitude( ), cartesianPosition_  );
}

std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors( const double geocentricLatitude,
                                                              const double geocentricLongitude )
{
    Eigen::Matrix3d toPlanetFixedFrameMatrix =
            Eigen::Matrix3d( reference_frames::getEnuLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                                 geocentricLongitude, geocentricLatitude ) );

    std::vector< Eigen::Vector3d > geocentricUnitVectors;
    geocentricUnitVectors.reserve( 3 );
    geocentricUnitVectors[ 0 ] = toPlanetFixedFrameMatrix.block( 0, 0, 3, 1 );
    geocentricUnitVectors[ 1 ] = toPlanetFixedFrameMatrix.block( 0, 1, 3, 1 );
    geocentricUnitVectors[ 2 ] = toPlanetFixedFrameMatrix.block( 0, 2, 3, 1 );
    return geocentricUnitVectors;
}


//! Function to calculate the rotation from a body-fixed to a topocentric frame.
Eigen::Quaterniond getRotationQuaternionFromBodyFixedToTopocentricFrame(
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodyShapeModel,
        const double geocentricLatitude,
        const double geocentricLongitude,
        const Eigen::Vector3d localPoint )
{
    // Declare unit vectors of topocentric frame, to be calculated.
    std::vector< Eigen::Vector3d > topocentricUnitVectors;

    bool isSurfaceModelRecognized = 1;

    // Identify type of body shape model
    if( boost::dynamic_pointer_cast< basic_astrodynamics::SphericalBodyShapeModel >( bodyShapeModel ) != NULL )
    {
        // For a sphere the topocentric and geocentric frames are equal.
        topocentricUnitVectors = getGeocentricLocalUnitVectors(
                    geocentricLatitude, geocentricLongitude );
    }
    else if( boost::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel ) != NULL )
    {
        boost::shared_ptr< basic_astrodynamics::OblateSpheroidBodyShapeModel > oblateSphericalShapeModel =
                boost::dynamic_pointer_cast< basic_astrodynamics::OblateSpheroidBodyShapeModel >( bodyShapeModel );

        // Calculate geodetic latitude.
        double flattening = oblateSphericalShapeModel->getFlattening( );
        double equatorialRadius = oblateSphericalShapeModel->getEquatorialRadius( );
        double geodeticLatitude = basic_astrodynamics::calculateGeodeticLatitude( localPoint, equatorialRadius, flattening, 1.0E-4 );

        // Calculte unit vectors of topocentric frame.
        topocentricUnitVectors = getGeocentricLocalUnitVectors( geodeticLatitude, geocentricLongitude );
    }    
    else
    {
        isSurfaceModelRecognized = 0;
        std::cerr<<"Error when making transformation to topocentric frame, shape model not recognized"<<std::endl;
    }

    // Create rotation matrix

    Eigen::Matrix3d bodyFixedToTopocentricFrame;

    if( isSurfaceModelRecognized == 1 )
    {
        bodyFixedToTopocentricFrame.block( 0, 0, 1, 3 ) = topocentricUnitVectors[ 0 ].transpose( );
        bodyFixedToTopocentricFrame.block( 1, 0, 1, 3 ) = topocentricUnitVectors[ 1 ].transpose( );
        bodyFixedToTopocentricFrame.block( 2, 0, 1, 3 ) = topocentricUnitVectors[ 2 ].transpose( );
    }
    else
    {
        bodyFixedToTopocentricFrame = Eigen::Matrix3d::Identity( );
    }

    // Convert to quaternion and return.
    return Eigen::Quaterniond( bodyFixedToTopocentricFrame );
}

}

