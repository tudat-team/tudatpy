#ifndef NOMINALGROUNDSTATIONSTATE_H
#define NOMINALGROUNDSTATIONSTATE_H

#include <vector>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"


namespace tudat
{

namespace ground_stations
{

//! Enum defininf available types of body deformations.
/*!
 *  Enum defininf available types of body deformations. Used for defining type of both body deformation and gravity field variation.
 */
enum BodyDeformationTypes
{
    basic_solid_body = 0,
};

Eigen::Vector3d dummyPositionVariationFunction( const double time );

//! Class containing settings for ground station position variations
/*!
 *  Class containing settings for ground station position variations, including functions for obtaining these variations and interpolation
 *  settings for (re)generating a ground station position variation function. Used as a member variable of NominalGroundStationState
 *  for updating position variation functions of ground station (if applicable & necessary).
 */
struct DirectGroundStationPositionVariationSettings
{
public:


    DirectGroundStationPositionVariationSettings(
            const double initialTime, const double finalTime, const double timeStep,
            const std::vector< boost::function< Eigen::Vector3d( const double ) > >& positionVariations =
            ( std::vector< boost::function< Eigen::Vector3d( const double ) > >( ) ),
            const std::map< int, boost::function< Eigen::Vector3d( ) > >& meanVariationsToSubtract =
            ( std::map< int, boost::function< Eigen::Vector3d( ) > >( ) ) ):
        initialTime_( initialTime ), finalTime_( finalTime ), timeStep_( timeStep ),
        positionVariations_( positionVariations ), meanVariationsToSubtract_( meanVariationsToSubtract )
    {
        subtractMeanValues_.reserve( positionVariations_.size( ) );

        for( unsigned i = 0; i < positionVariations_.size( ); i++ )
        {
            if( meanVariationsToSubtract_.count( i ) == 1 )
            {
                subtractMeanValues_[ i ] = 1;
            }
            else
            {
                subtractMeanValues_[ i ] = 0;
            }
        }
    }

    double initialTime_;
    double finalTime_;
    double timeStep_;
    std::vector< boost::function< Eigen::Vector3d( const double ) > > positionVariations_;

    std::vector< bool > subtractMeanValues_;

    std::map< int, boost::function< Eigen::Vector3d( ) > > meanVariationsToSubtract_;

    boost::function< Eigen::Vector3d( const double ) > createInterpolatedPositionVariationFunction( );
};

//! Class containing for storing of and calculations on the body-fixed position of a ground station
/*!
 *  Class containing for storing of and calculations on the body-fixed position of a ground station
 */
class NominalGroundStationState
{
public:
    //! Constructor taking cartesian position data.
    /*!
     *  Constructor taking cartesian position data.
     */
    NominalGroundStationState(
            const Eigen::Vector3d stationCartesianPosition,
            const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface = boost::shared_ptr< basic_astrodynamics::BodyShapeModel >( ),
            const double referenceJulianYear = basic_astrodynamics::JULIAN_DAY_ON_J2000 / physical_constants::JULIAN_YEAR_IN_DAYS,
            const bool setTransformation = 1 );

    virtual ~NominalGroundStationState( ){ }

    //! Function to obtain the Cartesian position of the state in the local frame at a given time.
    /*!
     *  Function to obtain the Cartesian position of the state in the local (body-fixed, not topocentrix) frame at a given time.
     *  Adds all position variations to the nominal state (at the requested time) and returns the state.
     *  \param secondsSinceEpoch Secons since reference epoch at which the position is to be retrieved.
     *  \param inputReferenceEpoch Reference epoch julian day
     *  \return Cartiesian position of station in local frame at requested time.
     */
    Eigen::Vector3d getCartesianPositionInTime(
            const double secondsSinceEpoch,
            const double inputReferenceEpoch = basic_astrodynamics::JULIAN_DAY_ON_J2000 );

    //! Function to return the nominal cartesian (unperturbed) position of the station
    /*!
     *  Function to return the nominal cartesian (unperturbed, i.e. not including linear drift, eccentricity, tidal variations, etc.)
     *  position of the station in body-fixed system. Typically used for applications require moderate to low precision of ground station position.
     *  \return Unperturbed position of the station
     */
    Eigen::Vector3d getNominalCartesianPosition( )
    {
        return cartesianPosition_;
    }

    //! Function to return the nominal spherical (unperturbed) position of the station
    /*!
     *  Function to return the nominal spherical (unperturbed, i.e. not including linear drift, eccentricity, tidal variations, etc.)
     *  position of the station in body-fixed system.
     *  \return Unperturbed position of the station
     */
    Eigen::Vector3d getNominalSphericalPosition( )
    {
        return sphericalPosition_;
    }


    //! Function to return current rotation from body-fixed to topocentric frame
    /*!
     *  Function to return current rotation from body-fixed to topocentric frame, with topocentric xyz-axes in ENU (Earth-North-Up) order.
     *  Note that no relativistic rescaling of vectors is performed, only a rotation is applied.
     *  \param time Time at which rotation quaternion is to be evaluated.
     *  \return Current rotation from body-fixed to topocentrix frame
     */
    Eigen::Quaterniond getRotationFromBodyFixedToTopocentricFrame( const double time )
    {        
        return bodyFixedToTopocentricFrameRotation_;
    }

    //! Function to return a set of unit vectors of 'geocentric-topocentric' frame.
    /*!
     *  Function to return a set of unit vectors of 'geocentric-topocentric' frame, expressed in body-fixed frame, in ENU (Earth-North-Up) order.
     *  \return Unit vectors of 'geocentric-topocentric' frame.
     */
    std::vector< Eigen::Vector3d > getEastNorthRadialGeocentricUnitVectors( )
    {
        return geocentricUnitVectors_;
    }

    Eigen::Vector3d getGeocentricUnitVector( int index )
    {
        return geocentricUnitVectors_[ index ];
    }


    void setPositionVariationsUpdateFunctions(
            boost::shared_ptr< DirectGroundStationPositionVariationSettings > settings )
    {
        directPositionVariationSettings_  = settings;
        updatePositionVariations( );
    }

    boost::shared_ptr< DirectGroundStationPositionVariationSettings > getGroundStationPositionVariationSettings( )
    {
        return directPositionVariationSettings_;
    }

    void updatePositionVariations( );

    virtual void resetGroundStationPositionAtEpoch( const Eigen::Vector3d cartesianPosition );

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > getBodySurface( )
    {
        return bodySurface_;
    }

    double getReferenceJulianYear( )
    {
        return referenceJulianYear_;
    }

protected:
    void setTransformationAndUnitVectors( );

    //! Cartesian position of station
    /*!
     *  Cartesian position of station, without variations (linear drift, eccentricity, tides, etc. ), in the body-fixed frame.
     */
    Eigen::Vector3d cartesianPosition_; //(x,y,z)

    //! Spherical position of station
    /*!
     *  Spherical position of station, without variations (linear drift, eccentricity, tides, etc. ), in the body-fixed frame.
     *  The order of the entries is: radius, colatitude, longitude.
     */
    Eigen::Vector3d sphericalPosition_; //(radius, co-geocentric latitude, longitude)

    //! Name of ground station of which this object is a member
    /*!
     *  Name of ground station of which this object is a member, i.e. the station of which this object describes the state.
     */
    std::string siteId_;

    //! Set of unit vectors of 'geocentric-topocentric' frame.
    /*!
     *  Set of unit vectors of 'geocentric-topocentric' frame, expressed in body-fixed frame, in ENU (Earth-North-Up) order.
     *  The up-vector is generated by drawing a line from teh center of the body to the (nominal) position of the station.
     *  For a sphere, this coincides with the true topocentric frame, for which the East-North vectors span a local tangent plane.
     *  The transformation from body-fixed to 'true' topocentric frame is described by the bodyFixedToTopocentricFrameRotation_ variable.
     */
    std::vector< Eigen::Vector3d > geocentricUnitVectors_;

    //! Rotation from body-fixed to topocentrix frame
    /*!
     *  Rotation from body-fixed to topocentrix frame, with topocentric xyz-axes in ENU (Earth-North-Up) order.
     *  In this frame the East-North vectors span a local tangent plane. The frame is time-independent and is based on the cartesianPosition_
     *  member, without variations.
     */
    Eigen::Quaterniond bodyFixedToTopocentricFrameRotation_;


    //! Reference year (01-01, 00:00:00) from which linear velocity displacement is calculated
    /*!
     *  Reference year (01-01, 00:00:00) from which linear velocity displacement is calculated.
     */
    double referenceJulianYear_;

    boost::function< Eigen::Vector3d( const double ) > positionVariationsFunction_;

    boost::shared_ptr< DirectGroundStationPositionVariationSettings > directPositionVariationSettings_;

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodySurface_;
};

std::vector< Eigen::Vector3d > getGeocentricLocalUnitVectors( const double geocentricLatitude,
                                                              const double geocentricLongitude );
/*!
 *  Function to calculate the rotation from a body-fixed to a topocentric frame at a given point (localPoint), based on the geocentric
 *  latitude and longitude (i.e. angles based on line from center of body to localPoint. The topocentric calculation differs per type of
 *  body shape wrt which the point is located. The rotation is calculated by concatenating the unit vectors of the topocentric frame,
 *  expresse in the body-fixed frame.
 *  \param bodyShapeModel Shape model of body on which point is located.
 *  \param geocentricLatitude Geocentric latitude of point, i.e. angle between equatorial planet and line from center of body to localPoint
 *  \param geocentricLongitude Geocentric longitude of point, i.e. angle between positive x-axis and line from center of body to localPoint,
 *  projected on equatorial plane, measured positive in +y direction
 *  \param localPoint Cartesian position of point in body-fixed frame.
 *  \return The rotation from the body-fixed to the topocentric frame at localPoint.
 */
Eigen::Quaterniond getRotationQuaternionFromBodyFixedToTopocentricFrame(
        const boost::shared_ptr< basic_astrodynamics::BodyShapeModel > bodyShapeModel,
        const double geocentricLatitude,
        const double geocentricLongitude,
        const Eigen::Vector3d localPoint );
}

}
#endif // NOMINALGROUNDSTATIONSTATE_H
