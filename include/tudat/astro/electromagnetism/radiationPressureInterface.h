/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */
//RP-OLD

#ifndef TUDAT_RADIATIONPRESSUREINTERFACE_H
#define TUDAT_RADIATIONPRESSUREINTERFACE_H

#include <vector>
#include <iostream>

#include <functional>
#include <boost/lambda/lambda.hpp>

#include <Eigen/Geometry>
#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace electromagnetism
{

//! Calculate radiation pressure at certain distance from a source.
/*!
 *  Calculate radiation pressure at certain distance from a source, in N/m^2.
 *  \param sourcePower Total power radiated by the source (isotropically) in W.
 *  \param distanceFromSource Distance from center of (spherical) source where radiation pressure
 *  is to be calculated.
 *  \return Radiation pressure at given distance from the source.
 */
double calculateRadiationPressure( const double sourcePower, const double distanceFromSource );

//! Class in which the properties of a solar radiation pressure acceleration model are stored.
/*!
 *  Class in which the properties of a solar radiation pressure acceleration model are stored and
 *  the current radiation pressure is calculated based on the source power and geometry. The
 *  current implementation is limited to a cannonball model.
 */
class RadiationPressureInterface{
public:

    //! Constructor.
    /*!
     *  Class construtor for radiation pressure interface.
     *  \param sourcePower Function returning the current total power (in W) emitted by the source
     *  body.
     *  \param sourcePositionFunction Function returning the current position of the source body.
     *  \param targetPositionFunction Function returning the current position of the target body.
     *  \param radiationPressureCoefficient Reflectivity coefficient of the target body.
     *  \param area Reflecting area of the target body.
     *  \param occultingBodyPositions List of functions returning the positions of the bodies
     *  causing occultations (default none) NOTE: Multiple concurrent occultations may currently
     *  result in slighlty underestimted radiation pressure.
     *  \param occultingBodyRadii List of radii of the bodies causing occultations (default none).
     *  \param sourceRadius Radius of the source body (used for occultation calculations) (default 0).
     */
    RadiationPressureInterface(
            const std::function< double( ) > sourcePower,
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::function< Eigen::Vector3d( ) > targetPositionFunction,
            const double radiationPressureCoefficient,
            const double area,
            const std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions =
            std::vector< std::function< Eigen::Vector3d( ) > >( ),
            const std::vector< double > occultingBodyRadii = std::vector< double > ( ),
            const double sourceRadius = 0.0,
            const std::vector< std::string > occultingBodies = std::vector< std::string >( ) ):
        sourcePower_( sourcePower ), sourcePositionFunction_( sourcePositionFunction ),
        targetPositionFunction_( targetPositionFunction ),
        radiationPressureCoefficient_( radiationPressureCoefficient ),
        radiationPressureCoefficientFunction_( [ = ]( const double ){ return radiationPressureCoefficient; } ),
        area_( area ),
        occultingBodyPositions_( occultingBodyPositions ),
        occultingBodyRadii_( occultingBodyRadii ),
        sourceRadius_( sourceRadius ),
        occultingBodies_( occultingBodies ),
        currentRadiationPressure_( TUDAT_NAN ),
        currentSolarVector_( Eigen::Vector3d::Zero( ) ),
        currentTime_( TUDAT_NAN ){ }

    RadiationPressureInterface(
            const std::function< double( ) > sourcePower,
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::function< Eigen::Vector3d( ) > targetPositionFunction,
            std::function< double( const double ) > radiationPressureCoefficientFunction,
            const double area,
            const std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions =
            std::vector< std::function< Eigen::Vector3d( ) > >( ),
            const std::vector< double > occultingBodyRadii = std::vector< double > ( ),
            const double sourceRadius = 0.0,
            const std::vector< std::string > occultingBodies = std::vector< std::string >( )  ):
        sourcePower_( sourcePower ), sourcePositionFunction_( sourcePositionFunction ),
        targetPositionFunction_( targetPositionFunction ),
        radiationPressureCoefficient_( TUDAT_NAN ),
        radiationPressureCoefficientFunction_( radiationPressureCoefficientFunction ),
        area_( area ),
        occultingBodyPositions_( occultingBodyPositions ),
        occultingBodyRadii_( occultingBodyRadii ),
        sourceRadius_( sourceRadius ),
        occultingBodies_( occultingBodies ),
        currentRadiationPressure_( TUDAT_NAN ),
        currentSolarVector_( Eigen::Vector3d::Zero( ) ),
        currentTime_( TUDAT_NAN ){ }

    //! Destructor
    virtual ~RadiationPressureInterface( ){ }

    //! Base class function to update the current properties of radiation pressure
    /*!
     *  Base class function to update the current properties of radiation pressure. This function is nominally called by the
     *  updateInterface function, which may be overridden by derived classes.
     *  \param currentTime Time at which acceleration model is to be updated.
     */
    void updateInterfaceBase(
            const double currentTime );

    //! Function to update the current value of the radiation pressure
    /*!
     *  Function to update the current value of the radiation pressure, based on functions returning
     *  the positions of the bodies involved and the source power.
     *  \param currentTime Time at which acceleration model is to be updated.
     */
    virtual void updateInterface( const double currentTime = TUDAT_NAN );

    //! Function to return the current radiation pressure due to source at target (in N/m^2).
    /*!
     *  Function to return the current radiation pressure due to source at target (in N/m^2).
     *  \return Current radiation pressure due to source at target (in N/m^2).
     */
    double getCurrentRadiationPressure( ) const
    {
        return currentRadiationPressure_;
    }

    //! Function to return the current vector from the target to the source.
    /*!
     *  Function to return the current vector from the target to the source.
     *  \return Current vector from the target to the source.
     */
    Eigen::Vector3d getCurrentSolarVector( ) const
    {
        return currentSolarVector_;
    }

    //! Function to return the function returning the current position of the source body.
    /*!
     *  Function to return the function returning the current position of the source body.
     *  \return The function returning the current position of the source body.
     */
    std::function< Eigen::Vector3d( ) > getSourcePositionFunction( ) const
    {
        return sourcePositionFunction_;
    }

    //! Function to return the function returning the current position of the target body.
    /*!
     *  Function to return the function returning the current position of the target body.
     *  \return The function returning the current position of the target body.
     */
    std::function< Eigen::Vector3d( ) > getTargetPositionFunction( ) const
    {
        return targetPositionFunction_;
    }

    //! Function to return the reflecting area of the target body.
    /*!
     *  Function to return the reflecting area of the target body.
     *  \return The reflecting area of the target body.
     */
    double getArea( ) const
    {
        return area_;
    }

    //! Function to return the radiation pressure coefficient of the target body.
    /*!
     *  Function to return the radiation pressure coefficient of the target body.
     *  \return The radiation pressure coefficient of the target body.
     */
    double getRadiationPressureCoefficient( ) const
    {
        return radiationPressureCoefficient_;
    }

    //! Function to reset a constant radiation pressure coefficient of the target body.
    /*!
     *  Function to reset a constant radiation pressure coefficient of the target body.
     *  \param radiationPressureCoefficient The new radiation pressure coefficient of the target body.
     */
    void resetRadiationPressureCoefficient( const double radiationPressureCoefficient )
    {
        radiationPressureCoefficient_ = radiationPressureCoefficient;
        radiationPressureCoefficientFunction_ = [ = ]( const double ){ return radiationPressureCoefficient; };
    }

    //! Function to reset the function to obtain the radiation pressure coefficient of the target body.
    /*!
     *  Function to reset the function to obtain the radiation pressure coefficient of the target body.
     *  \param radiationPressureCoefficientFunction New function to obtain the radiation pressure coefficient of the target body.
     */
    void resetRadiationPressureCoefficientFunction(
            const std::function< double( const double ) > radiationPressureCoefficientFunction )
    {
        radiationPressureCoefficientFunction_ = radiationPressureCoefficientFunction;
    }

    //! Function to return the function returning the current total power (in W) emitted by the
    //! source body.
    /*!
     *  Function to return the function returning the current total power (in W) emitted by the
     *  source body.
     *  \return  The function returning the current total power emitted by the source body.
     */
    std::function< double( ) > getSourcePowerFunction( ) const
    {
        return sourcePower_;
    }

    //! Function to return the current time of interface (i.e. time of last updateInterface call).
    /*!
     *  Function to return the current time of interface (i.e. time of last updateInterface call).
     *  \return Current time of interface (i.e. time of last updateInterface call).
     */
    double getCurrentTime( )
    {
        return currentTime_;
    }

    //! Function to return the list of functions returning the positions of the bodies causing
    //! occultations
    /*!
     *  Function to return the list of functions returning the positions of the bodies causing
     *  occultations
     *  \return List of functions returning the positions of the bodies causing
     *  occultations
     */
    std::vector< std::function< Eigen::Vector3d( ) > > getOccultingBodyPositions( )
    {
        return occultingBodyPositions_;
    }

    //! Function to return the list of radii of the bodies causing occultations.
    /*!
     *  Function to return the list of radii of the bodies causing occultations
     *  \return List of radii of the bodies causing occultations
     */
    std::vector< double > getOccultingBodyRadii( )
    {
        return occultingBodyRadii_;
    }

    //! Function to return the radius of the source body.
    /*!
     *  Function to return the source radius of the target body.
     *  \return The source radius of the target body.
     */
    double getSourceRadius( )
    {
        return sourceRadius_;
    }

    std::vector< std::string > getOccultingBodies( )
    {
        return occultingBodies_;
    }



protected:

    //! Function returning the current total power (in W) emitted by the source body.
    std::function< double( ) > sourcePower_;

    //! Function returning the current position of the source body.
    std::function< Eigen::Vector3d( ) > sourcePositionFunction_;

    //! Function returning the current position of the target body.
    std::function< Eigen::Vector3d( ) > targetPositionFunction_;

    //! Radiation pressure coefficient of the target body.
    double radiationPressureCoefficient_;

    //! Function to reset a constant radiation pressure coefficient of the target body.
    std::function< double( const double ) > radiationPressureCoefficientFunction_;

    //! Reflecting area of the target body.
    double area_;

    //! List of functions returning the positions of the bodies causing occultations
    std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions_;

    //! List of radii of the bodies causing occultations.
    std::vector< double > occultingBodyRadii_;

    //! Radius of the source body.
    double sourceRadius_;

    std::vector< std::string > occultingBodies_;


    //! Current radiation pressure due to source at target (in N/m^2).
    double currentRadiationPressure_;

    //! Current vector from the target to the source.
    Eigen::Vector3d currentSolarVector_;

    //! Current time of interface (i.e. time of last updateInterface call).
    double currentTime_;
};



//! Class in which the properties of a panelled radiation pressure acceleration model are stored.
/*!
 *  Class in which the properties of a panelled radiation pressure acceleration model are stored and in which the current radiation pressure
 *  is calculated based on the source power and geometry, as well as from the boxes-and-wings model used to describe the vehicle
 *  undergoing the acceleration.
 */
class PanelledRadiationPressureInterface: public RadiationPressureInterface
{
public:


    //! Constructor.
    /*!
     *  Class construtor for panelled radiation pressure interface.
     *  \param sourcePower Function returning the current total power (in W) emitted by the source
     *  body.
     *  \param sourcePositionFunction Function returning the current position of the source body.
     *  \param targetPositionFunction Function returning the current position of the target body.
     *  \param localFrameSurfaceNormals Vector containing the functions returning the surface normal
     *  in body-fixed frame for each panel.
     *  \param emissivities Vector containing the emissivity for each panel.
     *  \param areas Vector containing the area for each panel.
     *  \param diffusionCoefficients Vector containing the diffuse reflection coefficient for each panel.
     *  \param rotationFromLocalToPropagationFrame  Vector containing the functions that return the rotation
     *  from body-fixed to propagation frame for each panel.
     *  \param occultingBodyPositions List of functions returning the positions of the bodies
     *  causing occultations (default none) NOTE: Multiple concurrent occultations may currently
     *  result in slighlty underestimted radiation pressure.
     *  \param occultingBodyRadii List of radii of the bodies causing occultations (default none).
     *  \param sourceRadius Radius of the source body (used for occultation calculations) (default 0).
     */
    PanelledRadiationPressureInterface(
            const std::function< double( ) > sourcePower,
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::function< Eigen::Vector3d( ) > targetPositionFunction,
            const std::vector< std::function< Eigen::Vector3d( const double ) > > localFrameSurfaceNormals,
            const std::vector< double > emissivities,
            const std::vector< double > areas,
            const std::vector< double > diffusionCoefficients,
            const std::function< Eigen::Quaterniond( ) > rotationFromLocalToPropagationFrame,
            const std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions =
            std::vector< std::function< Eigen::Vector3d( ) > >( ),
            const std::vector< double > occultingBodyRadii = std::vector< double > ( ),
            const double sourceRadius = 0.0 ):
        RadiationPressureInterface(
            sourcePower, sourcePositionFunction, targetPositionFunction, TUDAT_NAN, TUDAT_NAN,
            occultingBodyPositions, occultingBodyRadii, sourceRadius ),
        localFrameSurfaceNormals_( localFrameSurfaceNormals ), emissivities_( emissivities ),
        areas_( areas ), diffusionCoefficients_( diffusionCoefficients ),
        rotationFromLocalToPropagationFrame_( rotationFromLocalToPropagationFrame )
    {
        surfaceNormalsInPropagationFrame_.resize( localFrameSurfaceNormals_.size( ) );
    }



    //! Base class function to update the current properties of the panelled radiation pressure acceleration model
    /*!
     *  Base class function to update the current properties of the panelled radiation pressure acceleration model.
     * This function is nominally called by the updateInterface function, which may be overridden by derived classes.
     *  \param currentTime Time at which acceleration model is to be updated.
     */
    void updateInterface( const double currentTime = TUDAT_NAN )
    {
        updateInterfaceBase( currentTime );

        Eigen::Quaterniond rotationToPropagationFrame = rotationFromLocalToPropagationFrame_( );
        for( unsigned int i = 0; i < surfaceNormalsInPropagationFrame_.size( ); i++ )
        {
            surfaceNormalsInPropagationFrame_[ i ] =
                    rotationToPropagationFrame * localFrameSurfaceNormals_.at( i )( currentTime );
        }
    }


    //! Function to return the panels emissivity vector.
    /*!
     *  Function to return the panels emissivity vector.
     *  \return Vector containing the panels emissivities.
     */
    std::vector< double > getEmissivities( )
    {
        return emissivities_;
    }

    //! Function to return the panels area vector.
    /*!
     *  Function to return the panels area vector.
     *  \return Vector containing the panels areas.
     */
    std::vector< double > getAreas( )
    {
        return areas_;
    }

    //! Function to return the panels diffuse reflection coefficients vector.
    /*!
     *  Function to return the panels diffuse reflection coefficients vector.
     *  \return Vector containing the panels diffuse reflection coefficients.
     */
    std::vector< double > getDiffusionCoefficients( )
    {
        return diffusionCoefficients_;
    }

    //! Function to return the area of a given panel.
    /*!
     *  Function to return the area of a given panel.
     *  \param index Index of the panel whose area is to be returned.
     *  \return Panel area
     */
    double getArea( const int index ) const
    {
        return areas_[ index ];
    }

    //! Function to return the emissivity of a given panel.
    /*!
     *  Function to return the emissivity coefficient of a given panel.
     *  \param index Index of the panel whose emissivity is to be returned.
     *  \return Panel emissivity
     */
    double getEmissivity( const int index ) const
    {
        return emissivities_[ index ];
    }

    //! Function to return the diffuse reflection coefficient of a given panel.
    /*!
     *  Function to return the diffuse reflection coefficient of a given panel.
     *  \param index Index of the panel whose diffuse reflection coefficient is to be returned.
     *  \return Panel diffuse reflection coefficient
     */
    double getDiffuseReflectionCoefficient( const int index ) const
    {
        return diffusionCoefficients_[ index ];
    }

    //! Function to return the current surface normal for a given panel.
    /*!
     *  Function to return the current surface normal for a given panel.
     *  \param index Index of the panel whose surface normal is to be returned.
     *  \return Panel current surface normal
     */
    Eigen::Vector3d getCurrentSurfaceNormal( const int index ) const
    {
        return surfaceNormalsInPropagationFrame_[ index ];
    }

    //! Function to return a vector containing the surface normal expressed in propagation
    //! frame for each panel.
    /*!
     *  Function to return a vector containing the surface normal expressed in propagation
     *  frame for each panel.
     *  \return Vector of surface normals in propagation frame
     */
    std::vector< Eigen::Vector3d > getSurfaceNormalsInPropagationFrame( )
    {
        return surfaceNormalsInPropagationFrame_;
    }

    //! Function to return the total number of panels.
    /*!
     *  Function to return the total number of panels.
     *  \return Number of panels
     */
    int getNumberOfPanels( )
    {
        return surfaceNormalsInPropagationFrame_.size( );
    }


private:

    //! Vector containing the time-dependent functions that return the surface normal for each panel,
    //! in body-fixed frame.
    std::vector< std::function< Eigen::Vector3d( const double ) > > localFrameSurfaceNormals_;

    //! Vector containing the emissivitie for all panels.
    std::vector< double > emissivities_;

    //! Vector containing the area for all panels.
    std::vector< double > areas_;

    //! Vector containing the diffuse reflection coefficient for all panels.
    std::vector< double > diffusionCoefficients_;

    //! Function returning the rotation from local to propagation frame.
    std::function< Eigen::Quaterniond( ) > rotationFromLocalToPropagationFrame_;

    //! Vector containing the surface normal expressed in propagation frame for each panel.
    std::vector< Eigen::Vector3d > surfaceNormalsInPropagationFrame_;

};



//! Class in which the properties of a solar sailing radiation pressure acceleration model are stored.
/*!
 *  Class in which the properties of a solar sailing radiation pressure acceleration model are stored and in which the current solar sailing
 *  radiation pressure is calculated based on the source power and geometry, as well as on the spacecraft solar sail model.
 */
class SolarSailingRadiationPressureInterface: public RadiationPressureInterface
{
public:

    //! Constructor
    /*!
     * Class constructor for radiation pressure interface (appropriate for solar sail model)
     *  \param sourcePower Function returning the current total power (in W) emitted by the source body.
     *  \param sourcePositionFunction Function returning the current position of the source body.
     *  \param targetPositionFunction Function returning the current position of the target body.
     *  \param targetVelocityFunction Function returning the current velocity of the target body w.r.t. central body.
     *  \param area Reflecting area of the target body.
     *  \param coneAngle Function returning the current cone angle of the target body.
     *  \param clockAngle Function returning the current clock angle of the target body.
     *  \param frontEmissivityCoefficient Front emissivity coefficient of the target body.
     *  \param backEmissivityCoefficient Back emissivity coefficient of the target body.
     *  \param frontLambertianCoefficient Front Lambertian coefficient of the target body.
     *  \param backLambertianCoefficient Back Lambertian coefficient of the target body.
     *  \param reflectivityCoefficient Reflectivity coefficient of the target body.
     *  \param specularReflectionCoefficient Specular reflection coefficient of the target body.
     *  \param occultingBodyPositions List of functions returning the positions of the bodies
     *  causing occultations (default none) NOTE: Multiple concurrent occultations may currently
     *  result in slighlty underestimted radiation pressure.
     *  \param centralBodyVelocity Function returning the current velocity of the central body.
     *  \param occultingBodyRadii List of radii of the bodies causing occultations (default none).
     *  \param sourceRadius Radius of the source body (used for occultation calculations) (default 0).
     */
    SolarSailingRadiationPressureInterface(
            const std::function< double( ) > sourcePower,
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::function< Eigen::Vector3d( ) > targetPositionFunction,
            const std::function< Eigen::Vector3d( ) > targetVelocityFunction,
            const double area,
            const std::function< double( const double ) > coneAngle,
            const std::function< double( const double  ) > clockAngle,
            const double frontEmissivityCoefficient,
            const double backEmissivityCoefficient,
            const double frontLambertianCoefficient,
            const double backLambertianCoefficient,
            const double reflectivityCoefficient,
            const double specularReflectionCoefficient,
            const std::vector< std::function< Eigen::Vector3d( ) > > occultingBodyPositions =
                std::vector< std::function< Eigen::Vector3d( ) > >( ),
            const std::function< Eigen::Vector3d( ) > centralBodyVelocity =
                std::function< Eigen::Vector3d( ) >( ),
            const std::vector< double > occultingBodyRadii = std::vector< double > ( ),
            const double sourceRadius = 0.0 ):
            RadiationPressureInterface( sourcePower, sourcePositionFunction, targetPositionFunction, TUDAT_NAN, area,
            occultingBodyPositions, occultingBodyRadii, sourceRadius ),
            targetVelocityFunction_( targetVelocityFunction ),
            coneAngleFunction_( coneAngle ),
            clockAngleFunction_( clockAngle ),
            frontEmissivityCoefficient_( frontEmissivityCoefficient ),
            backEmissivityCoefficient_( backEmissivityCoefficient ),
            frontLambertianCoefficient_( frontLambertianCoefficient ),
            backLambertianCoefficient_( backLambertianCoefficient ),
            reflectivityCoefficient_( reflectivityCoefficient ),
            specularReflectionCoefficient_( specularReflectionCoefficient ),
            centralBodyVelocity_( centralBodyVelocity ){ }


    //! Base class function to update the current properties of the solar sailing radiation pressure acceleration model
    /*!
     *  Base class function to update the current properties of the solar sailing radiation pressure acceleration model.
     * This function is nominally called by the updateInterface function, which may be overridden by derived classes.
     *  \param currentTime Time at which acceleration model is to be updated.
     */
    void updateInterface( const double currentTime = TUDAT_NAN )
    {
        updateInterfaceBase( currentTime );

        // Update cone and clock angles.
        currentConeAngle_ = coneAngleFunction_( currentTime );
        currentClockAngle_ = clockAngleFunction_( currentTime );

        // Update normalised velocity vector of the spacecraft w.r.t. central body.
        currentUnitVelocityVector_ = ( targetVelocityFunction_( ) - centralBodyVelocity_( ) ).normalized( );

    }


    //! Function to return the current normalised velocity of the target.
    /*!
     *  Function to return the current normalised velocity of the target.
     *  \return Current normalised velocity vector of the target [-].
     */
    Eigen::Vector3d getCurrentVelocityVector( ) const
    {
        return currentUnitVelocityVector_;
    }

    //! Function to return the function returning the current normalised velocity of the target body.
    /*!
     *  Function to return the function returning the current normalised velocity of the target body.
     *  \return Function returning the current normalised velocity of the target body.
     */
    std::function< Eigen::Vector3d( ) > getTargetVelocityFunction( ) const
    {
        return targetVelocityFunction_;
    }

    //! Function to return the current cone angle of the target body.
    /*!
     *  Function to return the current cone angle of the target body.
     *  \return Current cone angle of the target body [rad].
     */
    double getCurrentConeAngle( ) const
    {
        return currentConeAngle_;
    }

    //! Function to return the current clock angle of the target body.
    /*!
     *  Function to return the current clock angle of the target body.
     *  \return The current clock angle of the target body.
     */
    double getCurrentClockAngle( ) const
    {
        return currentClockAngle_;
    }

    //! Function to return the cone angle function
    /*!
     *  Function to return the cone angle function
     *  \return Cone angle function.
     */
    std::function< double( const double ) > getConeAngleFunction( ) const
    {
        return coneAngleFunction_;
    }

    //! Function to return the clock angle function
    /*!
     *  Function to return the clock angle function
     *  \return Clock angle function.
     */
    std::function< double( const double ) > getClockAngleFunction( ) const
    {
        return clockAngleFunction_;
    }

    //! Function to return the front emissivity coefficient of the target body.
    /*!
     *  Function to return the front emissivity coefficient of the target body.
     *  \return Front emissivity coefficient of the target body [-].
     */
    double getFrontEmissivityCoefficient( ) const
    {
        return frontEmissivityCoefficient_;
    }

    //! Function to return the back emissivity coefficient of the target body.
    /*!
     *  Function to return the back emissivity coefficient of the target body.
     *  \return Back emissivity coefficient of the target body [-].
     */
    double getBackEmissivityCoefficient( ) const
    {
        return backEmissivityCoefficient_;
    }

    //! Function to return the front Lambertian coefficient of the target body.
    /*!
     *  Function to return the front Lambertian coefficient of the target body.
     *  \return Front Lambertian coefficient of the target body [-].
     */
    double getFrontLambertianCoefficient( ) const
    {
        return frontLambertianCoefficient_;
    }

    //! Function to return the back Lambertian coefficient of the target body.
    /*!
     *  Function to return the back Lambertian coefficient of the target body.
     *  \return Back Lambertian coefficient of the target body [-].
     */
    double getBackLambertianCoefficient( ) const
    {
        return backLambertianCoefficient_;
    }

    //! Function to return the reflectivity coefficient of the target body.
    /*!
     *  Function to return the reflectivity coefficient of the target body.
     *  \return Reflectivity coefficient of the target body  [-].
     */
    double getReflectivityCoefficient( ) const
    {
        return reflectivityCoefficient_;
    }

    //! Function to return the specular reflection coefficient of the target body.
    /*!
     *  Function to return the specular reflection coefficient of the target body.
     *  \return Specular reflection coefficient of the target body [-].
     */
    double getSpecularReflectionCoefficient( ) const
    {
        return specularReflectionCoefficient_;
    }

    //! Returns the function returning the velocity of the central body.
    /*!
     * Returns the function returning the velocity of the central body.
     *  \return Function returning the velocity of the central body.
     */
    std::function< Eigen::Vector3d( ) > getCentralBodyVelocity( )
    {
        return centralBodyVelocity_;
    }


private:

    //! Function returning the current position of the source body.
    std::function< Eigen::Vector3d( ) > targetVelocityFunction_;

    //! Function returning the current cone angle of the target body.
    std::function< double( const double  ) > coneAngleFunction_;

    //! Function returning the current clock angle of the target body.
    std::function< double( const double  ) > clockAngleFunction_;

    //! Current cone angle of the body (in rad).
    double currentConeAngle_;

    //! Current clock angle of the body (in rad).
    double currentClockAngle_;

    //! Front emissivity coefficient of the target body.
    double frontEmissivityCoefficient_;

    //! Back emissivity coefficient of the target body.
    double backEmissivityCoefficient_;

    //! Front Lambertian coefficient of the target body.
    double frontLambertianCoefficient_;

    //! Back Lambertian coefficient of the target body.
    double backLambertianCoefficient_;

    //! Reflectivity coefficient of the target body.
    double reflectivityCoefficient_;

    //! Specular reflection coefficient of the target body.
    double specularReflectionCoefficient_;

    //! Function returning the velocity of the central body.
    std::function< Eigen::Vector3d( ) > centralBodyVelocity_;

    //! Current vector of the target's normalised velocity.
    Eigen::Vector3d currentUnitVelocityVector_;

};


} // namespace electromagnetism

} // namespace tudat

#endif // TUDAT_RADIATIONPRESSUREINTERFACE_H
