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

} // namespace electromagnetism

} // namespace tudat

#endif // TUDAT_RADIATIONPRESSUREINTERFACE_H
