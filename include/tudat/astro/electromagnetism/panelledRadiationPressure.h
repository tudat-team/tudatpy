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

#ifndef TUDAT_PANELLEDRADIATIONPRESSURE_H
#define TUDAT_PANELLEDRADIATIONPRESSURE_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>

#include <boost/function.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"



namespace tudat
{

namespace electromagnetism
{

//! Function to calculate radiation pressure force on a single partially reflecting panel
/*!
 *  Function to calculate radiation pressure force on a single partially reflecting panel. Model for calculation is
 *  according to Eq. (3.72) of Montenbruck and Gill (2000)
 *  \param normalizedVectorToSource Current normalized vector from accelerated body to source
 *  \param panelSurfaceNormal Panel surface normal vector, in the same frame as the normalizedVectorToSource vector
 *  \param panelArea Area of panel that is considered
 *  \param panelEmissivitty Emissivity of panel that is considered
 *  \param panelDiffuseReflectionCoefficient Diffuse reflection coefficient of panel that is considered
 *  \return The radiation pressure force on a single partially reflecting panel
 */
Eigen::Vector3d computeSinglePanelNormalizedRadiationPressureForce(
        const Eigen::Vector3d& normalizedVectorToSource, const Eigen::Vector3d& panelSurfaceNormal,
        const double panelArea, const double panelEmissivitty, const double panelDiffuseReflectionCoefficient );

//! Class for calculating the radiation pressure acceleration on a panelled body
/*!
 *  Class for calculating the radiation pressure acceleration on a panelled body, with the force due to each panel calculated from
 *  its area, orientation and emissivity. The emissivity determines the fraction that is absorbed (modelled as a force in line with
 *  the vector to th source) and one from reflection (modelled as a force normal to the panel).
 */
class PanelledRadiationPressureAcceleration: public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:

    //! Constructor for setting up the acceleration model, with separate input variables for all required data.
    /*!
     *  Constructor for setting up the acceleration model, with separate input variables for all required data.
     *  \param sourcePositionFunction Function providing current position for the source body (i.e. the body from which the
     *  radiation originates)
     *  \param acceleratedBodyPositionFunction Function providing current position for the body on which the force is acting
     *  \param panelEmissivittyFunctions Vector of functions returning emissivities for all panels
     *  \param panelDiffuseReflectionCoefficientFunctions  Vector of functions returning diffuse reflection coefficients
     *  for all panels
     *  \param panelSurfaceNormalFunctions Vector of functions returning panel surface normal function, in the same frame as the
     *  position functions of the accelerated and radiating bodies.
     *  \param panelAreaFunctions Vector of functions returning areas for all panels
     *  \param radiationPressureFunction Function returning the current radiation pressure (i.e. incident flux, in W/m^2, divided by
     *  speed of light)
     *  \param massFunction Function returning the current mass of the body being accelerated
     */
    PanelledRadiationPressureAcceleration(
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::function< Eigen::Vector3d( ) > acceleratedBodyPositionFunction,
            const std::vector< std::function< double( ) > >& panelEmissivittyFunctions,
            const std::vector< std::function< double( ) > >& panelDiffuseReflectionCoefficientFunctions,
            const std::vector< std::function< Eigen::Vector3d( ) > >& panelSurfaceNormalFunctions,
            const std::vector< std::function< double( ) > >& panelAreaFunctions,
            const std::function< double( ) > radiationPressureFunction,
            const std::function< double( ) > massFunction ):
        sourcePositionFunction_( sourcePositionFunction ),
        acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ),
        radiationPressureFunction_( radiationPressureFunction ),
        panelEmissivittyFunctions_( panelEmissivittyFunctions ),
        panelDiffuseReflectionCoefficientFunctions_( panelDiffuseReflectionCoefficientFunctions ),
        panelSurfaceNormalFunctions_( panelSurfaceNormalFunctions ),
        panelAreaFunctions_( panelAreaFunctions ),
        massFunction_( massFunction )
    {
        // Set number of panels and resize vector of panel accelerations.
        currentPanelAccelerations_.resize( panelEmissivittyFunctions_.size( ) );
        numberOfPanels_ = panelEmissivittyFunctions.size( );
    }

    //! Constructor for setting up the acceleration model, with input RadiationPressureInterface (and massFunction)
    /*!
     *  Constructor for setting up the acceleration model, with input RadiationPressureInterface (and massFunction). The massFunction
     *  is only used to convert force to acceleration. All radiation pressure properties are taken from the
     *  RadiationPressureInterface object.
     *  \param radiationPressureInterface Object in which radiation pressure properties of accelerated body due to radiation
     *  from body causing acceleration is stored.
     *  \param massFunction Function returning the current mass of the body being accelerated
     */
    PanelledRadiationPressureAcceleration(
            const std::shared_ptr< PanelledRadiationPressureInterface > radiationPressureInterface,
            const std::function< double( ) > massFunction );

    //! Update member variables used by the radiation pressure acceleration model.
    /*!
     * Updates member variables used by the acceleration model.
     * This function evaluates all dependent variables to the 'current'
     * values of these parameters. Only these current values, not the function-pointers are then
     * used by the getAcceleration( ) function.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Calculate current normal vector and distance from accelerated body to source
            currentNormalizedVectorToSource_ = ( sourcePositionFunction_( ) - acceleratedBodyPositionFunction_( ) );
            currentDistanceToSource_ = currentNormalizedVectorToSource_.norm( );
            currentNormalizedVectorToSource_.normalize( );

            // Retrieve current radiation pressure and mass
            currentRadiationPressure_ = radiationPressureFunction_( );
            currentMass_ = massFunction_( );

            // Iterate over all panels and calculate acceleration due to radiation pressure on panel
            currentAcceleration_.setZero( );
            if( currentRadiationPressure_ > 0.0 )
            {
                for( int i = 0; i < numberOfPanels_; i++ )
                {
                    currentPanelAccelerations_[ i ] = currentRadiationPressure_ / massFunction_( ) * computeSinglePanelNormalizedRadiationPressureForce(
                                currentNormalizedVectorToSource_ , panelSurfaceNormalFunctions_[ i ]( ),  panelAreaFunctions_[ i ]( ),
                                panelEmissivittyFunctions_[ i ]( ), panelDiffuseReflectionCoefficientFunctions_[ i ]( ) );
                    currentAcceleration_ += currentPanelAccelerations_[ i ];
                }
            }
            this->currentTime_ = currentTime;
        }
    }

    //! Returns the function returning the current mass of the body being accelerated.
    /*!
     *  Returns the function returning the current mass of the body being accelerated.
     *  \return Function returning the current mass of the body being accelerated.
     */
    std::function< double( ) > getMassFunction( )
    {
        return massFunction_;
    }

    //! Returns the current normalized vector from the accelerated body to the source body
    /*!
     *  Returns the current normalized vector from the accelerated body to the source body, as set by the last call to the updateMembers function
     *  \return The current normalized vector from the accelerated body to the source body
     */
    Eigen::Vector3d getCurrentVectorToSource( )
    {
        return currentNormalizedVectorToSource_;
    }

    //! Returns the current distance from the accelerated body to the source body
    /*!
     *  Returns the current distance from the accelerated body to the source body, as set by the last call to the updateMembers function
     *  (i.e. norm of the current vector from accelerated body to the source)
     *  \return The current distance from the accelerated body to the source body
     */
    double getCurrentDistanceToSource( )
    {
        return currentDistanceToSource_;
    }

    //! Returns the current acceleration due to the radiation pressure on a single panel
    /*!
     *  Returns the current acceleration due to the radiation pressure on a single panel, as calculated by the last call to
     *  the updateMembers function
     *  \param panelIndex Index of panel for which acceleration is to be retrieved
     *  \return The current acceleration due to the radiation pressure on a single panel
     */
    Eigen::Vector3d getCurrentPanelAcceleration( const int panelIndex )
    {
        return currentPanelAccelerations_[ panelIndex ];
    }

    //! Returns the current surface normal in propagation frame of a single panel
    /*!
     *   Returns the current surface normal in propagation frame of a single panel, as calculated by the last call to
     *  the updateMembers function
     *  \param panelIndex Index of panel for which acceleration is to be retrieved
     *  \return The current surface normal in propagation frame of a single panels
     */
    Eigen::Vector3d getCurrentPanelSurfaceNormalInPropagationFrame( const int panelIndex )
    {
        return panelSurfaceNormalFunctions_[ panelIndex ]( );
    }

    //! Returns the current radiation pressure at the accelerated body
    /*!
     *  Returns the current radiation pressure at the accelerated body, as set by the last call to the updateMembers function
     *  (i.e. incident flux, in W/m^2, divided by speed of light)
     *  \return The current radiation pressure at the accelerated body
     */
    double getCurrentRadiationPressure( )
    {
        return currentRadiationPressure_;
    }

    //! Returns the current mass of the accelerated body
    /*!
     *  Returns the current mass of the accelerated body, as set by the last call to the updateMembers function
     *  \return The current mass of the accelerated body
     */
    double getCurrentMass( )
    {
        return currentMass_;
    }



private:

    //! Function pointer returning position of source.
    /*!
     *  Function pointer returning position of source.
     */
    std::function< Eigen::Vector3d( ) > sourcePositionFunction_;

    //! Function pointer returning position of accelerated body.
    /*!
     *  Function pointer returning position of accelerated body.
     */
    std::function< Eigen::Vector3d( ) > acceleratedBodyPositionFunction_;

    //! Function pointer returning radiation pressure.
    /*!
     *  Function pointer returning radiation pressure.
     */
    std::function< double( ) > radiationPressureFunction_;

    //! Vector of functions returning emissivities for all panels
    std::vector< std::function< double( ) > > panelEmissivittyFunctions_;

    //! Vector of functions returning diffuse reflection coefficients for all panels
    std::vector< std::function< double( ) > > panelDiffuseReflectionCoefficientFunctions_;

    //! Vector of functions returning panel surface normal functions
    /*!
     *  Vector of functions returning panel surface normal function, in the same frame as the
     *  position functions of the accelerated and radiating bodies.
     */
    std::vector< std::function< Eigen::Vector3d( ) > > panelSurfaceNormalFunctions_;

    //! Vector of functions returning areas for all panels
    std::vector< std::function< double( ) > > panelAreaFunctions_;

    //! Function pointer returning mass of accelerated body.
    /*!
     *  Function pointer returning mass of accelerated body.
     */
    std::function< double( ) > massFunction_;

    //! The current accelerations due to the radiation pressure on all single panel
    /*!
     *  The current accelerations due to the radiation pressure on all single panel, as calculated by the last call to
     *  the updateMembers function
     */
    std::vector< Eigen::Vector3d > currentPanelAccelerations_;

    //! Number of panels used in acceleration model
    int numberOfPanels_;

    //! Current radiation pressure.
    /*!
     *  Current radiation pressure, as set by the last call to the updateMembers function
     */
    double currentRadiationPressure_;

    //! Current vector from accelerated body to source.
    /*!
     *  Current vector from accelerated body to source, as calculated by the last call to the updateMembers function
     */
    Eigen::Vector3d currentNormalizedVectorToSource_;

    //! The current distance from the accelerated body to the source body
    /*!
     *  The current distance from the accelerated body to the source body, as set by the last call to the updateMembers function
     *  (i.e. norm of the current vector from accelerated body to the source)
     */
    double currentDistanceToSource_;

    //! Current mass of accelerated body.
    /*!
     *  Current mass of accelerated body, as set by the last call to the updateMembers function
     */
    double currentMass_;


};

}

}
#endif // TUDAT_PANELLEDRADIATIONPRESSURE_H
