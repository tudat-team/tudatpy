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

#ifndef TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
#define TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H

#include <iostream>
#include <functional>
#include <memory>
#include <boost/lambda/lambda.hpp>

#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/basic_astro/accelerationModel.h"

namespace tudat
{

namespace gravitation
{

//! Function to compute the acceleration acting on a satellite due to tidal deformation caused by this satellite on host planet.
/*!
 * Function to compute the acceleration acting on a satellite due to tidal deformation caused by this satellite on host planet.
 * The computation is according to Lainey et al. (2007, 2009, 2012, etc.) and uses a single (degree 2) real Love number and a
 * single time lag to include the dissipation in the host planet.
 * \param relativeStateOfBodyExertingTide State of satellite w.r.t. host planet.
 * \param planetAngularVelocityVector Angular velocity vector ogf host planet
 * \param currentTidalAccelerationMultiplier Scalar multiplier that is common to all vector terms (term outside of brackets in
 * Eq. (3) of Lainey et al. (2007)).
 * \param timeLag Time lag of tidal bulge
 * \param includeDirectRadialComponent True if term independent of time lag is to be included, false otherwise
 * \return Value of acceleration vector according to tidal model
 */
Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnPlanet(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d planetAngularVelocityVector,
        const double currentTidalAccelerationMultiplier, const double timeLag,
        const bool includeDirectRadialComponent = true );

//! Function to compute the acceleration acting on a satellite due to tidal deformation caused in this satellite by host planet.
/*!
 * Function to compute the acceleration acting on a satellite due to tidal deformation caused in this satellite by host planet.
 * The computation is according to Lainey et al. (2007, 2009, 2012, etc.) and uses a single (degree 2) real Love number and a
 * single time lag to include the dissipation in the host planet.
 * \param relativeStateOfBodyExertingTide State of satellite w.r.t. host planet.
 * \param currentTidalAccelerationMultiplier Scalar multiplier that is common to all vector terms (term outside of brackets in
 * Eq. (3) of Lainey et al. (2007), modified according to tide on satellite).
 * \param timeLag Time lag of tidal bulge
 * \param includeDirectRadialComponent True if term independent of time lag is to be included, false otherwise
 * \return Value of acceleration vector according to tidal model
 */
Eigen::Vector3d computeDirectTidalAccelerationDueToTideOnSatellite(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double currentTidalAccelerationMultiplier,
        const double timeLag, const bool includeDirectRadialComponent );

//! Class to compute the acceleration due to tidal dissipation caused by/on this satellite on/by host planet.
/*!
 *  Class to compute the acceleration due to tidal dissipation caused by/on this satellite on/by host planet. Depending on the
 *  constructor input, an object uses either the tide raised on the satellite, or the tide raised on the planet.
 *  The computation is according to Lainey et al. (2007, 2009, 2012, etc.) and uses a single (degree 2) real Love number and a
 *  single time lag to include the dissipation in the host planet. The derivation of the model is given by Mignard (1981).
 *  For the tide raised on the satellite, the term that includes the angular velocity vector of the deformed body is omitted.
 *  Instead, the radial dissipation term is derived by Murray & Dermott (1999, p. 172) to be 3/4 of the tangential dissipation
 *  term, see also Caudal et al. (2017).
 */
class DirectTidalDissipationAcceleration: public basic_astrodynamics::AccelerationModel3d
{
public:

    //! Constructor for tide raised on planet
    /*!
     * Constructor for tide raised on planet
     * \param stateFunctionOfBodyExertingTide Function returning the state of the satellite
     * \param stateFunctionOfBodyUndergoingTide Function returning the state of the planet
     * \param gravitationalParameterFunctionOfBodyUndergoingTide Function returning the gravitational parameterof the satellite
     * \param angularVelocityVectorOfBodyUndergoingTide Function returning the angular velocity vector of the planet
     * \param k2LoveNumber Static k2 Love number of the planet
     * \param timeLag Time lag of tidal bulge on planet
     * \param equatorialRadiusOfBodyUndergoingTide Reference (equatorial) radius of planet associated with Love numner
     * \param includeDirectRadialComponent  True if term independent of time lag is to be included, false otherwise
     */
    DirectTidalDissipationAcceleration(
            const std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide,
            const std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide,
            const std::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTide,
            const std::function< Eigen::Vector3d( ) > angularVelocityVectorOfBodyUndergoingTide,
            const double k2LoveNumber,
            const double timeLag,
            const double equatorialRadiusOfBodyUndergoingTide,
            const bool includeDirectRadialComponent ):
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ),
        stateFunctionOfBodyUndergoingTide_( stateFunctionOfBodyUndergoingTide ),
        gravitationalParameterFunctionOfBodyUndergoingTide_( gravitationalParameterFunctionOfBodyUndergoingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( angularVelocityVectorOfBodyUndergoingTide ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ),
        equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), modelTideOnPlanet_( true ),
        explicitLibraionalTideOnSatellite_( false )
    {
        equatorialRadiusToFifthPower_ =
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_;
    }


    DirectTidalDissipationAcceleration(
            const std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide,
            const std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide,
            const std::function< double( ) > gravitationalParameterFunctionOfBodyExertingTide,
            const std::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTide,
            const double k2LoveNumber,
            const double timeLag,
            const double equatorialRadiusOfBodyUndergoingTide,
            const bool includeDirectRadialComponent ):
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ),
        stateFunctionOfBodyUndergoingTide_( stateFunctionOfBodyUndergoingTide ),
        gravitationalParameterFunctionOfBodyExertingTide_( gravitationalParameterFunctionOfBodyExertingTide ),
        gravitationalParameterFunctionOfBodyUndergoingTide_( gravitationalParameterFunctionOfBodyUndergoingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( [ = ]( ){ return Eigen::Vector3d::Constant( TUDAT_NAN ); } ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ),
        equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), modelTideOnPlanet_( false ),
        explicitLibraionalTideOnSatellite_( false )
    {
        equatorialRadiusToFifthPower_ =
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_;
    }

    //! Constructor for tide raised on satellite
    /*!
     * Constructor for tide raised on satellite
     * \param stateFunctionOfBodyExertingTide Function returning the state of the satellite
     * \param stateFunctionOfBodyUndergoingTide Function returning the state of the planet
     * \param gravitationalParameterFunctionOfBodyExertingTide Function returning the gravitational parameterof the satellite
     * \param gravitationalParameterFunctionOfBodyUndergoingTide Function returning the gravitational parameterof the planet
     * \param k2LoveNumber Static k2 Love number of the satellite
     * \param timeLag Time lag of tidal bulge on satellite
     * \param equatorialRadiusOfBodyUndergoingTide Reference (equatorial) radius of satellite associated with Love numner
     * \param includeDirectRadialComponent  True if term independent of time lag is to be included, false otherwise
     */
    DirectTidalDissipationAcceleration(
            const std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide,
            const std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide,
            const std::function< double( ) > gravitationalParameterFunctionOfBodyExertingTide,
            const std::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTide,
            const std::function< Eigen::Vector3d( ) > angularVelocityVectorOfBodyUndergoingTide,
            const double k2LoveNumber,
            const double timeLag,
            const double equatorialRadiusOfBodyUndergoingTide,
            const bool includeDirectRadialComponent ):
        stateFunctionOfBodyExertingTide_( stateFunctionOfBodyExertingTide ),
        stateFunctionOfBodyUndergoingTide_( stateFunctionOfBodyUndergoingTide ),
        gravitationalParameterFunctionOfBodyExertingTide_( gravitationalParameterFunctionOfBodyExertingTide ),
        gravitationalParameterFunctionOfBodyUndergoingTide_( gravitationalParameterFunctionOfBodyUndergoingTide ),
        angularVelocityVectorOfBodyUndergoingTide_( angularVelocityVectorOfBodyUndergoingTide ),
        k2LoveNumber_( k2LoveNumber ), timeLag_( timeLag ),
        equatorialRadiusOfBodyUndergoingTide_( equatorialRadiusOfBodyUndergoingTide ),
        includeDirectRadialComponent_( includeDirectRadialComponent ), modelTideOnPlanet_( false ),
        explicitLibraionalTideOnSatellite_( true )
    {
        equatorialRadiusToFifthPower_ =
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_ * equatorialRadiusOfBodyUndergoingTide_ *
                equatorialRadiusOfBodyUndergoingTide_;
    }

    //! Destructor
    ~DirectTidalDissipationAcceleration( ){ }

    //! Update member variables used by the acceleration model.
    /*!
    * Updates member variables used by the acceleration model.
    * Function pointers to retrieve the current values of quantities from which the
    * acceleration is to be calculated are set by constructor. This function calls
    * them to update the associated variables to their current state.
    * \param currentTime Time at which acceleration model is to be updated.
    */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        if( !( this->currentTime_ == currentTime ) )
        {
            // Update relative state
            currentRelativeState_ = stateFunctionOfBodyExertingTide_( ) - stateFunctionOfBodyUndergoingTide_( );
            double distance = currentRelativeState_.segment( 0, 3 ).norm( );
            double distanceSquared = distance * distance;
            double distanceToEighthPower = distanceSquared * distanceSquared * distanceSquared * distanceSquared;

            // Check if tide on planet or satellite is to be used
            if( modelTideOnPlanet_ )
            {
                // Retrieve/compute planet-specific models
                currentTidalAccelerationMultiplier_ =
                        - 3.0 * gravitationalParameterFunctionOfBodyUndergoingTide_( ) * equatorialRadiusToFifthPower_ /
                        distanceToEighthPower * k2LoveNumber_;
                currentAngularVelocityVectorOfBodyUndergoingTide_ = angularVelocityVectorOfBodyUndergoingTide_( );

                currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnPlanet(
                            currentRelativeState_, currentAngularVelocityVectorOfBodyUndergoingTide_,
                            currentTidalAccelerationMultiplier_, timeLag_, includeDirectRadialComponent_ );
            }
            else
            {
                if( !explicitLibraionalTideOnSatellite_ )
                {
                    // Retrieve/compute satellite-specific models
                    currentTidalAccelerationMultiplier_ =
                            - 3.0 * gravitationalParameterFunctionOfBodyExertingTide_( ) * equatorialRadiusToFifthPower_ /
                            distanceToEighthPower * k2LoveNumber_;
                    currentTidalAccelerationMultiplier_ *= gravitationalParameterFunctionOfBodyExertingTide_( ) /
                            gravitationalParameterFunctionOfBodyUndergoingTide_( );

                    currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnSatellite(
                                currentRelativeState_, currentTidalAccelerationMultiplier_,
                                timeLag_, includeDirectRadialComponent_ );
                }
                else
                {
                    // Retrieve/compute satellite-specific models
                    currentTidalAccelerationMultiplier_ =
                            - 3.0 * gravitationalParameterFunctionOfBodyExertingTide_( ) * equatorialRadiusToFifthPower_ /
                            distanceToEighthPower * k2LoveNumber_;
                    currentTidalAccelerationMultiplier_ *= gravitationalParameterFunctionOfBodyExertingTide_( ) /
                            gravitationalParameterFunctionOfBodyUndergoingTide_( );
                    currentAngularVelocityVectorOfBodyUndergoingTide_ = angularVelocityVectorOfBodyUndergoingTide_( );

                    currentAcceleration_ = computeDirectTidalAccelerationDueToTideOnPlanet(
                                currentRelativeState_, currentAngularVelocityVectorOfBodyUndergoingTide_, currentTidalAccelerationMultiplier_,
                                timeLag_, includeDirectRadialComponent_ );
                }
            }

            currentTime_ = currentTime;
        }
    }

    //! Function to retrieve the current state of satellite w.r.t. planet
    /*!
     * Function to retrieve the current state of satellite w.r.t. planet
     * \return Current state of satellite w.r.t. planet
     */
    Eigen::Vector6d getCurrentRelativeState( )
    {
        return currentRelativeState_;
    }

    //! Function to retrieve the current angular velocity vector of the planet
    /*!
     * Function to retrieve the current angular velocity vector of the planet
     * \return Current angular velocity vector of the planet
     */
    Eigen::Vector3d getCurrentAngularVelocityVectorOfBodyUndergoingTide( )
    {
        return currentAngularVelocityVectorOfBodyUndergoingTide_;
    }


    //! Function to retrieve the function returning the gravitational parameter of the planet
    /*!
     * Function to retrieve the function returning the gravitational parameter of the planet
     * \return Function returning the gravitational parameter of the planet
     */
    std::function< double( ) > getGravitationalParameterFunctionOfBodyExertingTide( )
    {
        return gravitationalParameterFunctionOfBodyExertingTide_;
    }

    //! Function to retrieve the function returning the gravitational parameter of the satellite
    /*!
     * Function to retrieve the function returning the gravitational parameter of the satellite
     * \return
     */
    std::function< double( ) > getGravitationalParameterFunctionOfBodyUndergoingTide( )
    {
        return gravitationalParameterFunctionOfBodyUndergoingTide_;
    }



    //! Function to retrieve the static k2 Love number of the planet/satellite
    /*!
     * Function to retrieve the static k2 Love number of the planet/satellite
     * \return Static k2 Love number of the planet/satellite
     */
    double getK2LoveNumber( )
    {
        return k2LoveNumber_;
    }

    //! Function to retrieve the time lag of tidal bulge on planet/satellite
    /*!
     * Function to retrieve the time lag of tidal bulge on planet/satellite
     * \return Time lag of tidal bulge on planet/satellite
     */
    double getTimeLag( )
    {
        return timeLag_;
    }

    //! Function to reset the time lag of tidal bulge on planet/satellite
    /*!
     * Function to reset the time lag of tidal bulge on planet/satellite
     * \param timeLag New time lag of tidal bulge on planet/satellite
     */
    void resetTimeLag( const double timeLag )
    {
        timeLag_ = timeLag;
    }

    //! Function to retrieve the reference (equatorial) radius of planet/satellite associated with Love numner
    /*!
     * Function to retrieve the reference (equatorial) radius of planet/satellite associated with Love numner
     * \return Reference (equatorial) radius of planet/satellite associated with Love numner
     */
    double getEquatorialRadiusOfBodyUndergoingTide( )
    {
        return equatorialRadiusOfBodyUndergoingTide_;
    }

    //! Function to retrieve whether the term independent of time lag is to be included
    /*!
     * Function to retrieve whether the term independent of time lag is to be included
     * \return True if term independent of time lag is to be included, false otherwise
     */
    bool getIncludeDirectRadialComponent( )
    {
        return includeDirectRadialComponent_;
    }

    //! Function to retrieve whether the object models tidal bulge on plane
    /*!
     * Function to retrieve whether the object models tidal bulge on plane
     * \return True if object models tidal bulge on planet, false if on satellite
     */
    bool getModelTideOnPlanet( )
    {
        return modelTideOnPlanet_;
    }

    //! Function to retrieve the scalar multiplier that is common to all vector terms in acceleration model
    /*!
     * Function to retrieve the scalar multiplier that is common to all vector terms in acceleration model
     * \return Scalar multiplier that is common to all vector terms (term outside of brackets in Eq. (3) of Lainey et al. (2007),
     *  modified according to tide on satellite if modelTideOnPlanet_ is false).
     */
    double getCurrentTidalAccelerationMultiplier( )
    {
        return currentTidalAccelerationMultiplier_;
    }

private:

    //! Current state of satellite w.r.t. planet
    Eigen::Vector6d currentRelativeState_;

    //! Current angular velocity vector of the planet
    Eigen::Vector3d currentAngularVelocityVectorOfBodyUndergoingTide_;


    //! Function returning the state of the satellite
    std::function< Eigen::Vector6d( ) > stateFunctionOfBodyExertingTide_;

    //! Function returning the state of the planet
    std::function< Eigen::Vector6d( ) > stateFunctionOfBodyUndergoingTide_;

    //! Function returning the gravitational parameter of the planet
    std::function< double( ) > gravitationalParameterFunctionOfBodyExertingTide_;

    //! Function returning the gravitational parameter of the satellite
    std::function< double( ) > gravitationalParameterFunctionOfBodyUndergoingTide_;

    //! Function returning the angular velocity vector of the planet
    std::function< Eigen::Vector3d( ) > angularVelocityVectorOfBodyUndergoingTide_;

    //! Static k2 Love number of the planet/satellite
    double k2LoveNumber_;

    //! Time lag of tidal bulge on planet/satellite
    double timeLag_;

    //! Reference (equatorial) radius of planet/satellite associated with Love numner
    double equatorialRadiusOfBodyUndergoingTide_;

    //! Reference radius to fifth power
    double equatorialRadiusToFifthPower_;

    //! True if term independent of time lag is to be included, false otherwise
    bool includeDirectRadialComponent_;

    //! True if object models tidal bulge on planet, false if on satellite
    bool modelTideOnPlanet_;

    bool explicitLibraionalTideOnSatellite_;

    //! Scalar multiplier that is common to all vector terms in acceleration model
    /*!
     *  Scalar multiplier that is common to all vector terms (term outside of brackets in Eq. (3) of Lainey et al. (2007),
     *  modified according to tide on satellite if modelTideOnPlanet_ is false).
     */
    double currentTidalAccelerationMultiplier_;

};

//! Function to retrieve all DirectTidalDissipationAcceleration from an AccelerationMap, for specific deformed/deforming bodies
/*!
 * Function to retrieve all DirectTidalDissipationAcceleration from an AccelerationMap, for specific deformed/deforming bodies
 * \param accelerationModelList Total list of acceleration models for which DirectTidalDissipationAcceleration objects are to be
 * retrieved
 * \param bodyBeingDeformed Name of body being deformed, for which DirectTidalDissipationAcceleration objects are to be
 * retrieved
 * \param bodiesCausingDeformation List of bodies causing deformation, for which DirectTidalDissipationAcceleration objects
 * are to be retrieved. If this list is empty, the selection is made only on the basis of the bodyBeingDeformed input.
 * \return All DirectTidalDissipationAcceleration from an accelerationModelList, for specific bodyBeingDeformed and
 * bodiesCausingDeformation
 */
std::vector< std::shared_ptr< DirectTidalDissipationAcceleration > > getTidalDissipationAccelerationModels(
        const basic_astrodynamics::AccelerationMap accelerationModelList, const std::string bodyBeingDeformed,
        const std::vector< std::string >& bodiesCausingDeformation );


} // namespace gravitation

} // namespace tudat

#endif // TUDAT_DIRECTTIDALDISSIPATIONACCELERATION_H
