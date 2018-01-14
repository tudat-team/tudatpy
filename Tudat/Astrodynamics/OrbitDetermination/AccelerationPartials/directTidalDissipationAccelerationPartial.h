/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIRECTTIDALDISSIPATIONACCELERATIONPARTIAL_H
#define TUDAT_DIRECTTIDALDISSIPATIONACCELERATIONPARTIAL_H

#include "Tudat/Astrodynamics/Gravitation/directTidalDissipationAcceleration.h"

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial derivative of direct tidal acceleration due to tide on planet w.r.t. position of satellite
/*!
 * Function to compute partial derivative of direct tidal acceleration due to tide on planet w.r.t. position of satellite (equal
 * to -1 * partial w.r.t. position of planet)
 * \param relativeStateOfBodyExertingTide State of satellite w.r.t. host planet.
 * \param planetAngularVelocityVector Angular velocity vector ogf host planet
 * \param currentTidalAccelerationMultiplier Scalar multiplier that is common to all vector terms (term outside of brackets in
 * Eq. (3) of Lainey et al. (2007)).
 * \param timeLag Time lag of tidal bulge
 * \param includeDirectRadialComponent True if term independent of time lag is to be included, false otherwise
 * \return Partial derivative of direct tidal acceleration due to tide on planet w.r.t. position of satellite
 */
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnPlanetWrtPosition(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d planetAngularVelocityVector,
        const double currentTidalAccelerationMultiplier, const double timeLag, const bool includeDirectRadialComponent );

//! Function to compute partial derivative of direct tidal acceleration due to tide on planet w.r.t. velocity of satellite
/*!
 * Function to compute partial derivative of direct tidal acceleration due to tide on planet w.r.t. velocity of satellite (equal
 * to -1 * partial w.r.t. velocity of planet)
 * \param relativeStateOfBodyExertingTide State of satellite w.r.t. host planet.
 * \param currentTidalAccelerationMultiplier Scalar multiplier that is common to all vector terms (term outside of brackets in
 * Eq. (3) of Lainey et al. (2007)).
 * \param timeLag Time lag of tidal bulge
 * \return Partial derivative of direct tidal acceleration due to tide on planet w.r.t. velocity of satellite
 */
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnPlanetWrtVelocity(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const double currentTidalAccelerationMultiplier,
        const double timeLag );

//! Function to compute partial derivative of direct tidal acceleration due to tide on satellite w.r.t. position of satellite
/*!
 * Function to compute partial derivative of direct tidal acceleration due to tide on satellite w.r.t. position of satellite (equal
 * to -1 * partial w.r.t. position of planet)
 * \param relativeStateOfBodyExertingTide State of satellite w.r.t. host planet.
 * \param currentTidalAccelerationMultiplier Scalar multiplier that is common to all vector terms (term outside of brackets in
 * Eq. (3) of Lainey et al. (2007)).
 * \param timeLag Time lag of tidal bulge
 * \param includeDirectRadialComponent True if term independent of time lag is to be included, false otherwise
 * \return Partial derivative of direct tidal acceleration due to tide on satellite w.r.t. position of satellite
 */
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnSatelliteWrtPosition(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double currentTidalAccelerationMultiplier, const double timeLag, const bool includeDirectRadialComponent );

//! Function to compute partial derivative of direct tidal acceleration due to tide on satellite w.r.t. velocity of satellite
/*!
 * Function to compute partial derivative of direct tidal acceleration due to tide on satellite w.r.t. velocity of satellite (equal
 * to -1 * partial w.r.t. velocity of planet)
 * \param relativeStateOfBodyExertingTide State of satellite w.r.t. host planet.
 * \param currentTidalAccelerationMultiplier Scalar multiplier that is common to all vector terms (term outside of brackets in
 * Eq. (3) of Lainey et al. (2007)).
 * \param timeLag Time lag of tidal bulge
 * \return Partial derivative of direct tidal acceleration due to tide on satellite w.r.t. velocity of satellite
 */
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnSatelliteWrtVelocity(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const double currentTidalAccelerationMultiplier,
        const double timeLag );

//! Class to calculate the partials of the direct tidal acceleration model (e.g. Lainey et al. 2007 ) w.r.t. parameters and states
class DirectTidalDissipationAccelerationPartial: public AccelerationPartial
{
public:

    //! Constructor.
    /*!
     *  Constructor.
     *  \param tidalAcceleration Acceleration model of which partial derivatives are to be computed
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    DirectTidalDissipationAccelerationPartial(
            const boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > tidalAcceleration,
            const std::string acceleratedBody,
            const std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::direct_tidal_dissipation_acceleration ),
        tidalAcceleration_( tidalAcceleration )
    { }

    //! Destructor
    ~DirectTidalDissipationAccelerationPartial( ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }


    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
    }

    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  No dependency is implemented, but a warning is provided if partial w.r.t. mass of body exerting acceleration
     *  (and undergoing acceleration if mutual attraction is used) is requested.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedNonTranslationalState(
                const std::pair< std::string, std::string >& stateReferencePoint,
                const propagators::IntegratedStateType integratedStateType )
    {
        if( ( stateReferencePoint.first == acceleratingBody_ ) ||
              ( stateReferencePoint.first == acceleratedBody_ ) )
        {
            throw std::runtime_error( "Warning, dependency of direct dissipation acceleration on non-translational states not yet implemented" );
        }
        return 0;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        boost::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    //! Function for updating partial w.r.t. the bodies' positions and velocities
    /*!
     *  Function for updating partial w.r.t. the bodies' positions and velocities
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN );

protected:

    //! Function to create a function returning the current partial w.r.t. a gravitational parameter.
    /*!
     * Function to create a function returning the current partial w.r.t. a gravitational parameter.
     * \param parameterId Identifier of parameter for which the partial is to be created.
     * \return Pair with partial function and paramater partial size.
     */
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getGravitationalParameterPartialFunction(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterId );

    //! Function to compute derivative w.r.t. gravitational parameter of planet
    void wrtGravitationalParameterOfPlanet( Eigen::MatrixXd& gravitationalParameterPartial );

    //! Function to compute derivative w.r.t. gravitational parameter of satellite
    void wrtGravitationalParameterOfSatellite( Eigen::MatrixXd& gravitationalParameterPartial );

    //! Function to compute derivative w.r.t. tidal time lag parameter.
    void wrtTidalTimeLag( Eigen::MatrixXd& gravitationalParameterPartial );

    //! Acceleration model of which partial derivatives are to be computed
    boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > tidalAcceleration_;

    //! Current state of satellite w.r.t. planet
    Eigen::Vector6d currentRelativeBodyState_;

    //! Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
    /*!
     *  Current partial of central gravity acceleration w.r.t. position of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. position of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentPartialWrtPosition_;


    //! Current partial of central gravity acceleration w.r.t. velocity of body undergoing acceleration
    /*!
     *  Current partial of central gravity acceleration w.r.t. velocity of body undergoing acceleration
     * ( = -partial of central gravity acceleration w.r.t. velocity of body exerting acceleration),
     *  calculated and set by update( ) function.
     */
    Eigen::Matrix3d currentPartialWrtVelocity_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_DIRECTTIDALDISSIPATIONACCELERATIONPARTIAL_H
