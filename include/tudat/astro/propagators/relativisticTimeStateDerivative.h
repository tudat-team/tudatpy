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

/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTICTIMESTATEDERIVATIVE_H
#define TUDAT_RELATIVISTICTIMESTATEDERIVATIVE_H

#include <map>
#include <vector>

#include <functional>

#include <string>
#include <stdexcept>

#include "tudat/math/basic/linearAlgebra.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"
#include "tudat/astro/ephemerides/customEphemeris.h"
#include "tudat/simulation/propagation_setup/createEnvironmentUpdater.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"
#include "tudat/astro/relativity/relativisticPotentials.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/astro/relativity/relativisticTimeConversion.h"
#include "tudat/astro/relativity/relativisticEquationsOfMotion.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{


//! Base class for relativistic time difference calculations.
/*!
 *  Base class for relativistic time-derivative propagation models.
 *  Derived classes implement specific models for:
 *  - barycentric <-> bodycentric coordinate-time conversion (Soffel et al., 2003, Eq. (58));
 *  - bodycentric <-> topocentric/proper-time conversion (Turyshev et al., 2013, Eq. (22));
 *  - direct proper-time rate from a metric model.
 *
 *  The propagated state is the relative time quantity (typically \f$\tau-t\f$ or
 *  barycentric/bodycentric coordinate-time difference) and the state derivative equals the
 *  corresponding model equation evaluated at the current environment state.
 */

template< typename StateScalarType = double, typename TimeType = double >
class RelativisticTimeStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;

public:

    //! Constructor.
    /*!
     *  \param referencePointId Pair identifying body and optional reference point.
     *  \param relativisticStateDerivativeType Relativistic time-derivative model type.
     */
    RelativisticTimeStateDerivative(
            const std::pair< std::string, std::string >& referencePointId,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >( propagators::proper_time ),
        relativisticStateDerivativeType_( relativisticStateDerivativeType ), referencePointId_( referencePointId )
    { }

    virtual ~RelativisticTimeStateDerivative( ){ }


    //! Convert internal propagated state to global representation.
    /*!
     *  \param internalSolution Internal propagated state.
     *  \param time Current time.
     *  \param currentCartesianLocalSolution Output global representation block.
     */
    void convertCurrentStateToGlobalRepresentation(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution,
        const TimeType& time,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution ) override
    {
        currentCartesianLocalSolution = internalSolution;
    }


    //! Convert output solution back to internal representation.
    /*!
     *  \param outputSolution Output solution history.
     *  \param time Current time.
     *  \return Converted solution history.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time ) override
    {
        return outputSolution;
    }

    //! Convert internal solution to output representation.
    /*!
     *  \param internalSolution Internal propagated state.
     *  \param time Current time.
     *  \param currentCartesianLocalSolution Output solution block.
     */
    void convertToOutputSolution(
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
        const TimeType& time,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution ) override
    {
        currentCartesianLocalSolution = internalSolution;
    }


    //! Get propagated state size.
    /*!
     *  \return State size.
     */
    virtual int getStateSize( )
    {
        return 1;
    }

    //! Get conventional state size.
    /*!
     *  \return Conventional state size.
     */
    int getConventionalStateSize( ) override
    { 
        return 1;
    }

    //! Clear internal model state (no-op).
    void clearStateDerivativeModel( ) override { }


    //! Get central body name.
    /*!
     *  \return Central body name.
     */
    std::string getCentralBody( )
    {
        return referencePointId_.first;
    }

    //! Get reference point identifier.
    /*!
     *  \return Pair of body and optional reference-point name.
     */
    std::pair< std::string, std::string > getReferencePoint( )
    {
        return referencePointId_;
    }

    //! Get relativistic derivative type.
    /*!
     *  \return Relativistic derivative model type.
     */
    propagators::RelativisticTimeStateDerivativeType getRelativisticStateDerivativeType( )
    {
        return relativisticStateDerivativeType_;
    }

protected:

    propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType_;

    std::pair< std::string, std::string > referencePointId_;
};

template< typename StateScalarType = double, typename TimeType = double >
class DirectProperTimeRateStateDerivative: public RelativisticTimeStateDerivative< StateScalarType, TimeType >
{
public:
    //! Constructor.
    /*!
     *  \param spaceTimeMetric Space-time metric used in proper-time-rate evaluation.
     *  \param referencePointId Reference point of propagated proper-time state.
     *  \param bodies System of bodies providing point-state functions.
     */
    DirectProperTimeRateStateDerivative(
            const std::shared_ptr< relativity::Metric > spaceTimeMetric,
            const std::pair< std::string, std::string > referencePointId,
            const simulation_setup::SystemOfBodies bodies ):
            RelativisticTimeStateDerivative< StateScalarType, TimeType >( referencePointId, propagators::direct_from_metric ),
        spaceTimeMetric_( spaceTimeMetric )
    {
        if( bodies.getBody( referencePointId.first ) == 0 )
        {
            throw std::runtime_error(
                        "Error when making direct proper time state derivative, could not find body " +
                        referencePointId.first );
        }
        else
        {
            if( referencePointId.second == "" )
            {
                referencePointStateFunction_ = std::bind( &simulation_setup::Body::getState, bodies.getBody( referencePointId.first ) );

            }
            else
            {
                referencePointStateFunction_ = simulation_setup::getLinkEndCompleteCurrentStateFunctionDuringPropagation< double >(
                            referencePointId, bodies );
            }
        }
    }

    //! Evaluate proper-time state derivative.
    /*!
     *  \param time Current evaluation time.
     *  \param stateOfSystemToBeIntegrated Current propagated time state.
     *  \param stateDerivative Output state derivative.
     */
    void calculateSystemStateDerivative(
        const TimeType time,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative( 0, 0 ) = currentProperTimeDerivative_;
    }


    //! Update environment-dependent quantities and current derivative value.
    /*!
     *  \param currentTime Current evaluation time.
     */
    void updateStateDerivativeModel( const TimeType currentTime )
    {
        if( spaceTimeMetric_ == nullptr )
        {
            throw std::runtime_error( "DirectProperTimeRateStateDerivative: spaceTimeMetric_ is nullptr" );
        }
        if( !referencePointStateFunction_ )
        {
            throw std::runtime_error( "DirectProperTimeRateStateDerivative: referencePointStateFunction_ is empty" );
        }

        const Eigen::Vector6d referenceState = referencePointStateFunction_( );
        spaceTimeMetric_->update( referenceState, currentTime, 1, 0  );
        currentProperTimeDerivative_ = evaluateProperTimeEquation( spaceTimeMetric_, referenceState, 2 );
    }

    //! Get metric object used by this model.
    /*!
     *  \return Shared pointer to metric object.
     */
    std::shared_ptr< relativity::Metric > getSpaceTimeMetric( )
    {
        return spaceTimeMetric_;
    }

    //! Get most recently computed proper-time derivative.
    /*!
     *  \return Current proper-time derivative.
     */
    double getCurrentProperTimeDerivative( )
    {
        return currentProperTimeDerivative_;
    }

    //! Function returning the reference-point Cartesian state used at each evaluation.
    /*!
     *  Wraps the function bound at construction (either the body's own state or, for a ground
     *  station reference, the topocentric state). Used by GR/SR-split dependent variables to
     *  read out the BCRS velocity that drives the kinematic time-dilation term.
     */
    std::function< Eigen::Vector6d( ) > getReferencePointStateFunction( ) const
    {
        return referencePointStateFunction_;
    }


protected:
    std::shared_ptr< relativity::Metric > spaceTimeMetric_;

    std::pair< std::string, std::string > referencePointId_;

    std::function< Eigen::Vector6d( ) > referencePointStateFunction_;

    double currentProperTimeDerivative_;
};

template< typename StateScalarType = double, typename TimeType = double >
class PostNewtonianRelativisticTimeStateDerivative: public RelativisticTimeStateDerivative< StateScalarType, TimeType >
{
public:

    //! Base class constructor, sets common member variables for derived classes.
    /*!
     *  Base class constructor, sets common member variables for derived classes.
     *  \param bodies Map of body objects.
     *  \param centralBody Name of central body
     *  \param externalBodies List of external body names
     */
    PostNewtonianRelativisticTimeStateDerivative(
            const simulation_setup::SystemOfBodies& bodies,
            const std::pair< std::string, std::string >& referencePoint,
            const std::vector< std::string >& externalBodies,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType,
            const std::function< TimeType( const TimeType& ) > timeVariableConversionFunction =
                []( const TimeType& inputTime ){ return inputTime; },
            const double distanceScalingFactor = 1.0 ):
        RelativisticTimeStateDerivative< StateScalarType, TimeType >( referencePoint, relativisticStateDerivativeType ),
        timeVariableConversionFunction_( timeVariableConversionFunction ),distanceScalingFactor_( distanceScalingFactor )
    {
        // Check if central body exists.
        if( referencePoint.first == "SSB" )
        {
            centralBodyStateFunction_ = []( const TimeType& ){ return Eigen::Vector6d::Zero( ); };
        }
        else if( !bodies.doesBodyExist( referencePoint.first ) )
        {
            throw std::runtime_error( "Error, could not find centeal body in body map when making 1st order relativistic time converter" );
        }
        else
        {
            // Set central body state function.
            const std::shared_ptr< simulation_setup::Body > centralBody = bodies.getBody( referencePoint.first );
            centralBodyStateFunction_ = [=]( const TimeType& evaluationTime )
            {
                return centralBody->getStateInBaseFrameFromEphemeris< double, TimeType >( evaluationTime );
            };
        }

        // Iterate over all external bodies and set their state and gravitational parameter functions.
        for( unsigned int i = 0; i < externalBodies.size( ); i++ )
        {
            // Check if current body exists.
            if( !bodies.doesBodyExist( externalBodies.at( i ) ) )
            {
                throw std::runtime_error(
                            "Error, could not find body " + externalBodies.at( i ) +
                            " in body map when making 1st order relativistic time converter" );
            }
            // Check if current body is a CelestialBody
            else if( std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( externalBodies.at( i ) ) ) == nullptr )
            {
                throw std::runtime_error(
                            "Error, could not find celestial body " + externalBodies.at( i ) +
                            " in body map when making 1st order relativistic time converter" );
            }
            // Check if current body has a gravity field.
            else if( std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( externalBodies.at( i ) ) )->getGravityFieldModel( ) == nullptr )
            {
                throw std::runtime_error(
                            "Error, celestial body " + externalBodies.at( i ) +
                            " has no gravity field model when making 1st order relativistic time converter" );
            }
            else
            {
                // Set state and gravitational parameter functions of current body.
                std::shared_ptr< simulation_setup::Body > body =
                        std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( externalBodies.at( i ) ) );
                this->externalBodyStateFunctions_.push_back(
                    [=]( const TimeType& evaluationTime )
                    {
                        return body->getStateInBaseFrameFromEphemeris< double, TimeType >( evaluationTime );
                    } );
                this->externalBodyGravitationalParameterFunctions_.push_back( std::bind( &gravitation::GravityFieldModel::getGravitationalParameter,
                                                                                    body->getGravityFieldModel( ) ) );
            }
        }
    }

    //!  Destructor
    /*!
       *   Destructor
       */
    virtual ~PostNewtonianRelativisticTimeStateDerivative( ){ }


    //! Get optional time-conversion function used prior to environment updates.
    /*!
     *  \return Time-conversion function.
     */
    std::function< TimeType( const TimeType& ) > getTimeVariableConversionFunction( )
    {
        return timeVariableConversionFunction_;
    }

protected:


    //! State function of central body.
    std::function< Eigen::Vector6d( const TimeType& ) > centralBodyStateFunction_;

    //! List of functions returning gravitational parameters of bodies influencing time conversion.
    std::vector< std::function< double( ) > > externalBodyGravitationalParameterFunctions_;

    //! List of functions returning Cartesian states as function of baseFrameTime (see updateStateDerivativeModel) of bodies influencing time conversion.
    std::vector< std::function< Eigen::Vector6d( const TimeType& ) > > externalBodyStateFunctions_;

    //! Bodies for which spherical-harmonic corrections are applied (index matches external body list).
    std::map< int, std::shared_ptr< simulation_setup::Body > > sphericalHarmonicBodies_;

    std::shared_ptr< propagators::EnvironmentUpdater< double, double > > environmentUpdater_;

    std::function< TimeType( const TimeType& ) > timeVariableConversionFunction_;

    double distanceScalingFactor_;
};

//! First-order barycentric->bodycentric time-derivative model (Soffel et al., Eq. (58), \f$c^{-2}\f$ part).
template< typename StateScalarType = double, typename TimeType = double >
class FirstOrderBarycentricToBodyCentricTimeStateDerivative: public PostNewtonianRelativisticTimeStateDerivative< StateScalarType, TimeType >
{
public:

    //! Class contructor, sets required properties of bodies.
    /*!
     *  Class contructor, sets required properties of bodies, such as state functions and gravitational parameters.
     *  \param bodies Map of all bodies, from which required functions will be retrieved and set as member variables in this class
     *  \param centralBody Name of central body used in conversion
     *  \param externalBodies List of bodies causing gravitational variation of time conversion, from which 1st order contributions will
     *  be calculated and combined.
     *  \param sphericalHarmonicGravityExpansions Optional spherical-harmonic degree/order settings per external body.
     *  \param timeVariableConversionFunction Optional time conversion applied before environment updates.
     *  \param distanceScalingFactor Optional scale factor applied to positions.
     *  \param relativisticStateDerivativeType Internal derivative type identifier.
     */
    FirstOrderBarycentricToBodyCentricTimeStateDerivative(
            const simulation_setup::SystemOfBodies& bodies,
            const std::string& centralBody,
            const std::vector< std::string >& externalBodies,
            const std::map< std::string, std::pair< int, int > > sphericalHarmonicGravityExpansions = ( std::map< std::string, std::pair< int, int > >( ) ),
            const std::function< TimeType( const TimeType& ) > timeVariableConversionFunction =
                []( const TimeType& inputTime ){ return inputTime; },
            const double distanceScalingFactor = 1.0,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType =
            propagators::first_order_barycentric_to_bodycentric):
                PostNewtonianRelativisticTimeStateDerivative< StateScalarType, TimeType >(
                    bodies, std::make_pair( centralBody, "" ), externalBodies, relativisticStateDerivativeType,
                    timeVariableConversionFunction, distanceScalingFactor ), externalBodies_( externalBodies )
        {
            // Set sizes of vectors which are to hold the current states of the system.
            this->currentExternalBodyDistances_.resize( this->externalBodyStateFunctions_.size( ) );
            this->currentExternalBodyGravitationalParameters_.resize( this->externalBodyGravitationalParameterFunctions_.size( ) );
            this->currentExternalBodyStates_.resize( this->externalBodyStateFunctions_.size( ) );

            for( unsigned int i = 0; i < externalBodies.size( ); i++ )
            {
                if( sphericalHarmonicGravityExpansions.count( externalBodies.at( i ) ) > 0 )
                {
                    std::shared_ptr< simulation_setup::Body > currentCelestialBody = std::dynamic_pointer_cast< simulation_setup::Body >(
                                bodies.getBody( externalBodies.at( i ) ) );
                    if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( currentCelestialBody->getGravityFieldModel( ) ) == nullptr )
                    {
                        throw std::runtime_error(
                                    "Error when making sh correction for relativistic time conversion, body " +
                                    externalBodies.at( i ) + " has no sh gravity field." );
                    }
                    else
                    {

                        auto sphericalHarmonicsCache = std::make_shared< basic_mathematics::SphericalHarmonicsCache >(
                            sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).first,
                            sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).second );
                        std::shared_ptr< gravitation::SphericalHarmonicsGravityField > shGravityField =
                            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( currentCelestialBody->getGravityFieldModel( ) );
                        this->sphericalHarmonicBodies_[ i ] = currentCelestialBody;

                        higherOrderGravityFieldPotentialFunctions_[ i ] = [=]( const Eigen::Vector3d& position ) {
                            return shGravityField->getGravitationalPotentialFromInertialPosition(
                                    position,
                                    currentCelestialBody->getCurrentRotationToLocalFrame( ),
                                    static_cast< double >( sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).first ),
                                    static_cast< double >( sphericalHarmonicGravityExpansions.at( externalBodies.at( i ) ).second ),
                                    sphericalHarmonicsCache );
                        };
                    }
                }
            }
        }

    //! Destructor.
    virtual ~FirstOrderBarycentricToBodyCentricTimeStateDerivative( ){ }

    //! Function for evaluating the integrand which is required for the conversion.
    /*!
     *  Evaluates the first-order barycentric/body-centric conversion contribution
     *  from the Soffel et al. (2003) framework (Eq. (58), first-order terms).
     *  The sign convention is set by the selected conversion direction.
     *  \param currentGlobalTime Current time, in frame from which conversion is done (i.e. independent variable of equations of motion of bodies)
     *  \param integratedValue Value of integral up to current time.
     *  \return Value of integrand at current time.
     */
    void calculateSystemStateDerivative(
        const TimeType currentGlobalTime,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative( 0, 0 ) = relativity::calculateFirstOrderTcbToTcgIntegrand(
            this->currentVelocity_, this->currentExternalPotential_ );
    }

    //! Get last-evaluated central-body BCRS speed |v_E|.
    /*!
     *  Returns the central-body barycentric speed used in the most recent integrand evaluation.
     *  Set by updateStateDerivativeModel(); valid after at least one update. Inherited unchanged
     *  by the second-order subclass.
     */
    double getCurrentCentralBodySpeed( ) const
    {
        return this->currentVelocity_;
    }

    //! Get last-evaluated external scalar potential U_ext at the central body.
    /*!
     *  Returns Sigma_i (mu_i / r_iE) over the configured external bodies (with optional
     *  higher-order spherical-harmonic gravity contributions if requested at construction time).
     *  Used together with the central-body speed to split the first-order proper-time-rate
     *  integrand into kinematic (-v^2/(2 c^2)) and potential (-U/c^2) contributions for
     *  diagnostics. Inherited unchanged by the second-order subclass.
     */
    double getCurrentExternalScalarPotential( ) const
    {
        return this->currentExternalPotential_;
    }

    //! Function to update all environment variables to current time.
    /*!
     *  Function to update environment variables (central body velocity and external potential) to current time.
     *  \param baseFrameTime Current time, in frame from which conversion is done (i.e. independent variable of equations of motion of bodies)
     */
    void updateStateDerivativeModel( const TimeType baseFrameTime )
    {
        // Update variables relevant for both this class and its derived class (2nd order conversion)
        updateBaseVariables( baseFrameTime );
        const TimeType convertedTime = this->timeVariableConversionFunction_( baseFrameTime );
        // Calculate state, and distance to central body, of each external body.
        for( unsigned int i = 0; i < this->externalBodyStateFunctions_.size( ); i++ )
        {
            this->currentExternalBodyStates_[ i ] = this->externalBodyStateFunctions_[ i ]( convertedTime );
            this->currentExternalBodyStates_[ i ].segment( 0, 3 ) = this->currentExternalBodyStates_[ i ].segment( 0, 3 ) * this->distanceScalingFactor_;

            this->currentExternalBodyDistances_[ i ] = ( this->currentCentralBodyState_.segment( 0, 3 ) - this->currentExternalBodyStates_[ i ].segment( 0, 3 ) ).norm( );
        }
        // Calculate 1st order external potential.
        this->currentExternalPotential_ = 0.0;
        for( unsigned int i = 0; i < this->externalBodyStateFunctions_.size( ); i++ )
        {
            const Eigen::Vector3d relativePosition =
                    this->currentCentralBodyState_.segment( 0, 3 ) - this->currentExternalBodyStates_[ i ].segment( 0, 3 );

            if( higherOrderGravityFieldPotentialFunctions_.count( i ) > 0 )
            {
                this->currentExternalPotential_ += higherOrderGravityFieldPotentialFunctions_.at( i )( relativePosition );
            }
            else
            {
                this->currentExternalPotential_ += this->currentExternalBodyGravitationalParameters_[ i ] / relativePosition.norm( );
            }
        }
    }


    double getCurrentExternalBodyDistance( const int index )
    {
        return currentExternalBodyDistances_.at( index );
    }

    double getCurrentExternalBodyGravitationalParameter( const int index )
    {
        return currentExternalBodyGravitationalParameters_.at( index );
    }

    Eigen::Vector6d getCurrentExternalBodyStates( const int index )
    {
        return currentExternalBodyStates_.at( index );
    }

    std::vector< std::string > getExternalBodies( )
    {
        return externalBodies_;
    }

    Eigen::Vector6d getCurrentCentralBodyState( )
    {
        return currentCentralBodyState_;
    }

    //! Function to update all environment variables common to this class, and its (2nd order) derived class to current time.
    /*!
     *  Function to update all environment variables common to this class, and its (2nd order) derived class to current time, namely
     *  body states and gravitational parameters.
     *  \param baseFrameTime Current time, in frame from which conversion is done (i.e. independent variable of equations of motion of bodies)
     */
    void updateBaseVariables( const TimeType& baseFrameTime )
    {
        const TimeType convertedTime = this->timeVariableConversionFunction_( baseFrameTime );

        std::map< propagators::IntegratedStateType, Eigen::VectorXd > currentState;
        currentState[ propagators::proper_time ] = ( Eigen::VectorXd( 1 ) << 0.0 ).finished( );

        // Ensure body-fixed frames used in spherical harmonics are evaluated at the current epoch.
        for( const auto& bodyEntry : this->sphericalHarmonicBodies_ )
        {
            bodyEntry.second->setCurrentRotationalStateToLocalFrameFromEphemeris( convertedTime );
        }
        //environmentUpdater_->updateEnvironment( convertedTime, currentState );

        // Calculate state of central body.
        this->currentCentralBodyState_ = this->centralBodyStateFunction_( convertedTime );
        this->currentCentralBodyState_.segment( 0, 3 ) = this->currentCentralBodyState_.segment( 0, 3 ) * this->distanceScalingFactor_;

        this->currentVelocity_ = this->currentCentralBodyState_.segment( 3, 3 ).norm( );
        // Retrieve gravitational parameters of all bodies.
        for( unsigned int i = 0; i < this->externalBodyStateFunctions_.size( ); i++ )
        {
            this->currentExternalBodyGravitationalParameters_[ i ] = this->externalBodyGravitationalParameterFunctions_[ i ]( ) * this->distanceScalingFactor_;
        }
    }

protected:

    //! Current speed of central body
    double currentVelocity_;

    //! Current state of central body
    Eigen::Vector6d currentCentralBodyState_;

    //! Current distances from external bodies to central body.
    std::vector< double > currentExternalBodyDistances_;

    //! Current gravitational parameters of external bodies.
    std::vector< double > currentExternalBodyGravitationalParameters_;

    //! Current states of external bodies.
    std::vector< Eigen::Vector6d > currentExternalBodyStates_;

    std::vector< std::string > externalBodies_;

    //! Current scalar potential due to external bodies.
    double currentExternalPotential_;

    std::map< int, std::function< double( Eigen::Vector3d ) > > higherOrderGravityFieldPotentialFunctions_;// only dergee > 0 term contributions, added for 1st order variation.

    std::map< int, std::function< double( Eigen::Vector3d ) > >::iterator shPotentialIterator_;

    std::map< int, double > fullShPotentials_;
};

//! Second-order barycentric->bodycentric time-derivative model (Soffel et al., Eq. (58), up to \f$c^{-4}\f$).
template< typename StateScalarType = double, typename TimeType = double >
class SecondOrderBarycentricToBodyCentricTimeStateDerivative: public FirstOrderBarycentricToBodyCentricTimeStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor.
    /*!
     *  \param bodies System of bodies.
     *  \param centralBody Name of central body used in conversion.
     *  \param externalBodies External bodies used in conversion.
     *  \param sphericalHarmonicGravityExpansions Optional spherical-harmonic degree/order settings per external body.
     *  \param timeVariableConversionFunction Optional time conversion applied before environment updates.
     *  \param distanceScalingFactor Optional scale factor applied to positions.
     *  \param angularMomentumBodies Bodies for which angular-momentum terms are included.
     */
    SecondOrderBarycentricToBodyCentricTimeStateDerivative(
            const simulation_setup::SystemOfBodies& bodies,
            const std::string centralBody,
            const std::vector< std::string >& externalBodies,
            const std::map< std::string, std::pair< int, int > > sphericalHarmonicGravityExpansions = ( std::map< std::string, std::pair< int, int > >( ) ),
            const std::function< TimeType( const TimeType& ) > timeVariableConversionFunction =
                []( const TimeType& inputTime ){ return inputTime; },
            const double distanceScalingFactor = 1.0,
            const std::vector< std::string >& angularMomentumBodies = std::vector< std::string >( ) ): 
                FirstOrderBarycentricToBodyCentricTimeStateDerivative< StateScalarType, TimeType >(
                    bodies, centralBody, externalBodies, sphericalHarmonicGravityExpansions, timeVariableConversionFunction,
                    distanceScalingFactor, propagators::second_order_barycentric_to_bodycentric )
            {
                for( unsigned int i = 0; i < angularMomentumBodies.size( ); i++ )
                {
                    if( angularMomentumBodies.at( i ) == "Sun" )
                    {
                        angularMomentumFunctions_[ i ] = []() { return ( Eigen::Vector3d() << 0.0, 0.0, 1.92E41 ).finished(); };
                    }
                    else if( angularMomentumBodies.at( i ) == "Jupiter" )
                    {
                        angularMomentumFunctions_[ i ] = []() { return ( Eigen::Vector3d() << 0.0, 0.0, 6.9E38 ).finished(); };
                    }
                    else if( angularMomentumBodies.at( i ) == "Saturn" )
                    {
                        angularMomentumFunctions_[ i ] = []() { return ( Eigen::Vector3d() << 0.0, 0.0, 1.42E38 ).finished(); };
                    }
                    //std::shared_ptr< simulation_setup::Body > celestialBody =
                    //        std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( externalBodies.at( i ) ) );
                    //this->externalBodyStateFunctions_.push_back( std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris, bodies.getBody( externalBodies.at( i ) ), std::placeholders::_1 ) );
                }
                currentBodySpeeds_.resize( this->externalBodyStateFunctions_.size( ) );
            }

    ~SecondOrderBarycentricToBodyCentricTimeStateDerivative( ){ }

    //! Evaluate first+second-order Soffel Eq. (58) integrand.
    /*!
     *  \param currentGlobalTime Current global coordinate time.
     *  \param stateOfSystemToBeIntegrated Current propagated time state.
     *  \param stateDerivative Output state derivative.
     */
    void calculateSystemStateDerivative(
        const TimeType currentGlobalTime,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative( 0, 0 ) =
            relativity::calculateFirstOrderTcbToTcgIntegrand( this->currentVelocity_, this->currentExternalPotential_ ) +
            relativity::calculateSecondOrderTcbToTcgIntegrand(
                this->currentVelocity_,
                this->currentExternalPotential_,
                this->currentCentralBodyState_.segment( 3, 3 ),
                currentExternalVectorPotential_,
                currentSecondOrderExternalPotentialCorrection_ );
    }


    //! Update all environment-dependent quantities used in second-order integrand.
    /*!
     *  \param baseFrameTime Current global coordinate time.
     */
    void updateStateDerivativeModel( const TimeType baseFrameTime )
    {
        this->updateBaseVariables( baseFrameTime );
        const TimeType convertedTime = this->timeVariableConversionFunction_( baseFrameTime );

        for( std::map< int, std::function< Eigen::Vector3d( ) > >::iterator momentaIterator = angularMomentumFunctions_.begin( );
            momentaIterator != angularMomentumFunctions_.end( ); momentaIterator++ )
        {
            currentAngularMomenta_[ momentaIterator->first ] = momentaIterator->second( );
        }

        for( unsigned int i = 0; i < this->externalBodyStateFunctions_.size( ); i++ )
        {
            this->currentExternalBodyStates_[ i ] = this->externalBodyStateFunctions_[ i ]( convertedTime );
            this->currentExternalBodyStates_[ i ].segment( 0, 3 ) = this->currentExternalBodyStates_[ i ].segment( 0, 3 ) * this->distanceScalingFactor_;

            this->currentExternalBodyDistances_[ i ] = ( this->currentCentralBodyState_.segment( 0, 3 ) - this->currentExternalBodyStates_[ i ].segment( 0, 3 ) ).norm( );

            currentBodySpeeds_[ i ] = this->currentExternalBodyStates_[ i ].segment( 3, 3 ).norm( );
        }

        this->currentExternalPotential_ = 0.0;
        for( unsigned int i = 0; i < this->externalBodyStateFunctions_.size( ); i++ )
        {
            const Eigen::Vector3d relativePosition =
                    this->currentCentralBodyState_.segment( 0, 3 ) - this->currentExternalBodyStates_[ i ].segment( 0, 3 );

            if( this->higherOrderGravityFieldPotentialFunctions_.count( i ) > 0 )
            {
                this->currentExternalPotential_ += this->higherOrderGravityFieldPotentialFunctions_.at( i )( relativePosition );
            }
            else
            {
                this->currentExternalPotential_ += this->currentExternalBodyGravitationalParameters_[ i ] / relativePosition.norm( );
            }
        }

        currentSecondOrderExternalPotentialCorrection_ = relativity::calculateSecondOrderExternalScalarPotentialCorrection(
                    this->currentExternalBodyGravitationalParameters_, this->currentExternalBodyDistances_, currentBodySpeeds_, this->currentExternalBodyStates_,
                    this->currentCentralBodyState_.segment( 3, 3 ) );
        currentExternalVectorPotential_ = relativity::calculateExternalVectorPotential( this->currentExternalBodyGravitationalParameters_, this->currentExternalBodyDistances_,
                                                                            this->currentExternalBodyStates_, currentAngularMomenta_ );
    }

private:

    std::map< int, std::function< Eigen::Vector3d( ) > > angularMomentumFunctions_;

    std::vector< double > currentBodySpeeds_;

    std::map< int, Eigen::Vector3d > currentAngularMomenta_;

    Eigen::Vector3d currentExternalVectorPotential_;

    double currentSecondOrderExternalPotentialCorrection_; //Delta term in Sofel et al.

};

//! First-order bodycentric->topocentric/proper-time model (Turyshev et al., 2013, Eq. (22)).
/*!
 *  Implements, for a station-fixed point \f$\mathbf{y}\f$:
 *  \f[
 *      \frac{d\tau}{dt_C} = 1-\frac{1}{c^2}\left[
 *      \frac{1}{2}v_0^2 + U_E(\mathbf{y}) +
 *      \sum_{b\neq E}\frac{GM_b}{2r_{bE}^3}\left(3(\mathbf{n}_{bE}\cdot\mathbf{y})^2-\mathbf{y}^2\right) +
 *      \mathbf{a}_E\cdot\mathbf{y}
 *      \right] + \mathcal{O}(y^3,c^{-4}).
 *  \f]
 */
template< typename StateScalarType = double, typename TimeType = double >
class FirstOrderBodyCentricToTopoCentricTimeCalculator: public PostNewtonianRelativisticTimeStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor of bodycentric<->topocentric time calculator for a single ground station on a celestial body.
    /*!
     *  Constructor of bodycentric<->topocentric time calculator for a single ground station on a celestial body. The class is
     *  basically a wrapper for automatically calculating and setting the result of Eq. (22) of Turyshev et al. (2013). The conversion
     *  to/from a barycentric frame is done by a separate instance of the FirstOrderBarycentricToBodyCentricTimeStateDerivative or
     *  SecondOrderBarycentricToBodyCentricTimeStateDerivative class.
     *  Function to add ground stations to converter, for inclusion of topocentic<->bodycentric times. For all added ground stations,
     *  automatic conversion between any of the topocentric, bodycentric and barycentric frames can be achieved.
     *  \param bodies List of all bodies in simulation.
     *  \param centralBody Name of body on which ground station is located.
     *  \param groundStation Name of ground station on centralBody for which topocentric frame is to be defined.
     *  \param maximumSphericalHarmonicDegree Value denoting maximum sh gravity degree of central body to use in evaluation of
     *  gravitational potential at ground station.
     *  \param externalBodies List of bodies for which the tidal contribution to the potential is to be included in the body<->topo-
     *  centric conversion.
     *  \param useAccelerationTerm Boolean denoting whether to use acceleration term in Eq. (22) of Turyshev et al.(2013)
     *  \param useTimeDependentBodyFixedPosition Boolean denoting whether use time-dependency of ground station position, or use its nominal
     *  position at all times.
     */
    FirstOrderBodyCentricToTopoCentricTimeCalculator( const simulation_setup::SystemOfBodies& bodies,
                                                      const std::string& centralBody,
                                                      const std::vector< std::string >& externalBodies,
                                                      const std::string& groundStation,
                                                      const int maximumSphericalHarmonicDegree = 0,
                                                      const bool useAccelerationTerm = 0,
                                                      const bool useTimeDependentBodyFixedPosition = 0 ): 
            PostNewtonianRelativisticTimeStateDerivative< StateScalarType, TimeType >(
            bodies, std::make_pair( centralBody, groundStation ), externalBodies, propagators::first_order_bodycentric_to_topocentric )
    {
        // Check if central body has a rotation to body-fixed frame.
        if( bodies.getBody( centralBody )->getRotationalEphemeris( ) == nullptr )
        {
            throw std::runtime_error(
                        "Error, could not find rotation model of central body when making PCRS to TPRS time transformation." );
        }
        else
        {
            // Get rotation vector and rotation matrix functions.
            centralBodyRotationVector_ = std::bind( &ephemerides::RotationalEphemeris::getRotationalVelocityVectorInBaseFrame,
                                                    bodies.getBody( centralBody )->getRotationalEphemeris( ), std::placeholders::_1 );
            toInertialFrameTransformation_ = std::bind( &ephemerides::RotationalEphemeris::getRotationToBaseFrame,
                                                        bodies.getBody( centralBody )->getRotationalEphemeris( ), std::placeholders::_1 );
        }

        // Check if central body is a celestial body
        if( std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( centralBody ) ) == nullptr )
        {
            throw std::runtime_error(
                        "Error, central body " + centralBody +
                        " is not a celestial body when making PCRS to TPRS time transformation." );
        }
        else
        {
            std::shared_ptr< simulation_setup::Body > centralCelestialBody =
                    std::dynamic_pointer_cast< simulation_setup::Body >( bodies.getBody( centralBody ) );

            if( centralCelestialBody->getGroundStationMap( ).count( groundStation ) == 0 )
            {
                throw std::runtime_error(
                            "Error, station " + groundStation + " not found on " + centralBody +
                            " when making PCRS to TPRS time transformation." );
            }
            else
            {
                // Get ground station position function.
                std::shared_ptr< ground_stations::GroundStationState > groundStationState = centralCelestialBody->getGroundStation(
                            groundStation )->getNominalStationState( );
                if( !useTimeDependentBodyFixedPosition )
                {
                    pointPositionFunctionInPcrs_ = std::bind( &ground_stations::GroundStationState::getNominalCartesianPosition, groundStationState );
                }
                else
                {            
                    pointPositionFunctionInPcrs_ = [=]( double t ) { return groundStationState->getCartesianPositionInTime( t, centralBody );};
                }

                // Get central body potential function.
                if( centralCelestialBody->getGravityFieldModel( ) == nullptr )
                {
                    throw std::runtime_error(
                                "Error, could not find gravity field model of central body when making PCRS to TPRS time transformation." );
                }
                else
                {
                    // If spherical harmonic expansion is requested, check if required field exists.
                    if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( centralCelestialBody->getGravityFieldModel( ) ) == nullptr
                            && maximumSphericalHarmonicDegree > 0 )
                    {
                        throw std::runtime_error(
                                    "Error, requested spherical harmonic terms when making PCRS to TPRS time transformation, but no such model found." );
                    }
                    // If no spherical harmonic expansion is request, check if field is purely central, and set function.
                    else if( std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >( centralCelestialBody->getGravityFieldModel( ) ) ==
                            nullptr && maximumSphericalHarmonicDegree == 0 )
                    {
                        localCentralBodyPotentialFunction_ = std::bind( &gravitation::GravityFieldModel::getGravitationalPotential,
                                                                        centralCelestialBody->getGravityFieldModel( ), std::placeholders::_1 );
                    }
                    else
                    {
                        // Retrieve spherical harmonic potential function from gravtiy field object.
                        
                        auto sphericalHarmonicField =
                            std::dynamic_pointer_cast< gravitation::SphericalHarmonicsGravityField >(
                                centralCelestialBody->getGravityFieldModel() );

                        localCentralBodyPotentialFunction_ = [=]( const Eigen::Vector3d& position ) {
                            return sphericalHarmonicField->getGravitationalPotential(
                                position,
                                maximumSphericalHarmonicDegree,
                                maximumSphericalHarmonicDegree,
                                0,
                                0 );
                        };

                    }

                }
            }
        }

        // Add acceleration function, if requested.
        if( useAccelerationTerm )
        {
            throw std::runtime_error(
                        "Error, acceleration term not yet implemented for PCRS to TPRS time transformation." );
        }
        else
        {
            centralBodyAccelerationFunction_ = []( double ) { return Eigen::Vector3d::Zero( ); };
        }

        this->currentExternalBodyGravitationalParameters_.resize( this->externalBodyStateFunctions_.size( ) );
        currentExternalBodyRelativePositions_.resize( this->externalBodyStateFunctions_.size( ) );
        this->currentExternalBodyDistances_.resize( this->externalBodyStateFunctions_.size( ) );
    }

    //! Destructor.
    ~FirstOrderBodyCentricToTopoCentricTimeCalculator( ){ }

    //! Function for evaluating the integrand which is required for the conversion.
    /*!
     *  Evaluates the Turyshev et al. (2013), Eq. (22) contribution for the currently selected
     *  body-fixed point and environment state.
     *  \param currentGlobalTime Current time, in frame from which conversion is done (i.e. bodycentric coordinate time)
     *  \param integratedValue Value of integral up to current time.
     *  \return Value of integrand at current time.
     */
    void calculateSystemStateDerivative(
        const TimeType currentGlobalTime,
        const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
        Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative( 0, 0 ) = relativity::calculateFirstOrderPlanetocentricToTopocentricConversion(
            currentPointPositionInPcrs_,
            currentPointVelocityInPcrs_,
            currentLocalPotential_,
            currentBarycentricAccelerationOfCentralBody_,
            currentExternalBodyRelativePositions_,
            this->currentExternalBodyGravitationalParameters_ );
    }



    Eigen::Vector3d getGroundStationPositionInBodyCenteredInertialFrame( const double time )
    {
        return toInertialFrameTransformation_( time ) * pointPositionFunctionInPcrs_( time );
    }

    //! Get the latest cached reference-point speed in the body-centred inertial frame.
    /*!
     *  Returns \f$\lVert\mathbf{v}_0\rVert\f$, with \f$\mathbf{v}_0=\boldsymbol{\omega}\times\mathbf{y}\f$
     *  the inertial velocity of the topocentric reference point relative to the body centre as
     *  used in Turyshev et al. 2013 Eq. (22). Set by ``updateStateDerivativeModel``; valid after
     *  at least one update.
     */
    double getCurrentReferencePointSpeed( ) const
    {
        return currentPointVelocityInPcrs_.norm( );
    }

    //! Get the latest cached scalar potential at the topocentric reference point.
    /*!
     *  Returns \f$U_E(\mathbf{y}) + \sum_{i\neq E}\frac{GM_i}{2 r_{iE}^3}(3(\hat{\mathbf{n}}_{iE}\cdot
     *  \mathbf{y})^2-\mathbf{y}^2) + \mathbf{a}_E\cdot\mathbf{y}\f$, the local body potential plus
     *  the third-body tidal sum and (optional) acceleration term that together form the
     *  potential-like contribution to the Turyshev et al. 2013 Eq. (22) integrand.
     *  Reproduces the term that ``calculateFirstOrderPlanetocentricToTopocentricConversion``
     *  multiplies by \f$-1/c^2\f$.
     */
    double getCurrentLocalPotentialAndTidalContribution( ) const
    {
        double total = currentLocalPotential_;
        total += currentBarycentricAccelerationOfCentralBody_.dot( currentPointPositionInPcrs_ );
        for( unsigned int i = 0; i < currentExternalBodyRelativePositions_.size( ); ++i )
        {
            const double distance = currentExternalBodyRelativePositions_[ i ].norm( );
            if( distance == 0.0 )
            {
                continue;
            }
            const Eigen::Vector3d unitVector = currentExternalBodyRelativePositions_[ i ] / distance;
            const double dotProduct = currentPointPositionInPcrs_.dot( unitVector );
            const double pointPositionMagnitude = currentPointPositionInPcrs_.norm( );
            total += this->currentExternalBodyGravitationalParameters_[ i ] /
                    ( 2.0 * distance * distance * distance ) *
                    ( 3.0 * dotProduct * dotProduct - pointPositionMagnitude * pointPositionMagnitude );
        }
        return total;
    }


    //! Function to update all environment variables to current time.
    /*!
     *  Function to update environment variables to current time (central body potential, station state and central body state).
     *  \param baseFrameTime Current time, in frame from which conversion is done (i.e. bodycentric coordinate time).
     */
    void updateStateDerivativeModel( const TimeType baseFrameTime )
    {
        const TimeType convertedTime = this->timeVariableConversionFunction_( baseFrameTime );
        std::map< propagators::IntegratedStateType, Eigen::VectorXd > currentState;
        currentState[ propagators::proper_time ] = ( Eigen::VectorXd( 1 ) << 0.0 ).finished( );
        //environmentUpdater_->updateEnvironment( baseFrameTime, currentState );
        // Get central body state
        Eigen::Matrix< double, 6, 1 > currentCentralBodyState = this->centralBodyStateFunction_( convertedTime );
        Eigen::Vector3d currentCentralBodyPosition = currentCentralBodyState.segment( 0, 3 );
        currentCentralBodyBarycentricVelocity_ = currentCentralBodyState.segment( 3, 3 );

        // Get gravitational parameter and positions of all external bodies.
        for( unsigned int i = 0; i < this->externalBodyStateFunctions_.size( ); i++ )
        {
            this->currentExternalBodyGravitationalParameters_[ i ] = this->externalBodyGravitationalParameterFunctions_[ i ]( );
            currentExternalBodyRelativePositions_[ i ] =
                this->externalBodyStateFunctions_[ i ]( convertedTime ).segment( 0, 3 ) - currentCentralBodyPosition;
            this->currentExternalBodyDistances_[ i ] = currentExternalBodyRelativePositions_[ i ].norm( );
        }

        // Get current position of ground station in body-fixed frame.
        currentPointPositionInPcrs_ = pointPositionFunctionInPcrs_( baseFrameTime );

        // Evaluate gravitational potential of central body at current postion of ground station.
        currentLocalPotential_ = localCentralBodyPotentialFunction_( currentPointPositionInPcrs_ );

        // Convert ground station position to non-co-rotating frame.
        currentPointPositionInPcrs_ = toInertialFrameTransformation_( baseFrameTime ) * currentPointPositionInPcrs_;

        // Get current velocity of ground station in non-co-rotating frame.
        currentPointVelocityInPcrs_ = centralBodyRotationVector_( baseFrameTime ).cross( currentPointPositionInPcrs_ );

        // Calculate barycentric acceleration of central body.
        currentBarycentricAccelerationOfCentralBody_ = centralBodyAccelerationFunction_( baseFrameTime );

    }

private:


    std::function< Eigen::Vector3d( const double ) > pointPositionFunctionInPcrs_;

    std::function< Eigen::Vector3d( const double ) > centralBodyRotationVector_;

    std::function< double( Eigen::Vector3d ) > localCentralBodyPotentialFunction_;

    std::function< Eigen::Vector3d( const double ) > centralBodyAccelerationFunction_;

    std::function< Eigen::Quaterniond( const double ) > toInertialFrameTransformation_;


    std::vector< double > currentExternalBodyDistances_;

    std::vector< double > currentExternalBodyGravitationalParameters_;

    std::vector< Eigen::Vector3d > currentExternalBodyRelativePositions_;


    Eigen::Vector3d currentCentralBodyBarycentricVelocity_;

    Eigen::Vector3d currentBarycentricAccelerationOfCentralBody_;

    Eigen::Vector3d currentPointPositionInPcrs_;

    Eigen::Vector3d currentPointVelocityInPcrs_;

    double currentLocalPotential_;
};

}
#endif // TUDAT_RELATIVISTICTIMESTATEDERIVATIVE_H
