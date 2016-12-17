/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THRUSTSETTINGS_H
#define TUDAT_THRUSTSETTINGS_H

#include <Eigen/Geometry>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/Propulsion/thrustGuidance.h"
#include "Tudat/Astrodynamics/Propulsion/thrustMagnitudeWrapper.h"
#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! List of available types of thust direction guidance
enum ThrustDirectionGuidanceTypes
{
    colinear_with_state_segment_thrust_direction,
    thrust_direction_from_existing_body_orientation,
    custom_thrust_direction,
    custom_thrust_orientation

};

//! Class defining settings for the thrust direction
/*!
 *  Class for providing settings the thrust direction of a single thrust model. This class is a functional (base) class for
 *  settings of thrust direction that require no information in addition to their type.
 *  Classes defining settings for thrust direction requiring additional information must be derived from this class.
 */
class ThrustDirectionGuidanceSettings
{
public:

    //! Constructor
    /*!
    * Constructor
    * \param thrustDirectionType Type of thrust direction that is to be used.
    * \param relativeBody Body relative to which thrust guidance algorithm is defined (empty if N/A).
    */
    ThrustDirectionGuidanceSettings(
            const ThrustDirectionGuidanceTypes thrustDirectionType,
            const std::string relativeBody ):
        thrustDirectionType_( thrustDirectionType ), relativeBody_( relativeBody ){ }

    //! Destructor.
    virtual ~ThrustDirectionGuidanceSettings( ){ }

    //! Type of thrust direction that is to be used.
    ThrustDirectionGuidanceTypes thrustDirectionType_;

    //! Body relative to which thrust guidance algorithm is defined.
    std::string relativeBody_;
};

//! Thrust guidance settings for thrust that is colinear with position/velocity vector
class ThrustDirectionFromStateGuidanceSettings: public ThrustDirectionGuidanceSettings
{
public:

    //! Constructor
    /*!
    * Constructor
    * \param centralBody Body w.r.t. which the state of the body undergoing thust is computed. This state is then
    * used to directly set the thrust direction.
    * \param isColinearWithVelocity Boolean denoting whether thrust is colinear with velocity (if true) or position
    * (if false)
    * \param directionIsOppositeToVector Boolean denoting whether the thrust is in the direction of position/velocity
    * of x_{thrusting body}-x_{central body}, or opposite this direction
    */
    ThrustDirectionFromStateGuidanceSettings(
            const std::string& centralBody,
            const bool isColinearWithVelocity,
            const bool directionIsOppositeToVector ):
        ThrustDirectionGuidanceSettings( colinear_with_state_segment_thrust_direction, centralBody ),
        isColinearWithVelocity_( isColinearWithVelocity ),
        directionIsOppositeToVector_( directionIsOppositeToVector ){ }

    //! Destructor
    ~ThrustDirectionFromStateGuidanceSettings( ){ }

    //! Boolean denoting whether thrust is colinear with velocity (if true) or position (if false)
    bool isColinearWithVelocity_;

    //! Boolean denoting whether the thrust is in the direction of position/velocity of x_{thrusting body}-x_{central body}/
    bool directionIsOppositeToVector_;

};

//! Class for defining custom thrust direction (i.e. predefined thrust function of time)
class CustomThrustDirectionSettings: public ThrustDirectionGuidanceSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustDirectionFunction Function returning thrust direction unit vector as function fo time.
     */
    CustomThrustDirectionSettings(
            const boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction ):
        ThrustDirectionGuidanceSettings( custom_thrust_direction, "" ),
        thrustDirectionFunction_( thrustDirectionFunction ){ }

    //! Destructor.
    ~CustomThrustDirectionSettings( ){ }

    //! Function returning thrust direction unit vector as function fo time.
    boost::function< Eigen::Vector3d( const double ) > thrustDirectionFunction_;
};

//! Class for defining custom orientation of thrust (i.e. predefined body-fixed-to-propagation rotation as function of time)
/*!
 *  Class for defining custom orientation of thrust (i.e. predefined body-fixed-to-propagation rotation as function of time).
 *  Thrust is then computed from body-fixed direction of thrust (defined in ThrustEngineSettings).
 */
class CustomThrustOrientationSettings: public ThrustDirectionGuidanceSettings
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param thrustOrientationFunction Custom orientation of thrust (i.e. predefined body-fixed-to-propagation rotation
     * as function of time)
     */
    CustomThrustOrientationSettings(
            const boost::function< Eigen::Quaterniond( const double ) > thrustOrientationFunction ):
        ThrustDirectionGuidanceSettings( custom_thrust_orientation, "" ),
        thrustOrientationFunction_( thrustOrientationFunction ){ }

    //! Destructor.
    ~CustomThrustOrientationSettings( ){ }

    //! Custom orientation of thrust (i.e. predefined body-fixed-to-propagation rotation as function of time.
    boost::function< Eigen::Quaterniond( const double ) > thrustOrientationFunction_ ;
};

//! Function to create the object determining the direction of the thrust acceleration.
/*!
 * Function to create the object determining the direction of the thrust acceleration.
 * \param thrustDirectionGuidanceSettings Settings for thrust direction gudiance.
 * \param bodyMap List of pointers to body objects defining the full simulation environment.
 * \param nameOfBodyWithGuidance Name of body for which thrust guidane is to be created.
 * \param bodyFixedThrustOrientation Thrust direction in body-fixed frame.
 * \param magnitudeUpdateSettings Settings for the required updates to the environment during propagation. List is
 * extended by this function as needed.
 * \return Function determining the thrust direction in the propagation frame according to given requirements.
 */
boost::shared_ptr< propulsion::BodyFixedForceDirectionGuidance  > createThrustGuidanceModel(
        const boost::shared_ptr< ThrustDirectionGuidanceSettings > thrustDirectionGuidanceSettings,
        const NamedBodyMap& bodyMap,
        const std::string& nameOfBodyWithGuidance,
        const boost::function< Eigen::Vector3d( ) > bodyFixedThrustOrientation,
        std::map< propagators::EnvironmentModelsToUpdate, std::vector< std::string > >& magnitudeUpdateSettings );

//! List of available types of thust magnitude types
enum ThrustMagnitudeTypes
{
    constant_thrust_magnitude,
    from_engine_properties_thrust_magnitude,
    thrust_magnitude_from_time_function,
    thrust_magnitude_from_dependent_variables
};

//! Class defining settings for the thrust magnitude
/*!
 *  Class for providing settings the thrust magnitude of a single thrust model. This class is a functional (base) class for
 *  settings of thrust magnitude that require no information in addition to their type.
 *  Classes defining settings for thrust magnitude requiring additional information must be derived from this class.
 */
class ThrustEngineSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitudeGuidanceType Type of thrust magnitude guidance that is to be used
     * \param thrustOriginId Reference id of thrust origin that is to be used (empty if N/A).
     */
    ThrustEngineSettings(
            const ThrustMagnitudeTypes thrustMagnitudeGuidanceType,
            const std::string& thrustOriginId ):
        thrustMagnitudeGuidanceType_( thrustMagnitudeGuidanceType ),
        thrustOriginId_( thrustOriginId ){ }

    //! Destructor
    virtual ~ThrustEngineSettings( ){ }

    //! Type of thrust magnitude guidance that is to be used
    ThrustMagnitudeTypes thrustMagnitudeGuidanceType_;

    //! Reference id of thrust origin that is to be used (empty if N/A).
    std::string thrustOriginId_;
};

//! Class to define settigns for constant thrust settings.
class ConstantThrustEngineSettings: public ThrustEngineSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitude Constant thrust magnitude that is to be used.
     * \param specificImpulse Constant specific impulse that is to be used
     * \param bodyFixedThrustDirection Direction of thrust force in body-fixed frame (along longitudinal axis by default).
     */
    ConstantThrustEngineSettings(
            const double thrustMagnitude,
            const double specificImpulse,
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
        ThrustEngineSettings( constant_thrust_magnitude, "" ),
        thrustMagnitude_( thrustMagnitude ), specificImpulse_( specificImpulse ),
        bodyFixedThrustDirection_( bodyFixedThrustDirection ){ }

    //! Destructor
    ~ConstantThrustEngineSettings( ){ }

    //! Constant thrust magnitude that is to be used.
    double thrustMagnitude_;

    //! Constant specific impulse that is to be used
    double specificImpulse_;

    //! Direction of thrust force in body-fixed frame
    Eigen::Vector3d bodyFixedThrustDirection_;
};

//! Class to define thrust magnitude  to be taken directly from an engine model
class FromBodyThrustEngineSettings: public ThrustEngineSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param useAllEngines Boolean denoting whether all engines of the associated body are to be combined into a
     * single thrust magnitude.
     * \param thrustOrigin Name of engine model from which thrust is to be derived (must be empty if
     * useAllThrustModels is set to true)
     */
    FromBodyThrustEngineSettings(
            const bool useAllEngines = 1,
            const std::string& thrustOrigin = "" ):
        ThrustEngineSettings( from_engine_properties_thrust_magnitude, thrustOrigin ),
        useAllEngines_( useAllEngines ){ }

    //! Boolean denoting whether all engines of the associated body are to be combined into a single thrust magnitude
    bool useAllEngines_;
};

//! Class to define custom settings for thrust magnitude/specific impulse.
/*!
 * Class to define custom settings for thrust magnitude/specific impulse. Using this function, the thrust magnitude and
 * specific impulse are defined by arbitrary functions of time, the implementation for which is completely open to the user.
 * Also, a reset-function can be added, which is used to signal that a new time step is being computed (if applicable).
 * Note that for the definition of thrust and specific impulse as a function of a number independent variables with
 * clear physical meaning (e.g. dynamic pressure, Mach number, freestream density, etc.), the
 * ParameterizedThrustMagnitudeSettings settings object can be used.
 */
class FromFunctionThrustEngineSettings: public ThrustEngineSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param thrustMagnitudeFunction Function returning thrust magnitude as a function of time.
     * \param specificImpulseFunction Function returning specific impulse as a function of time.
     * \param isEngineOnFunction Function returning boolean denoting whether thrust is to be used (thrust and mass rate
     * set to zero if false, regardless of output of thrustMagnitudeFunction).
     * \param bodyFixedThrustDirection Direction of thrust force in body-fixed frame (along longitudinal axis by default).
     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
     * being started (empty by default)
     */
    FromFunctionThrustEngineSettings(
            const boost::function< double( const double ) > thrustMagnitudeFunction,
            const boost::function< double( const double ) > specificImpulseFunction,
            const boost::function< bool( const double ) > isEngineOnFunction = boost::lambda::constant( true ),
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ),
            const boost::function< void( const double ) > customThrustResetFunction = boost::function< void( const double ) >( ) ):
        ThrustEngineSettings( thrust_magnitude_from_time_function, "" ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        isEngineOnFunction_( isEngineOnFunction ),
        bodyFixedThrustDirection_( bodyFixedThrustDirection ),
        customThrustResetFunction_( customThrustResetFunction ){ }

    //! Destructor.
    ~FromFunctionThrustEngineSettings( ){ }

    //! Function returning thrust magnitude as a function of time.
    boost::function< double( const double) > thrustMagnitudeFunction_;

    //! Function returning specific impulse as a function of time.
    boost::function< double( const double) > specificImpulseFunction_;

    //! Function returning boolean denoting whether thrust is to be used.
    boost::function< bool( const double ) > isEngineOnFunction_;

    //! Direction of thrust force in body-fixed frame
    Eigen::Vector3d bodyFixedThrustDirection_;

    boost::function< void( const double ) > customThrustResetFunction_;
};

//! Interface function to multiply a maximum thrust by a multiplier to obtain the actual thrust
/*!
 * Interface function to multiply a maximum thrust by a multiplier to obtain the actual thrust
 * \param maximumThrustFunction Function returning the maxumum thrust as a function of a number of independent variables
 * \param maximumThrustMultiplier Function returning a value by which the output of maximumThrustFunction is to be multiplied
 * to obtain the actual thrust.
 * \param maximumThrustIndependentVariables List of variables to be passed as input to maximumThrustFunction
 * \return
 */
double multiplyMaximumThrustByScalingFactor(
        const boost::function< double( const std::vector< double >& ) > maximumThrustFunction,
        const boost::function< double( ) > maximumThrustMultiplier,
        const std::vector< double >& maximumThrustIndependentVariables );

//! Class to define the thrust magnitude and specific impulse as an interpolated function of N independent variables
/*!
 *  Class to define the thrust magnitude and specific impulse as an interpolated function of N independent variables.
 *  The physical meaning of the variables must be defined here, selecting from the options in the ThrustDependentVariables
 *  enum, and they are automatically retrieved from the relevant environment models during the propagation.
 *  Note that any number of user-specific functions may be included, as a  guidance_input_dependent_thrust type or
 *  maximum_thrust_multiplier. In the case of the maximum_thrust_multiplier, the thrustMagnitudeInterpolator input variable
 *  is assumed to define the maximum possible thrust, which is then multiplied by the function defining the
 *  maximum_thrust_multiplier.
 */
class ParameterizedThrustMagnitudeSettings: public ThrustEngineSettings
{
public:

    //! Constructor for parameterized thrust and specific impulse.
    /*!
     * Constructor, defines the interpolators for thrust and specific impulse, as well as the physical meaning of each of the
     * independent variables.
     * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
     * thrustDependentVariables contains an maximum_thrust_multiplier entry) as a function of the independent variables.
     * \param thrustDependentVariables List of identifiers for the physical meaning of each of the entries of the input to
     * the 'interpolate' function of thrustMagnitudeInterpolator.
     * \param specificImpulseInterpolator  Interpolator returning the current specific impulse as a function of the
     * independent variables.
     * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
     * input to the 'interpolate' function of specificImpulseInterpolator.
     * \param thrustGuidanceInputVariables List of functions returning user-defined guidance input variables for the thrust
     * (default none). The order of the functions in this vector is passed to the thrustMagnitudeInterpolator
     * in the order of the maximum_thrust_multiplier and guidance_input_dependent_thrust in the thrustDependentVariables
     * vector
     * \param specificImpulseGuidanceInputVariables List of functions returning user-defined guidance input variables
     * for the specific impulse (default none). The order of the functions in this vector is passed to the
     * specificImpulseInterpolator in the order of the and guidance_input_dependent_thrust in the thrustDependentVariables
     * vector
     * \param bodyFixedThrustDirection Direction of the thrust vector in the body-fixed frame (default in x-direction; to
     * vehicle front).
     */
    ParameterizedThrustMagnitudeSettings(
            const boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
            const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
            const boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
            const std::vector< propulsion::ThrustDependentVariables > specificImpulseDependentVariables,
            const std::vector< boost::function< double( ) > > thrustGuidanceInputVariables =
            std::vector< boost::function< double( ) > >( ),
            const std::vector< boost::function< double( ) > > specificImpulseGuidanceInputVariables =
            std::vector< boost::function< double( ) > >( ),
            const boost::function< void( const double) > inputUpdateFunction = boost::function< void( const double) >( ),
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
        ThrustEngineSettings( thrust_magnitude_from_dependent_variables, "" ),
        thrustMagnitudeFunction_( boost::bind( &interpolators::Interpolator< double, double >::interpolate,
                                               thrustMagnitudeInterpolator, _1 ) ),
        specificImpulseFunction_( boost::bind( &interpolators::Interpolator< double, double >::interpolate,
                                               specificImpulseInterpolator, _1 ) ),
        thrustDependentVariables_( thrustDependentVariables ),
        specificImpulseDependentVariables_( specificImpulseDependentVariables ),
        thrustGuidanceInputVariables_( thrustGuidanceInputVariables ),
        specificImpulseGuidanceInputVariables_( specificImpulseGuidanceInputVariables ),
        inputUpdateFunction_( inputUpdateFunction ),
        bodyFixedThrustDirection_( bodyFixedThrustDirection )
    {
        parseInputDataAndCheckConsistency( thrustMagnitudeInterpolator, specificImpulseInterpolator );
    }

    //! Constructor for parameterized thrust and constant specific impulse.
    /*!
     * Constructor, defines a constant thrust and an interpolator for thrust, as well as the physical meaning of each of the
     * independent variables.
     * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
     * thrustDependentVariables contains an maximum_thrust_multiplier entry) as a function of the independent variables.
     * \param thrustDependentVariables List of identifiers for the physical meaning of each of the entries of the input to
     * the 'interpolate' function of thrustMagnitudeInterpolator.
     * \param constantSpecificImpulse Constant specific impulse
     * \param thrustGuidanceInputVariables List of functions returning user-defined guidance input variables for the thrust
     * (default none). The order of the functions in this vector is passed to the thrustMagnitudeInterpolator
     * in the order of the maximum_thrust_multiplier and guidance_input_dependent_thrust in the thrustDependentVariables
     * vector
     * \param bodyFixedThrustDirection Direction of the thrust vector in the body-fixed frame (default in x-direction; to
     * vehicle front).
     */
    ParameterizedThrustMagnitudeSettings(
            const boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
            const std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables,
            const double constantSpecificImpulse,
            const std::vector< boost::function< double( ) > > thrustGuidanceInputVariables =
            std::vector< boost::function< double( ) > >( ),
            const boost::function< void( const double ) > inputUpdateFunction = boost::function< void( const double) >( ),
            const Eigen::Vector3d bodyFixedThrustDirection = Eigen::Vector3d::UnitX( ) ):
        ThrustEngineSettings( thrust_magnitude_from_dependent_variables, "" ),
        thrustMagnitudeFunction_( boost::bind( &interpolators::Interpolator< double, double >::interpolate,
                                               thrustMagnitudeInterpolator, _1 ) ),
        specificImpulseFunction_( boost::lambda::constant( constantSpecificImpulse ) ),
        thrustDependentVariables_( thrustDependentVariables ),
        thrustGuidanceInputVariables_( thrustGuidanceInputVariables ),
        inputUpdateFunction_( inputUpdateFunction ),
        bodyFixedThrustDirection_( bodyFixedThrustDirection )
    {
        parseInputDataAndCheckConsistency(
                    thrustMagnitudeInterpolator, boost::shared_ptr< interpolators::Interpolator< double, double > >( ) );
    }

    //! Function returning the thrust as a function of the independent variables.
    boost::function< double( const std::vector< double >& ) > thrustMagnitudeFunction_;

    //! Function returning the specific impulse as a function of the independent variables.
    boost::function< double( const std::vector< double >& ) > specificImpulseFunction_;

    //! List of identifiers for the physical meaning of each of the entries of the input to thrustMagnitudeFunction_.
    std::vector< propulsion::ThrustDependentVariables > thrustDependentVariables_;

    //! List of identifiers for the physical meaning of each of the entries of the input to specificImpulseDependentVariables_.
    std::vector< propulsion::ThrustDependentVariables > specificImpulseDependentVariables_;

    //! List of functions returning user-defined guidance input variables for the thrust
    std::vector< boost::function< double( ) > > thrustGuidanceInputVariables_;

    //! List of functions returning user-defined guidance input variables for the specific impulse
    std::vector< boost::function< double( ) > > specificImpulseGuidanceInputVariables_;

    boost::function< void( const double ) > inputUpdateFunction_;

    //! Direction of the thrust vector in the body-fixed frame
    Eigen::Vector3d bodyFixedThrustDirection_;

private:

    //! Function to check the validity of the input data, and process the maximum thrust multiplier if provided
    /*!
     *  Function to check the validity of the input data, and process the maximum thrust multiplier if provided.
     *  \param thrustMagnitudeInterpolator Interpolator for the (maximum) thrust provided to the constructor
     *  \param specificImpulseInterpolator Interpolator for the specific impulse provided to the constructor
     */
    void parseInputDataAndCheckConsistency(
            const boost::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
            const boost::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator );
};

} // namespace simulation_setup

} // namespace tudat
#endif // TUDAT_THRUSTSETTINGS_H
