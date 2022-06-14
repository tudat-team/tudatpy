/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *
 *    Kluever (2010), Low-Thrust Trajectory Optimization Using Orbital Averaging and Control Parameterization, In: Conway,
 *    (editor) Spacecraft trajectory optimization. Cambridge University Press, 2010.
 *    Boudestijn (2014), DEVELOPMENT OF A LOW -THRUST EARTH-CENTERED TRANSFER OPTIMIZER FOR THE PRELIMINARY MISSION DESIGN PHASE,
 *    M.Sc. Thesis, Delft University of Technology
 */

#ifndef TUDAT_THRUSTSETTINGS_H
#define TUDAT_THRUSTSETTINGS_H

#include <Eigen/Geometry>

#include <boost/bind/bind.hpp>
using namespace boost::placeholders;

#include <functional>

#include "tudat/astro/propulsion/thrustGuidance.h"
#include "tudat/astro/propulsion/thrustMagnitudeWrapper.h"
#include "tudat/astro/propagators/environmentUpdateTypes.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/astro/low_thrust/lowThrustLeg.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createFlightConditions.h"
#include "tudat/basics/deprecationWarnings.h"

namespace tudat
{

namespace simulation_setup
{

// List of available types of thrust direction guidance
//! @get_docstring(ThrustDirectionTypes.__docstring__)
enum ThrustDirectionTypes
{
    colinear_with_state_segment_thrust_direction,
    thrust_direction_from_existing_body_orientation,
    custom_thrust_direction,
    custom_thrust_orientation,
    mee_costate_based_thrust_direction

}; 

// Function to create a list of functions that (compute and) return independent variables for thrust
/*
 * Function to create a list of functions that (compute and) return independent variables for thrust and/or specific impulse.
 * This parameterization is used in the thrust mangitude type is thrust_magnitude_from_dependent_variables. This function
 * retrieves all input functions from the environment and a list of user-defined functions.
 * \param bodyWithGuidance Name of body for which the propulsion settings are to be retrieved.
 * \param independentVariables List of variables for which function returning them are to be created. Note that the number
 * of guidance_input_dependent_thrust entries must be equal to the size of guidanceInputFunctions. No entries of type
 * throttle_dependent_thrust are allowed.
 * \param guidanceInputFunctions Functions returning user-defined variables on which the thrust/specific impulse depends
 * \return List of functions that (compute and) return independent variables for thrust
 */
std::vector< std::function< double( ) > > getPropulsionInputVariables(
        const std::shared_ptr< Body > bodyWithGuidance,
        const std::vector< propulsion::ThrustIndependentVariables > independentVariables,
        const std::vector< std::function< double( ) > > guidanceInputFunctions =
        std::vector< std::function< double( ) > >( ) );


class ThrustDirectionSettings
{
public:

    ThrustDirectionSettings(
            const ThrustDirectionTypes thrustDirectionType,
            const std::string relativeBody = "" ):
        thrustDirectionType_( thrustDirectionType ), relativeBody_( relativeBody )
    {
        utilities::printDeprecationError(
                    "tudatpy.numerical_simulation.propagation_setup.thrust.ThrustDirectionSettings",
                    "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );

    }

    virtual ~ThrustDirectionSettings( ){ }
    ThrustDirectionTypes thrustDirectionType_;
    std::string relativeBody_;

};


class ThrustDirectionFromStateGuidanceSettings: public ThrustDirectionSettings
{
public:


    ThrustDirectionFromStateGuidanceSettings(
            const std::string& centralBody,
            const bool isColinearWithVelocity,
            const bool directionIsOppositeToVector ):
            ThrustDirectionSettings(colinear_with_state_segment_thrust_direction, centralBody ),
            isColinearWithVelocity_( isColinearWithVelocity ),
            directionIsOppositeToVector_( directionIsOppositeToVector )
    {
        utilities::printDeprecationError(
                    "tudatpy.numerical_simulation.propagation_setup.thrust.ThrustDirectionFromStateGuidanceSettings",
                    "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    }

    ~ThrustDirectionFromStateGuidanceSettings( ){ }

    bool isColinearWithVelocity_;
    bool directionIsOppositeToVector_;

};


class CustomThrustDirectionSettings: public ThrustDirectionSettings
{
public:


    CustomThrustDirectionSettings(
            const std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction ):
            ThrustDirectionSettings(custom_thrust_direction, "" ),
            thrustDirectionFunction_( thrustDirectionFunction )
    {
        utilities::printDeprecationError(
                    "tudatpy.numerical_simulation.propagation_setup.thrust.CustomThrustDirectionSettings",
                    "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    }

    ~CustomThrustDirectionSettings( ){ }
    std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction_;

};

class CustomThrustOrientationSettings: public ThrustDirectionSettings
{
public:

    CustomThrustOrientationSettings(
            const std::function< Eigen::Quaterniond( const double ) > thrustOrientationFunction ):
            ThrustDirectionSettings(custom_thrust_orientation, "" ),
            thrustOrientationFunction_( thrustOrientationFunction )
    {
        {
            utilities::printDeprecationError(
                        "tudatpy.numerical_simulation.propagation_setup.thrust.CustomThrustOrientationSettings",
                        "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
        }
    }

    CustomThrustOrientationSettings(
            const std::function< Eigen::Matrix3d( const double ) > thrustOrientationFunction ):
            ThrustDirectionSettings(custom_thrust_orientation, "" ),
            thrustOrientationFunction_(
            [=]( const double time ){ return Eigen::Quaterniond( thrustOrientationFunction( time ) ); } )
    {
        {
            utilities::printDeprecationError(
                        "tudatpy.numerical_simulation.propagation_setup.thrust.CustomThrustOrientationSettings",
                        "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
        }
    }

    ~CustomThrustOrientationSettings( ){ }
    std::function< Eigen::Quaterniond( const double ) > thrustOrientationFunction_;

};

//// Class for defining settings for MEE-costate based thrust direction guidance
///*
// *  Class for defining settings for MEE-costate based thrust direction guidance. Model details can be found in Kluever (2010) and
// *  Boudestijn (2014). The MEE-costates are provided for the five slow elements, as a function of time. Constructors for
// *  constant costates, and costates from an interpolator, are also provided.
// */
////! @get_docstring(MeeCostateBasedThrustDirectionSettings.__docstring__)
//class MeeCostateBasedThrustDirectionSettings: public ThrustDirectionSettings
//{
//public:

//    // Constructor with costate function
//    /*
//     * Constructor with costate function
//     * \param vehicleName Name of vehicle under thrust
//     * \param centralBodyName Name of central body (w.r.t. which MEE are calculated)
//     * \param costateFunction Function returning the 5 costates as a function of time
//     */
//    MeeCostateBasedThrustDirectionSettings(
//            const std::string& vehicleName,
//            const std::string& centralBodyName,
//            const std::function< Eigen::VectorXd( const double ) > costateFunction ):
//            ThrustDirectionSettings(mee_costate_based_thrust_direction, centralBodyName ),
//            vehicleName_( vehicleName ), costateFunction_( costateFunction ){ }

//    // Constructor with costate function
//    /*
//     * Constructor with costate function
//     * \param vehicleName Name of vehicle under thrust
//     * \param centralBodyName Name of central body (w.r.t. which MEE are calculated)
//     * \param costateInterpolator Interpolator returning the 5 costates as a function of time
//     */
//    MeeCostateBasedThrustDirectionSettings(
//            const std::string& vehicleName,
//            const std::string& centralBodyName,
//            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > costateInterpolator ):
//            ThrustDirectionSettings(mee_costate_based_thrust_direction, centralBodyName ),
//            vehicleName_( vehicleName ),
//            costateFunction_(
//            std::bind( static_cast< Eigen::VectorXd( interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd >::* )
//                       ( const double ) >( &interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd >::interpolate ),
//                       costateInterpolator, std::placeholders::_1 ) ){ }

//    // Constructor with costate function
//    /*
//     * Constructor with costate function
//     * \param vehicleName Name of vehicle under thrust
//     * \param centralBodyName Name of central body (w.r.t. which MEE are calculated)
//     * \param constantCostates The 5 costates, which will be used as constants in time
//     */
//    MeeCostateBasedThrustDirectionSettings(
//            const std::string& vehicleName,
//            const std::string& centralBodyName,
//            const Eigen::VectorXd constantCostates ):
//            ThrustDirectionSettings(mee_costate_based_thrust_direction, centralBodyName ),
//            vehicleName_( vehicleName ), costateFunction_( [ = ]( const double ){ return constantCostates; } ){ }


//    // Destructor.
//    ~MeeCostateBasedThrustDirectionSettings( ){ }

//    // Name of vehicle under thrust
//    std::string vehicleName_;

//    // Function returning the 5 costates as a function of time
//    std::function< Eigen::VectorXd( const double ) > costateFunction_;

//};

//! @get_docstring(thrustDirectionFromStateGuidanceSettings)
inline std::shared_ptr< ThrustDirectionSettings > thrustDirectionFromStateGuidanceSettings(
        const std::string& centralBody,
        const bool isColinearWithVelocity,
        const bool directionIsOppositeToVector  )
{
    return std::make_shared< ThrustDirectionFromStateGuidanceSettings >(
                centralBody, isColinearWithVelocity, directionIsOppositeToVector );
}

//! @get_docstring(thrustFromExistingBodyOrientation)
inline std::shared_ptr< ThrustDirectionSettings > thrustFromExistingBodyOrientation(  )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.propagation_setup.thrust.thrust_from_existing_body_orientation",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    return nullptr;
}


//! @get_docstring(customThrustOrientationSettings, 2)
inline std::shared_ptr< ThrustDirectionSettings > customThrustOrientationSettings(
        const std::function< Eigen::Matrix3d( const double ) > thrustOrientationFunction  )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.propagation_setup.thrust.custom_thrust_orientation",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    return nullptr;

}

//! @get_docstring(customThrustDirectionSettings)
inline std::shared_ptr< ThrustDirectionSettings > customThrustDirectionSettings(
        const std::function< Eigen::Vector3d( const double ) > thrustDirectionFunction  )
{
    utilities::printDeprecationError(
                "tudatpy.numerical_simulation.propagation_setup.thrust.custom_thrust_direction",
                "https://docs.tudat.space/en/stable/_src_user_guide/state_propagation/environment_setup/thrust_refactor/thrust_refactor.html#thrust-acceleration" );
    return nullptr;

}

// List of available types of thrust magnitude types
//! @get_docstring(ThrustMagnitudeTypes.__docstring__)
enum ThrustMagnitudeTypes
{
    constant_thrust_magnitude,
    from_engine_properties_thrust_magnitude,
    thrust_magnitude_from_time_function,
    thrust_magnitude_from_dependent_variables,
    bang_bang_thrust_magnitude_from_mee_costates
};

// Class defining settings for the thrust magnitude
/*
 *  Class for providing settings the thrust magnitude of a single thrust model. This class is a functional (base) class for
 *  settings of thrust magnitude that require no information in addition to their type.
 *  Classes defining settings for thrust magnitude requiring additional information must be derived from this class.
 */
//! @get_docstring(ThrustMagnitudeSettings.__docstring__)
class ThrustMagnitudeSettings
{
public:

    // Constructor
    /*
     * Constructor
     * \param thrustMagnitudeType Type of thrust magnitude guidance that is to be used
     * \param thrustOriginId Reference id of thrust origin that is to be used (empty if N/A).
     */
    ThrustMagnitudeSettings(
            const ThrustMagnitudeTypes thrustMagnitudeType,
            const std::string& thrustOriginId ):
            thrustMagnitudeType_(thrustMagnitudeType ),
            thrustOriginId_( thrustOriginId ){ }

    // Destructor
    virtual ~ThrustMagnitudeSettings( ){ }

    // Type of thrust magnitude guidance that is to be used
    ThrustMagnitudeTypes thrustMagnitudeType_;

    // Reference id of thrust origin that is to be used (empty if N/A).
    std::string thrustOriginId_;

};

// Class to define settigns for constant thrust settings.
//! @get_docstring(ConstantThrustMagnitudeSettings.__docstring__)
class ConstantThrustMagnitudeSettings: public ThrustMagnitudeSettings
{
public:

    // Constructor
    /*
     * Constructor
     * \param thrustMagnitude Constant thrust magnitude that is to be used.
     * \param specificImpulse Constant specific impulse that is to be used
     * \param bodyFixedThrustDirection Direction of thrust force in body-fixed frame (along longitudinal axis by default).
     */
    ConstantThrustMagnitudeSettings(
            const double thrustMagnitude,
            const double specificImpulse ):
        ThrustMagnitudeSettings( constant_thrust_magnitude, "" ),
        thrustMagnitude_( thrustMagnitude ), specificImpulse_( specificImpulse ){ }

    // Destructor
    ~ConstantThrustMagnitudeSettings( ){ }


    // Constant thrust magnitude that is to be used.
    double thrustMagnitude_;

    // Constant specific impulse that is to be used
    double specificImpulse_;

};

//! @get_docstring(FromFunctionThrustMagnitudeSettings.__docstring__)
class CustomThrustMagnitudeSettings: public ThrustMagnitudeSettings
{
public:

    // Constructor
    /*
     * Constructor
     * \param thrustMagnitudeFunction Function returning thrust magnitude as a function of time.
     * \param specificImpulseFunction Function returning specific impulse as a function of time.
     * \param isEngineOnFunction Function returning boolean denoting whether thrust is to be used (thrust and mass rate
     * set to zero if false, regardless of output of thrustMagnitudeFunction).
     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
     * being started (empty by default)
     */
    CustomThrustMagnitudeSettings(
            const std::function< double( const double ) > thrustMagnitudeFunction,
            const std::function< double( const double ) > specificImpulseFunction,
            const bool inputIsForce = true ):
        ThrustMagnitudeSettings( thrust_magnitude_from_time_function, "" ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( specificImpulseFunction ),
        specificImpulseIsConstant_( false ),
        inputIsForce_( inputIsForce ){ }

    CustomThrustMagnitudeSettings(
            const std::function< double( const double ) > thrustMagnitudeFunction,
            const double specificImpulse,
            const bool inputIsForce = true ):
        ThrustMagnitudeSettings( thrust_magnitude_from_time_function, "" ),
        thrustMagnitudeFunction_( thrustMagnitudeFunction ),
        specificImpulseFunction_( [=](const double){return specificImpulse;} ),
        specificImpulseIsConstant_( true ),
        inputIsForce_( inputIsForce ){ }

    // Destructor.
    ~CustomThrustMagnitudeSettings( ){ }

    // Function returning thrust magnitude as a function of time.
    std::function< double( const double) > thrustMagnitudeFunction_;

    // Function returning specific impulse as a function of time.
    std::function< double( const double) > specificImpulseFunction_;

    bool specificImpulseIsConstant_;

    bool inputIsForce_;
};


//// Class for defining settings for bang-bang MEE-costate based thrust magnitude.
///*
// *  Class for defining settings for bang-bang MEE-costate based thrust magnitude. Model details can be found in Kluever (2010) and
// *  Boudestijn (2014). The MEE-costates are provided for the five slow elements, as a function of time. Constructors for
// *  constant costates, and costates from an interpolator, are also provided.
// */
//class FromMeeCostatesBangBangThrustMagnitudeSettings: public ThrustMagnitudeSettings
//{
//public:

//    // Constructor with costates function.
//    /*
//     * Constructor with costates function.
//     * \param thrustMagnitude Maximum thrust magnitude (if engine is on).
//     * \param specificImpulseFunction Function returning specific impulse as a function of time.
//     * \param costateFunction Function returning the MEE 5 costates as a function of time.
//     * \param vehicleName Name of vehicle under thrust.
//     * \param centralBodyName Name of central body (w.r.t. which MEE are calculated).
//     * \param bodyFixedThrustDirection Direction of thrust force in body-fixed frame (along longitudinal axis by default).
//     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
//     * being started (empty by default).
//     *
//     */
//    FromMeeCostatesBangBangThrustMagnitudeSettings(
//            const double thrustMagnitude,
//            const std::function< double( const double ) > specificImpulseFunction,
//            const std::function< Eigen::VectorXd( const double ) > costatesFunction,
//            const std::string& vehicleName,
//            const std::string& centralBodyName,
//            const std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = [ ]( ){ return  Eigen::Vector3d::UnitX( ); },
//            const std::function< void( const double ) > customThrustResetFunction = std::function< void( const double ) >( ) ):
//        ThrustMagnitudeSettings( bang_bang_thrust_magnitude_from_mee_costates, "" ),
//        maximumThrustMagnitude_( thrustMagnitude ), specificImpulseFunction_( specificImpulseFunction ),
//        costatesFunction_( costatesFunction ), vehicleName_( vehicleName ), centralBodyName_( centralBodyName ),
//        bodyFixedThrustDirection_( bodyFixedThrustDirection ), customThrustResetFunction_( customThrustResetFunction ){ }

//    // Constructor with costate function
//    /*
//     * Constructor with costate function
//     * \param thrustMagnitude Maximum thrust magnitude (if engine is on).
//     * \param specificImpulseFunction Function returning specific impulse as a function of time.
//     * \param costateInterpolator Interpolator returning the 5 MEE costates as a function of time.
//     * \param vehicleName Name of vehicle under thrust.
//     * \param centralBodyName Name of central body (w.r.t. which MEE are calculated).
//     * \param bodyFixedThrustDirection Direction of thrust force in body-fixed frame (along longitudinal axis by default).
//     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
//     * being started (empty by default).
//     */
//    FromMeeCostatesBangBangThrustMagnitudeSettings(
//            const double thrustMagnitude,
//            const std::function< double( const double ) > specificImpulseFunction,
//            const std::string& vehicleName,
//            const std::string& centralBodyName,
//            const std::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd > > costateInterpolator,
//            const std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = [ ]( ){ return  Eigen::Vector3d::UnitX( ); },
//            const std::function< void( const double ) > customThrustResetFunction = std::function< void( const double ) >( ) ):
//        ThrustMagnitudeSettings( bang_bang_thrust_magnitude_from_mee_costates, "" ),
//        maximumThrustMagnitude_( thrustMagnitude ), specificImpulseFunction_( specificImpulseFunction ),
//        costatesFunction_( std::bind( static_cast< Eigen::VectorXd( interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd >::* )
//                                      ( const double ) >( &interpolators::OneDimensionalInterpolator< double, Eigen::VectorXd >::interpolate ),
//                                      costateInterpolator, std::placeholders::_1 ) ),
//        vehicleName_( vehicleName ), centralBodyName_( centralBodyName ), bodyFixedThrustDirection_( bodyFixedThrustDirection ),
//        customThrustResetFunction_( customThrustResetFunction ){ }

//    // Constructor with costate function
//    /*
//     * Constructor with costate function
//     * \param thrustMagnitude Maximum thrust magnitude (if engine is on).
//     * \param specificImpulseFunction Function returning specific impulse as a function of time.
//     * \param costateInterpolator The MEE 5 costates, which will be used as constants in time.
//     * \param vehicleName Name of vehicle under thrust.
//     * \param centralBodyName Name of central body (w.r.t. which MEE are calculated).
//     * \param bodyFixedThrustDirection Direction of thrust force in body-fixed frame (along longitudinal axis by default).
//     * \param customThrustResetFunction Custom function that is to be called when signalling that a new time step is
//     * being started (empty by default).
//     */
//    FromMeeCostatesBangBangThrustMagnitudeSettings(
//            const double thrustMagnitude,
//            const std::function< double( const double ) > specificImpulseFunction,
//            const std::string& vehicleName,
//            const std::string& centralBodyName,
//            const Eigen::VectorXd constantCostates,
//            const std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection = [ ]( ){ return  Eigen::Vector3d::UnitX( ); },
//            const std::function< void( const double ) > customThrustResetFunction = std::function< void( const double ) >( ) ):
//        ThrustMagnitudeSettings( bang_bang_thrust_magnitude_from_mee_costates, "" ),
//        maximumThrustMagnitude_( thrustMagnitude ), specificImpulseFunction_( specificImpulseFunction ),
//        costatesFunction_( [ = ]( const double ){ return constantCostates; } ),
//        vehicleName_( vehicleName ), centralBodyName_( centralBodyName ), bodyFixedThrustDirection_( bodyFixedThrustDirection ),
//        customThrustResetFunction_( customThrustResetFunction ) { }


//    // Destructor.
//    ~FromMeeCostatesBangBangThrustMagnitudeSettings( ){ }

//    // Maximum thrust magnitude.
//    double maximumThrustMagnitude_;

//    // Function returning specific impulse as a function of time.
//    std::function< double( const double) > specificImpulseFunction_;

//    // Vector of Mee costates.
//    std::function< Eigen::VectorXd( const double ) > costatesFunction_;

//    // Name of the body under thrust.
//    std::string vehicleName_;

//    // Name of the central body.
//    std::string centralBodyName_;

//    // Direction of thrust force in body-fixed frame.
//    std::function< Eigen::Vector3d( ) > bodyFixedThrustDirection_;

//    std::function< void( const double ) > customThrustResetFunction_;
//};


//! @get_docstring(constantThrustMagnitudeSettings)
inline std::shared_ptr< ThrustMagnitudeSettings > constantThrustMagnitudeSettings(
        const double thrustMagnitude,
        const double specificImpulse )
{
    return std::make_shared< ConstantThrustMagnitudeSettings >(
                thrustMagnitude, specificImpulse );
}

//! @get_docstring(fromFunctionThrustMagnitudeSettings)
inline std::shared_ptr< ThrustMagnitudeSettings > fromFunctionThrustMagnitudeSettings(
        const std::function< double( const double ) > thrustMagnitudeFunction,
        const std::function< double( const double ) > specificImpulseFunction )
{
    return std::make_shared< CustomThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulseFunction );
}

inline std::shared_ptr< ThrustMagnitudeSettings > fromFunctionThrustMagnitudeFixedIspSettings(
        const std::function< double( const double ) > thrustMagnitudeFunction,
        const double specificImpulse )
{
    return std::make_shared< CustomThrustMagnitudeSettings >(
                thrustMagnitudeFunction, specificImpulse );
}

//! @get_docstring(fromFunctionThrustMagnitudeSettings)
inline std::shared_ptr< ThrustMagnitudeSettings > customThrustAccelerationMagnitudeSettings(
        const std::function< double( const double ) > thrustAccelerationMagnitudeFunction,
        const std::function< double( const double ) > specificImpulseFunction )
{
    return std::make_shared< CustomThrustMagnitudeSettings >(
                thrustAccelerationMagnitudeFunction, specificImpulseFunction, false );
}

//! @get_docstring(fromFunctionThrustMagnitudeSettings)
inline std::shared_ptr< ThrustMagnitudeSettings > customThrustAccelerationMagnitudeFixedIspSettings(
        const std::function< double( const double ) > thrustAccelerationMagnitudeFunction,
        const double specificImpulse )
{
    return std::make_shared< CustomThrustMagnitudeSettings >(
                thrustAccelerationMagnitudeFunction, specificImpulse, false );
}






// Interface function to multiply a maximum thrust by a multiplier to obtain the actual thrust
/*
         * Interface function to multiply a maximum thrust by a multiplier to obtain the actual thrust
         * \param maximumThrustFunction Function returning the maxumum thrust as a function of a number of independent variables
         * \param maximumThrustMultiplier Function returning a value by which the output of maximumThrustFunction is to be multiplied
         * to obtain the actual thrust.
         * \param maximumThrustIndependentVariables List of variables to be passed as input to maximumThrustFunction
         * \return
         */
double multiplyMaximumThrustByScalingFactor(
        const std::function< double( const std::vector< double >& ) > maximumThrustFunction,
        const std::function< double( ) > maximumThrustMultiplier,
        const std::vector< double >& maximumThrustIndependentVariables );

// Interface base class that can be defined by user to make the use of the ParameterizedThrustMagnitudeSettings class easier
/*
 *  For a specific user-defined guidance approach, a derived class must be defined, that computes the guidance/throttle
 *  commands at each time step/vehicle state. Using the createParameterizedThrustMagnitudeSettings functions, a
 *  ParameterizedThrustMagnitudeSettings can then be made.
 *  In the derived class, the user must define an updateGuidanceParameters function, which sets the entries of the
 *  currentThrustGuidanceParameters_ and currentSpecificImpulseParameters_ vectors (size must be compatible with
 *  ThrustInputParameterGuidance constructor). Note that the throttle setting is simply defined as one of the entries of
 *  currentThrustGuidanceParameters_.
 */
class ThrustInputParameterGuidance
{
public:

    // Constructor
    /*
     *  Constructor
     *  \param numberOfThrustInputParameters Number of guidance-input parameters that are computed for the thrust
     *  \param numberOfSpecificImpulseInputParameters Number of guidance-input parameters that are computed for the
     *  specific impulse
     *  \param includeThrottleSetting Boolean that denotes whether any of the guidance-input parameters for the thrust
     *  define a throttle setting
     *  \param throttleSettingIndex Index in list of thrust guidance parameters that denotes the throttle setting (if any)
     */
    ThrustInputParameterGuidance(
            const int numberOfThrustInputParameters,
            const int numberOfSpecificImpulseInputParameters,
            const bool includeThrottleSetting = false,
            const int throttleSettingIndex = -1 ):
        numberOfThrustInputParameters_( numberOfThrustInputParameters ),
        numberOfSpecificImpulseInputParameters_( numberOfSpecificImpulseInputParameters ),
        includeThrottleSetting_( includeThrottleSetting ), throttleSettingIndex_( throttleSettingIndex )
    {
        if( includeThrottleSetting && ( throttleSettingIndex < 0 || throttleSettingIndex >= numberOfThrustInputParameters ) )
        {
            throw std::runtime_error( "Error when creating ThrustInputParameterGuidance, inconsistent throttle input" );
        }
        else if( !includeThrottleSetting && throttleSettingIndex >= 0 )
        {
            throw std::runtime_error( "Error when creating ThrustInputParameterGuidance, no throttle input with positive throttle index" );
        }

        currentThrustGuidanceParameters_.resize( numberOfThrustInputParameters_ );
        currentSpecificImpulseParameters_.resize( numberOfSpecificImpulseInputParameters_ );
    }

    // Destructor
    virtual ~ThrustInputParameterGuidance( ){ }

    // Function to retrieve the number of guidance-input parameters that are computed for the thrust
    /*
     * Function to retrieve the number of guidance-input parameters that are computed for the thrust
     * \return Number of guidance-input parameters that are computed for the thrust
     */
    int getNumberOfThrustInputParameters( )
    {
        return numberOfThrustInputParameters_;
    }

    // Function to retrieve the number of guidance-input parameters that are computed for the specific impulse
    /*
     * Function to retrieve the number of guidance-input parameters that are computed for the specific impulse
     * \return Number of guidance-input parameters that are computed for the specific impulse
     */
    int getNumberOfSpecificImpulseInputParameters( )
    {
        return numberOfSpecificImpulseInputParameters_;
    }

    // Function to retrieve a guidance-input parameter that for the thrust
    /*
     * Function to retrieve a guidance-input parameter that for the thrust
     * \param inputParameterIndex Entry of guidance-input vector for thrust that is to be retrieved.
     * \return Guidance-input parameters that for the thrust at given index
     */
    double getThrustInputGuidanceParameter( const int inputParameterIndex )
    {
        return currentThrustGuidanceParameters_.at( inputParameterIndex );
    }

    // Function to retrieve a guidance-input parameter that for the specific impulse
    /*
     * Function to retrieve a guidance-input parameter that for the specific impulse
     * \param inputParameterIndex Entry of guidance-input vector for specific impulse that is to be retrieved.
     * \return Guidance-input parameters that for the specific impulse at given index
     */
    double getSpecificImpulseInputGuidanceParameter( const int inputParameterIndex )
    {
        return currentSpecificImpulseParameters_.at( inputParameterIndex );
    }

    // Function to retrieve whether any of the guidance-input parameters for the thrust define a throttle setting
    /*
     * Function to retrieve whether any of the guidance-input parameters for the thrust define a throttle setting
     * \return Boolean that denotes whether any of the guidance-input parameters for the thrust define a throttle setting.
     */
    bool getIncludeThrottleSetting( )
    {
        return includeThrottleSetting_;
    }

    // Function to retrieve the index in list of thrust guidance parameters that denotes the throttle setting
    /*
     * Function to retrieve the index in list of thrust guidance parameters that denotes the throttle setting (-1 if none)
     * \return Index in list of thrust guidance parameters that denotes the throttle setting
     */
    int getThrottleSettingIndex( )
    {
        return throttleSettingIndex_;
    }

    // Pure virtual function that updates the guidance algorithm to the current time/state
    virtual void updateGuidanceParameters( ) = 0;

    // Update function for the guidance algorithm
    /*
     * Update function for the guidance algorithm, calls the updateGuidanceParameters that is to be implemented in derived
     * class. Calling this function with NaN time signals that a new time step has commenced.
     * \param time Time to which guidance object is to be updated
     */
    void update( const double time )
    {
        if( time != time )
        {
            currentTime_ = time;
        }
        else
        {
            if( !( currentTime_ == time ) )
            {
                currentTime_ =  time;
                updateGuidanceParameters( );
            }
        }
    }

protected:

    // Vector of guidance-input parameters that are computed for the thrust
    std::vector< double > currentThrustGuidanceParameters_;

    // Vector of guidance-input parameters that are computed for the specific impulse
    std::vector< double > currentSpecificImpulseParameters_;

    // Number of guidance-input parameters that are computed for the thrust
    int numberOfThrustInputParameters_;

    // Number of guidance-input parameters that are computed for the specific impulse
    int numberOfSpecificImpulseInputParameters_;

    // Boolean that denotes whether any of the guidance-input parameters for the thrust define a throttle setting
    bool includeThrottleSetting_;

    // Index in list of thrust guidance parameters that denotes the throttle setting (if any)
    int throttleSettingIndex_;

    // Current time of interface object (e.g. time to which object was last updated).
    double currentTime_;

};

// Class to compute throttling law for the parameterized thrust magnitude, where the throttle is determined from a maximum
// axial g-load.
/*
 *  Class to compute throttling law for the parameterized thrust magnitude, where the throttle is determined from a maximum,
 *  the user must provide the interpolator for the (maximum) thrust, as well as the associated physical meaning of the
 *  independent variables. Also, a maximum axial acceleration must be provided. If the current thrust results in an
 *  acceleration less than this limit, the throttle is set to 1, if it is higher than the maximum, the throttle is set such
 *  that the acceleration is on this limit.
 */
class AccelerationLimitedThrottleGuidance: public ThrustInputParameterGuidance
{
public:

    // Constructor.
    /*
     * Constructor
     * \param bodies List of pointers to body objects defining the full simulation environment.
     * \param nameOfBodyWithGuidance Name of body for which thrust guidance is to be created.
     * \param nameOfCentralBody Name of body w.r.t. which thrust guidance is computed (e.g. Earth if the altitude from Earth
     * is used as an independent variable of the thrust).
     * \param independentVariables Physical meaning of each of the independent variables used as input to the
     * thrustInterpolator.
     * \param thrustInterpolator Interpolator that computes the maximum thrust as a function of the independent variables.
     * \param maximumAcceleration Maxmum allowable acceleration due to the thrust force.
     */
    AccelerationLimitedThrottleGuidance(
            const SystemOfBodies& bodies,
            const std::string nameOfBodyWithGuidance,
            const std::string nameOfCentralBody,
            const std::vector< propulsion::ThrustIndependentVariables > independentVariables,
            const std::shared_ptr< interpolators::Interpolator< double, double > > thrustInterpolator,
            const double maximumAcceleration ): ThrustInputParameterGuidance( 1, 0, true, 0 ),
        bodyWithGuidance_( bodies.at( nameOfBodyWithGuidance ) ), thrustInterpolator_( thrustInterpolator ),
        maximumAcceleration_( maximumAcceleration )
    {
        // Split independent variables into environmental/guidance.
        int numberOfThrottles = 0;
        std::vector< propulsion::ThrustIndependentVariables > guidanceFreeIndependentVariables;
        for( unsigned int i = 0; i < independentVariables.size( ); i++ )
        {
            if( independentVariables.at( i ) == propulsion::throttle_dependent_thrust )
            {
                numberOfThrottles ++;
            }
            else if( independentVariables.at( i ) == propulsion::guidance_input_dependent_thrust )
            {
                throw std::runtime_error( "Error in AccelerationLimitedThrottleGuidance, guidance input not supported" );
            }
            else
            {
                guidanceFreeIndependentVariables.push_back( independentVariables.at( i ) );
            }
        }

        // Check input consistency.
        if( numberOfThrottles > 1 )
        {
            throw std::runtime_error( "Error in AccelerationLimitedThrottleGuidance, multiple throttles detected" );
        }

        // Create thrust input.
        if( bodyWithGuidance_->getFlightConditions( ) == nullptr && nameOfCentralBody != "" )
        {
            bodyWithGuidance_->setFlightConditions(
                        createAtmosphericFlightConditions( bodyWithGuidance_,
                                                           bodies.at( nameOfCentralBody ),
                                                           nameOfBodyWithGuidance,
                                                           nameOfCentralBody ) );
        }
        thrustInputFunctions_ = getPropulsionInputVariables(
                    bodyWithGuidance_, guidanceFreeIndependentVariables );

        currentThrustInput_.resize( thrustInputFunctions_.size( ) );
    }

    // Function that updates the guidance algorithm to the current time/state: sets the throttle value based on axia
    /*
     *  Function that updates the guidance algorithm to the current time/state: sets the throttle value based on axial
     *  acceleration when using the maximum thrust.
     */
    void updateGuidanceParameters( )
    {
        // Compute environmental input variables for thrust.
        for( unsigned int i = 0; i < thrustInputFunctions_.size( ); i++ )
        {
            currentThrustInput_[ i ] = thrustInputFunctions_.at( i )( );
        }

        // Compute maximum thrust
        double currentThrust = thrustInterpolator_->interpolate( currentThrustInput_ );
        double currentMass = bodyWithGuidance_->getBodyMass( );

        // Set throttle (to < 1 if necessary)
        if( currentThrust / currentMass < maximumAcceleration_ )
        {
            currentThrustGuidanceParameters_[ 0 ] = 1.0;
        }
        else
        {
            currentThrustGuidanceParameters_[ 0 ] = maximumAcceleration_ * currentMass / currentThrust;
        }
    }

private:

    // Body for which thrust guidance is used.
    std::shared_ptr< simulation_setup::Body > bodyWithGuidance_;

    // Interpolator that computes the maximum thrust as a function of the independent variables.
    std::shared_ptr< interpolators::Interpolator< double, double > > thrustInterpolator_;

    // Maxmum allowable acceleration due to the thrust force.
    double maximumAcceleration_;

    // Functions to compute the current entries of the input to thrustInterpolator_.
    std::vector< std::function< double( ) > >  thrustInputFunctions_;

    // Pre-declared vector used as input to thrustInterpolator_.
    std::vector< double >  currentThrustInput_;

};

// Class to define the thrust magnitude and specific impulse as an interpolated function of N independent variables
/*
 *  Class to define the thrust magnitude and specific impulse as an interpolated function of N independent variables.
 *  The physical meaning of the variables must be defined here, selecting from the options in the ThrustIndependentVariables
 *  enum, and they are automatically retrieved from the relevant environment models during the propagation.
 *  Note that any number of user-specific functions may be included, as a  guidance_input_dependent_thrust type or
 *  throttle_dependent_thrust. In the case of the throttle_dependent_thrust, the thrustMagnitudeInterpolator input variable
 *  is assumed to define the maximum possible thrust, which is then multiplied by the function defining the
 *  throttle_dependent_thrust.
 */
class ParameterizedThrustMagnitudeSettings: public ThrustMagnitudeSettings
{
public:

    // Constructor for parameterized thrust and specific impulse.
    /*
     * Constructor, defines the interpolators for thrust and specific impulse, as well as the physical meaning of each of the
     * independent variables.
     * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
     * thrustIndependentVariables contains an throttle_dependent_thrust entry) as a function of the independent variables.
     * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
     * the 'interpolate' function of thrustMagnitudeInterpolator.
     * \param specificImpulseInterpolator  Interpolator returning the current specific impulse as a function of the
     * independent variables.
     * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
     * input to the 'interpolate' function of specificImpulseInterpolator.
     * \param thrustGuidanceInputVariables List of functions returning user-defined guidance input variables for the thrust
     * (default none). The order of the functions in this vector is passed to the thrustMagnitudeInterpolator
     * in the order of the throttle_dependent_thrust and guidance_input_dependent_thrust in the thrustIndependentVariables
     * vector
     * \param specificImpulseGuidanceInputVariables List of functions returning user-defined guidance input variables
     * for the specific impulse (default none). The order of the functions in this vector is passed to the
     * specificImpulseInterpolator in the order of the and guidance_input_dependent_thrust in the thrustIndependentVariables
     * vector
     * \param inputUpdateFunction Function that is called to update the user-defined guidance to the current time
     * (empty by default).
     */
    ParameterizedThrustMagnitudeSettings(
            const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
            const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
            const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
            const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
            const std::vector< std::function< double( ) > > thrustGuidanceInputVariables =
            std::vector< std::function< double( ) > >( ),
            const std::vector< std::function< double( ) > > specificImpulseGuidanceInputVariables =
            std::vector< std::function< double( ) > >( ),
            const std::function< void( const double) > inputUpdateFunction = std::function< void( const double) >( ) ):
        ThrustMagnitudeSettings( thrust_magnitude_from_dependent_variables, "" ),
        thrustMagnitudeFunction_( std::bind( &interpolators::Interpolator< double, double >::interpolate,
                                             thrustMagnitudeInterpolator, std::placeholders::_1 ) ),
        specificImpulseFunction_( std::bind( &interpolators::Interpolator< double, double >::interpolate,
                                             specificImpulseInterpolator, std::placeholders::_1 ) ),
        thrustIndependentVariables_( thrustIndependentVariables ),
        specificImpulseDependentVariables_( specificImpulseDependentVariables ),
        thrustGuidanceInputVariables_( thrustGuidanceInputVariables ),
        specificImpulseGuidanceInputVariables_( specificImpulseGuidanceInputVariables ),
        inputUpdateFunction_( inputUpdateFunction )
    {
        parseInputDataAndCheckConsistency( thrustMagnitudeInterpolator, specificImpulseInterpolator );
    }

    // Constructor for parameterized thrust and constant specific impulse.
    /*
     * Constructor, defines a constant thrust and an interpolator for thrust, as well as the physical meaning of each of the
     * independent variables.
     * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
     * thrustIndependentVariables contains an throttle_dependent_thrust entry) as a function of the independent variables.
     * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
     * the 'interpolate' function of thrustMagnitudeInterpolator.
     * \param constantSpecificImpulse Constant specific impulse
     * \param thrustGuidanceInputVariables List of functions returning user-defined guidance input variables for the thrust
     * (default none). The order of the functions in this vector is passed to the thrustMagnitudeInterpolator
     * in the order of the throttle_dependent_thrust and guidance_input_dependent_thrust in the thrustIndependentVariables
     * vector
     * \param inputUpdateFunction Function that is called to update the user-defined guidance to the current time
     * (empty by default).
     */
    ParameterizedThrustMagnitudeSettings(
            const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
            const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
            const double constantSpecificImpulse,
            const std::vector< std::function< double( ) > > thrustGuidanceInputVariables =
            std::vector< std::function< double( ) > >( ),
            const std::function< void( const double ) > inputUpdateFunction = std::function< void( const double) >( ) ):
        ThrustMagnitudeSettings( thrust_magnitude_from_dependent_variables, "" ),
        thrustMagnitudeFunction_( std::bind( &interpolators::Interpolator< double, double >::interpolate,
                                             thrustMagnitudeInterpolator, std::placeholders::_1 ) ),
        specificImpulseFunction_( [ = ]( const std::vector< double >& ){ return constantSpecificImpulse; } ),
        thrustIndependentVariables_( thrustIndependentVariables ),
        thrustGuidanceInputVariables_( thrustGuidanceInputVariables ),
        inputUpdateFunction_( inputUpdateFunction )
    {
        parseInputDataAndCheckConsistency(
                    thrustMagnitudeInterpolator, std::shared_ptr< interpolators::Interpolator< double, double > >( ) );
    }

    // Function returning the thrust as a function of the independent variables.
    std::function< double( const std::vector< double >& ) > thrustMagnitudeFunction_;

    // Function returning the specific impulse as a function of the independent variables.
    std::function< double( const std::vector< double >& ) > specificImpulseFunction_;

    // List of identifiers for the physical meaning of each of the entries of the input to thrustMagnitudeFunction_.
    std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables_;

    // List of identifiers for the physical meaning of each of the entries of the input to specificImpulseDependentVariables_.
    std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables_;

    // List of functions returning user-defined guidance input variables for the thrust
    std::vector< std::function< double( ) > > thrustGuidanceInputVariables_;

    // List of functions returning user-defined guidance input variables for the specific impulse
    std::vector< std::function< double( ) > > specificImpulseGuidanceInputVariables_;

    // Function that is called to update the user-defined guidance to the current time
    std::function< void( const double ) > inputUpdateFunction_;

private:

    // Function to check the validity of the input data, and process the maximum thrust multiplier if provided
    /*
     *  Function to check the validity of the input data, and process the maximum thrust multiplier if provided.
     *  \param thrustMagnitudeInterpolator Interpolator for the (maximum) thrust provided to the constructor
     *  \param specificImpulseInterpolator Interpolator for the specific impulse provided to the constructor
     */
    void parseInputDataAndCheckConsistency(
            const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
            const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator );
};

// Function to read a thrust or specific impulse interpolator from a file.
/*
 * Function to read a thrust or specific impulse interpolator from a file.
 * \param coefficientFile Filename containing data to be used as input for interpolator.
 * \return Interpolator set according to data in coefficientFile
 */
std::shared_ptr< interpolators::Interpolator< double, double > > readCoefficientInterpolatorFromFile(
        const std::string coefficientFile );

// Function to create thrust magnitude settings from guidance input and tables
/*
 * Function to create thrust magnitude settings from guidance input and tables
 * \param thrustInputParameterGuidance Object that computes all guidance-input parameters as a function of time/state
 * Note that the number of implemented parameters must be consistent with the numbet of associated entries in
 * thrustIndependentVariables and specificImpulseDependentVariables
 * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
 * thrustIndependentVariables contains an throttle_dependent_thrust entry) as a function of the independent variables.
 * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
 * the 'interpolate' function of thrustMagnitudeInterpolator.
 * \param specificImpulseInterpolator  Interpolator returning the current specific impulse as a function of the
 * independent variables.
 * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
 * input to the 'interpolate' function of specificImpulseInterpolator.
 * \return Thrust magnitude settings for given input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables );

// Function to create thrust magnitude settings from guidance input and tables
/*
 * Function to create thrust magnitude settings from guidance input and tables
 * \param thrustInputParameterGuidance Object that computes all guidance-input parameters as a function of time/state
 * Note that the number of implemented parameters must be consistent with the numbet of associated entries in
 * thrustIndependentVariables and specificImpulseDependentVariables
 * \param thrustMagnitudeDataFile File containing data for the thrust (or maximum thrust if
 * thrustIndependentVariables contains an throttle_dependent_thrust entry) and associated independent variables
 * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
 * the 'interpolate' function of thrustMagnitudeInterpolator.
 * \param specificImpulseDataFile File containing data for the specific impulse and associated independent variables.
 * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
 * input to the 'interpolate' function of specificImpulseInterpolator.
 * \return Thrust magnitude settings for given input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::string specificImpulseDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables );

// Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
/*
 * Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
 * \param thrustInputParameterGuidance Object that computes all guidance-input parameters as a function of time/state
 * Note that the number of implemented parameters must be consistent with the numbet of associated entries in
 * thrustIndependentVariables and specificImpulseDependentVariables
 * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
 * thrustIndependentVariables contains an throttle_dependent_thrust entry) as a function of the independent variables.
 * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
 * the 'interpolate' function of thrustMagnitudeInterpolator.
 * \param constantSpecificImpulse Specific impulse that is to be used at all times.
 * \return Thrust magnitude settings for given input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double constantSpecificImpulse );

// Function to create thrust magnitude settings from guidance input and tables, with constant specific impulse
/*
 * Function to create thrust magnitude settings from guidance input and tables
 * \param thrustInputParameterGuidance Object that computes all guidance-input parameters as a function of time/state
 * Note that the number of implemented parameters must be consistent with the numbet of associated entries in
 * thrustIndependentVariables and specificImpulseDependentVariables
 * \param thrustMagnitudeDataFile File containing data for the thrust (or maximum thrust if
 * thrustIndependentVariables contains an throttle_dependent_thrust entry) and associated independent variables
 * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
 * the 'interpolate' function of thrustMagnitudeInterpolator.
 * \param constantSpecificImpulse Specific impulse that is to be used at all times.
 * \return Thrust magnitude settings for given input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createParameterizedThrustMagnitudeSettings(
        const std::shared_ptr< ThrustInputParameterGuidance > thrustInputParameterGuidance,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double constantSpecificImpulse );

// Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
// maximum allowed axial acceleration (constant specific impulse).
/*
 * Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
 * maximum allowed axial acceleration.
 * \param bodies List of pointers to body objects defining the full simulation environment.
 * \param nameOfBodyWithGuidance Name of body for which thrust guidance is to be created.
 * \param maximumAcceleration Maxmum allowable acceleration due to the thrust force.
 * \param thrustMagnitudeInterpolator Interpolator that computes the maximum thrust as a function of the independent
 * variables.
 * \param thrustIndependentVariables Physical meaning of each of the independent variables used as input to the
 * thrustInterpolator.
 * \param specificImpulse Specific impulse of the propulsion system
 * \param nameOfCentralBody Name of body w.r.t. which thrust guidance is computed (e.g. Earth if the altitude from Earth
 * is used as an independent variable of the thrust).
 * \return Thrust magnitude settings according to input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double specificImpulse,
        const std::string nameOfCentralBody = "" );

// Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
// maximum allowed axial acceleration (constant specific impulse).
/*
 * Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
 * maximum allowed axial acceleration.
 * \param bodies List of pointers to body objects defining the full simulation environment.
 * \param nameOfBodyWithGuidance Name of body for which thrust guidance is to be created.
 * \param maximumAcceleration Maxmum allowable acceleration due to the thrust force.
 * \param thrustMagnitudeDataFile File containing table with maximum thrust values.
 * \param thrustIndependentVariables Physical meaning of each of the independent variables used as input to the
 * thrustInterpolator.
 * \param specificImpulse Specific impulse of the propulsion system
 * \param nameOfCentralBody Name of body w.r.t. which thrust guidance is computed (e.g. Earth if the altitude from Earth
 * is used as an independent variable of the thrust).
 * \return Thrust magnitude settings according to input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const double specificImpulse,
        const std::string nameOfCentralBody = "" );

// Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
// maximum allowed axial acceleration.
/*
 * Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
 * maximum allowed axial acceleration.
 * \param bodies List of pointers to body objects defining the full simulation environment.
 * \param nameOfBodyWithGuidance Name of body for which thrust guidance is to be created.
 * \param maximumAcceleration Maxmum allowable acceleration due to the thrust force.
 * \param thrustMagnitudeInterpolator Interpolator returning the current thrust (or maximum thrust if
 * thrustIndependentVariables contains an throttle_dependent_thrust entry) as a function of the independent variables.
 * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
 * the 'interpolate' function of thrustMagnitudeInterpolator.
 * \param specificImpulseInterpolator  Interpolator returning the current specific impulse as a function of the
 * independent variables.
 * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
 * input to the 'interpolate' function of specificImpulseInterpolator.
 * \param nameOfCentralBody Name of body w.r.t. which thrust guidance is computed (e.g. Earth if the altitude from Earth
 * is used as an independent variable of the thrust).
 * \return Thrust magnitude settings according to input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::shared_ptr< interpolators::Interpolator< double, double > > thrustMagnitudeInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::shared_ptr< interpolators::Interpolator< double, double > > specificImpulseInterpolator,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
        const std::string nameOfCentralBody = "" );

// Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
// maximum allowed axial acceleration.
/*
 * Function to create a thrust magnitude settings based on interpolated maximum thrust, with throttle determined by
 * maximum allowed axial acceleration.
 * \param bodies List of pointers to body objects defining the full simulation environment.
 * \param nameOfBodyWithGuidance Name of body for which thrust guidance is to be created.
 * \param maximumAcceleration Maxmum allowable acceleration due to the thrust force.
 * \param thrustMagnitudeDataFile File containing data for the thrust (or maximum thrust if
 * thrustIndependentVariables contains an throttle_dependent_thrust entry) and associated independent variables
 * \param thrustIndependentVariables List of identifiers for the physical meaning of each of the entries of the input to
 * the 'interpolate' function of thrustMagnitudeInterpolator.
 * \param specificImpulseDataFile File containing data for the specific impulse and associated independent variables.
 * \param specificImpulseDependentVariables List of identifiers for the physical meaning of each of the entries of the
 * input to the 'interpolate' function of specificImpulseInterpolator.
 * \param nameOfCentralBody Name of body w.r.t. which thrust guidance is computed (e.g. Earth if the altitude from Earth
 * is used as an independent variable of the thrust).
 * \return Thrust magnitude settings according to input.
 */
std::shared_ptr< ParameterizedThrustMagnitudeSettings > createAccelerationLimitedParameterizedThrustMagnitudeSettings(
        const SystemOfBodies& bodies,
        const std::string nameOfBodyWithGuidance,
        const double maximumAcceleration,
        const std::string thrustMagnitudeDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > thrustIndependentVariables,
        const std::string specificImpulseDataFile,
        const std::vector< propulsion::ThrustIndependentVariables > specificImpulseDependentVariables,
        const std::string nameOfCentralBody = "" );



} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_THRUSTSETTINGS_H
