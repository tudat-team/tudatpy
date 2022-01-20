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

#ifndef TUDAT_LOW_THRUST_LEG_SETTINGS_H
#define TUDAT_LOW_THRUST_LEG_SETTINGS_H

#include <Eigen/Geometry>

#include <boost/bind/bind.hpp>
#include <functional>

#if( TUDAT_BUILD_WITH_PAGMO )
//    #include "tudat/astro/low_thrust/hybridMethod.h"
    #include "tudat/astro/low_thrust/simsFlanagan.h"
#endif
#include "tudat/astro/low_thrust/shape_based/hodographicShaping.h"
#include "tudat/astro/low_thrust/shape_based/sphericalShaping.h"

#include "tudat/astro/low_thrust/lowThrustLeg.h"

using namespace boost::placeholders;

namespace tudat
{

namespace low_thrust_trajectories
{

//! List of available types of low thrust leg
enum LowThrustLegTypes
{
    hodographic_shaping_leg,
    spherical_shaping_leg
#if( TUDAT_BUILD_WITH_PAGMO )
    ,sims_flanagan_leg,
    hybrid_method_leg
#endif
};


//! Class defining settings for the low-thrust leg.
/*!
 *  Class defining settings for the low-thrust leg. This class is a functional (base) class for
 *  settings of low-thrust leg that require no information in addition to their type.
 *  Classes defining settings for low-thrust leg requiring additional information must be derived from this class.
 */
class LowThrustLegSettings
{
public:

    //! Constructor
    /*!
    * Constructor
    * \param lowThrustLegType Type of low-thrust leg that is to be used.
    */
    LowThrustLegSettings(
            const LowThrustLegTypes lowThrustLegType ):
        lowThrustLegType_( lowThrustLegType ){ }

    //! Destructor.
    virtual ~LowThrustLegSettings( ){ }

    //! Type of low-thrust leg that is to be used.
    LowThrustLegTypes lowThrustLegType_;

};

//! Low-thrust leg settings for hodographic shaping method.
class HodographicShapingLegSettings: public LowThrustLegSettings
{
public:

    //! Constructor
    /*!
    * Constructor
    * \param numberOfRevolutions Number of revolutions of the shape-based trajectory.
    * \param centralBodyGravitationalParameter Gravitational parameter of the central body.
    * \param radialVelocityFunctionComponents Base components of the radial velocity function.
    * \param normalVelocityFunctionComponents Base components of the normal velocity function.
    * \param axialVelocityFunctionComponents Base components of the axial velocity function.
    * \param radialVelocityFunctionComponents Coefficients vector for the components of the radial velocity function.
    * \param radialVelocityFunctionComponents Coefficients vector for the components of the normal velocity function.
    * \param radialVelocityFunctionComponents Coefficients vector for the components of the axial velocity function.
    */
    HodographicShapingLegSettings(
            const int numberOfRevolutions,
            const double centralBodyGravitationalParameter,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& radialVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& normalVelocityFunctionComponents,
            std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > >& axialVelocityFunctionComponents,
            const Eigen::VectorXd freeCoefficientsRadialVelocityFunction,
            const Eigen::VectorXd freeCoefficientsNormalVelocityFunction,
            const Eigen::VectorXd freeCoefficientsAxialVelocityFunction ):
        LowThrustLegSettings( hodographic_shaping_leg ),
        numberOfRevolutions_( numberOfRevolutions ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        radialVelocityFunctionComponents_( radialVelocityFunctionComponents ),
        normalVelocityFunctionComponents_( normalVelocityFunctionComponents ),
        axialVelocityFunctionComponents_( axialVelocityFunctionComponents ),
        freeCoefficientsRadialVelocityFunction_( freeCoefficientsRadialVelocityFunction ),
        freeCoefficientsNormalVelocityFunction_( freeCoefficientsNormalVelocityFunction ),
        freeCoefficientsAxialVelocityFunction_( freeCoefficientsAxialVelocityFunction ){ }

    //! Destructor
    ~HodographicShapingLegSettings( ){ }

    //! Number of revolutions for the shape-based trajectory.
    int numberOfRevolutions_;

    //! Gravitational parameter of the central body.
    bool centralBodyGravitationalParameter_;

    //! Base components for radial velocity function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > radialVelocityFunctionComponents_;

    //! Base components for normal velocity function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > normalVelocityFunctionComponents_;

    //! Base components for axial velocity function.
    std::vector< std::shared_ptr< shape_based_methods::BaseFunctionHodographicShaping > > axialVelocityFunctionComponents_;

    //! Vector containing the coefficients of the radial function.
    Eigen::VectorXd freeCoefficientsRadialVelocityFunction_;

    //! Vector containing the coefficients of the normal function.
    Eigen::VectorXd freeCoefficientsNormalVelocityFunction_;

    //! Vector containing the coefficients of the axial function.
    Eigen::VectorXd freeCoefficientsAxialVelocityFunction_;

};


//! Low-thrust leg settings for spherical shaping method.
class SphericalShapingLegSettings: public LowThrustLegSettings
{
public:

    //! Constructor
    /*!
    * Constructor
    * \param numberOfRevolutions Number of revolutions of the shape-based trajectory.
    * \param centralBodyGravitationalParameter Gravitational parameter of the central body.
    * \param initialValueFreeCoefficient Initial guess for the free coefficient.
    * \param rootFinderSettings Root finder settings to match the required time of flight.
    * \param boundsFreeCoefficient Boundary values for the free coefficient.
    */
    SphericalShapingLegSettings(
            const int numberOfRevolutions,
            const double centralBodyGravitationalParameter,
            const double initialValueFreeCoefficient,
            const std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings,
            const std::pair< double, double > boundsFreeCoefficient = std::make_pair( TUDAT_NAN, TUDAT_NAN ) ):
        LowThrustLegSettings( spherical_shaping_leg ),
        numberOfRevolutions_( numberOfRevolutions ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        initialValueFreeCoefficient_( initialValueFreeCoefficient ),
        rootFinderSettings_( rootFinderSettings ),
        boundsFreeCoefficient_( boundsFreeCoefficient ){ }

    //! Destructor
    ~SphericalShapingLegSettings( ){ }

    //! Number of revolutions for the shape-based trajectory.
    int numberOfRevolutions_;

    //! Gravitational parameter of the central body.
    double centralBodyGravitationalParameter_;

    //! Initial guess for the free coefficient (i.e. coefficient of the second order component of the radial inverse polynomial).
    double initialValueFreeCoefficient_;

    //! Root finder settings, to be used to find the free coefficient value that ensures the time of flight is correct.
    std::shared_ptr< root_finders::RootFinderSettings > rootFinderSettings_;

    //! Bounds for the free coefficient, to be used when trying to match the required time of flight.
    std::pair< double, double > boundsFreeCoefficient_;

};


#if( TUDAT_BUILD_WITH_PAGMO )
//! Low-thrust leg settings for Sims-Flanagan method.
class SimsFlanaganLegSettings: public LowThrustLegSettings
{
public:

    //! Constructor
    /*!
    * Constructor
    * \param maximumThrust Maximum allowed thrust.
    * \param specificImpulseFunction Function returning specific impulse as a function of time.
    * \param numberOfSegments Number of segments into which the trajectory is divided.
    * \param centralBody Central body of the trajectory.
    * \param optimisationAlgorithm Optimisation algorithm to be used.
    */
    SimsFlanaganLegSettings(
            const double maximumThrust,
            const double centralBodyGravitationalParameter,
            const double vehicleInitialMass,
            std::function< double( const double ) > specificImpulseFunction,
            const int numberOfSegments,
            const std::string centralBody,
            std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings ):
        LowThrustLegSettings( sims_flanagan_leg ),
        maximumThrust_( maximumThrust ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        vehicleInitialMass_( vehicleInitialMass ),

        specificImpulseFunction_( specificImpulseFunction ),
        numberSegments_( numberOfSegments ),
        centralBody_( centralBody ),
        optimisationSettings_( optimisationSettings ){ }

    //! Destructor
    ~SimsFlanaganLegSettings( ){ }

    //! Maximum allowed thrust.
    double maximumThrust_;

    double centralBodyGravitationalParameter_;

    double vehicleInitialMass_;

    //! Specific impulse function.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Name of the central body.
    std::string centralBody_;

    //! Optimisation settings.
    std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings_;

};


//! Low-thrust leg settings for hybrid method.
class HybridMethodLegSettings: public LowThrustLegSettings
{
public:

    //! Constructor
    HybridMethodLegSettings(
            const double maximumThrust,
            const double specificImpulse,
            const std::string centralBody,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings,
            std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings,
            const std::pair< double, double > initialAndFinalMEEcostatesBounds ):
        LowThrustLegSettings( hybrid_method_leg ),
        maximumThrust_( maximumThrust ),
        specificImpulse_( specificImpulse ),
        centralBody_( centralBody ),
        integratorSettings_( integratorSettings ),
        optimisationSettings_( optimisationSettings ){ }

    //! Destructor
    ~HybridMethodLegSettings( ){ }

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse.
    double specificImpulse_;

    //! Name of the central body.
    std::string centralBody_;

    //! Integrator settings.
    std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings_;

    //! Optimisation settings.
    std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings_;

    //! Lower and upper bounds for the initial and final MEE costates.
    std::pair< double, double > initialAndFinalMEEcostatesBounds_;

};
#endif


std::shared_ptr< low_thrust_trajectories::LowThrustLeg  > createLowThrustLeg(
        const std::shared_ptr< LowThrustLegSettings >& lowThrustLegSettings,
        const Eigen::Vector6d& stateAtDeparture,
        const Eigen::Vector6d& stateAtArrival,
        const double& timeOfFlight );

} // namespace low_thrust_trajectories

} // namespace tudat

#endif // TUDAT_LOW_THRUST_LEG_SETTINGS_H
