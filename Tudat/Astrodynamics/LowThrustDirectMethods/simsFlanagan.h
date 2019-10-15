/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef SIMSFLANAGAN_H
#define SIMSFLANAGAN_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "Tudat/Astrodynamics/LowThrustDirectMethods/simsFlanaganLeg.h"
#include "Tudat/Astrodynamics/LowThrustDirectMethods/optimisationSettings.h"
#include "Tudat/Astrodynamics/ShapeBasedMethods/shapeBasedMethodLeg.h"

namespace tudat
{
namespace low_thrust_direct_methods
{


//! Transform thrust model as a function of time into Sims Flanagan thrust model.
std::vector< Eigen::Vector3d > convertToSimsFlanaganThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime,
                                                                 const double maximumThrust,
                                                                 const double timeOfFlight, const int numberSegmentsForwardPropagation,
                                                                 const int numberSegmentsBackwardPropagation );

//! Convert a shape-based thrust profile into a possible initial guess for Sims-Flanagan.
std::function< Eigen::Vector3d( const double ) > getInitialGuessFunctionFromShaping(
        std::shared_ptr< shape_based_methods::ShapeBasedMethodLeg > shapeBasedLeg,
        const int numberSegmentsSimsFlanagan,
        const double timeOfFlight,
        std::function< double( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );



class SimsFlanagan : public transfer_trajectories::LowThrustLeg
{
public:

    //! Constructor.
    SimsFlanagan(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double maximumThrust,
            const std::function< double ( const double ) > specificImpulseFunction,
            const int numberSegments,
            const double timeOfFlight,
            simulation_setup::NamedBodyMap bodyMap,
            const std::string bodyToPropagate,
            const std::string centralBody,
            std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings ) :
        LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, bodyMap, bodyToPropagate, centralBody ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ),
        numberSegments_( numberSegments ),
        optimisationSettings_( optimisationSettings )
    {

        // Store initial spacecraft mass.
        initialSpacecraftMass_ = bodyMap_[ bodyToPropagate_ ]->getBodyMass();

        // Calculate number of segments for both the forward propagation (from departure to match point)
        // and the backward propagation (from arrival to match point).
        numberSegmentsForwardPropagation_ = ( numberSegments_ + 1 ) / 2;
        numberSegmentsBackwardPropagation_ = numberSegments_ / 2;

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the Sims-Flanagan method.
        if ( optimisationSettings_->initialGuessThrustModel_.first != nullptr )
        {
            initialGuessThrustModel_.first = convertToSimsFlanaganThrustModel( optimisationSettings_->initialGuessThrustModel_.first ,
                                                                               maximumThrust_, timeOfFlight_, numberSegmentsForwardPropagation_,
                                                                               numberSegmentsBackwardPropagation_ );
        }
        else
        {
            initialGuessThrustModel_.first = std::vector< Eigen::Vector3d >( );
        }
        initialGuessThrustModel_.second = optimisationSettings_->initialGuessThrustModel_.second;

        // Perform optimisation
        std::pair< std::vector< double >, std::vector< double > > bestIndividual = performOptimisation( );
        championFitness_ = bestIndividual.first;
        championDesignVariables_ = bestIndividual.second;

        // Transform best design variables into vector of throttles.
        std::vector< Eigen::Vector3d > throttles;
        for ( int i = 0 ; i < numberSegments ; i++ )
        {
            throttles.push_back( ( Eigen::Vector3d( ) << championDesignVariables_[ i * 3 ], championDesignVariables_[ i * 3 + 1 ],
                    championDesignVariables_[ i * 3 + 2 ] ).finished( ) );
        }

        bodyMap_[ bodyToPropagate_ ]->setConstantBodyMass( initialSpacecraftMass_ );

        // Create Sims-Flanagan leg from the best optimisation individual.
        simsFlanaganLeg_ = std::make_shared< SimsFlanaganLeg >( stateAtDeparture_, stateAtArrival_, maximumThrust_, specificImpulseFunction_,
                                            timeOfFlight_, bodyMap_, throttles, bodyToPropagate_, centralBody_ );

    }

    //! Default destructor.
    ~SimsFlanagan( ) { }

    //! Convert time to independent variable.
    double convertTimeToIndependentVariable( const double time )
    {
        return time;
    }

    //! Convert independent variable to time.
    double convertIndependentVariableToTime( const double independentVariable )
    {
        return independentVariable;
    }

    //! Perform optimisation.
    std::pair< std::vector< double >, std::vector< double > > performOptimisation( );

    //! Compute DeltaV.
    double computeDeltaV( )
    {
        return simsFlanaganLeg_->getTotalDeltaV( );
    }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
    {
        simsFlanaganLeg_->propagateTrajectory( epochsVector, propagatedTrajectory );
    }

    Eigen::Vector3d computeCurrentThrust( double time,
                                          std::function< double ( const double ) > specificImpulseFunction,
                                          std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude(
            double currentTime, std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return thrust acceleration profile.
    void getThrustAccelerationProfile(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::VectorXd >& thrustAccelerationProfile,
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Return best individual.
    std::vector< double > getBestIndividual( )
    {
        return championDesignVariables_;
    }

    //! Return fitness of best individual.
    std::vector< double > getBestIndividualFitness( )
    {
        return championFitness_;
    }


    //! Return best Sims-Flanagan leg after optimisation.
    std::shared_ptr< SimsFlanaganLeg > getOptimalSimsFlanaganLeg( )
    {
        return simsFlanaganLeg_;
    }


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            std::function< double ( const double ) > specificImpulseFunction,
            std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            basic_astrodynamics::AccelerationMap accelerationModelMap,
            std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );



protected:

    int convertTimeToLegSegment( double currentTime )
    {
        return simsFlanaganLeg_->convertTimeToLegSegment( currentTime );
    }

private:

    //! Maximum allowed thrust.
    double maximumThrust_;

    //! Specific impulse function.
    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< std::vector< Eigen::Vector3d >, double > initialGuessThrustModel_;

    //! Optimisation settings.
    std::shared_ptr< transfer_trajectories::OptimisationSettings > optimisationSettings_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Initial mass of the spacecraft.
    double initialSpacecraftMass_;

    //! Number of segments for the forward propagation from departure to match point.
    int numberSegmentsForwardPropagation_;

    //! Number of segments for the backward propagation from arrival to match point.
    int numberSegmentsBackwardPropagation_;

    //! DeltaV corresponding to best trajectory.
    double deltaV_;

    //! Sims-Flanagan leg corresponding to the best optimisation output.
    std::shared_ptr< SimsFlanaganLeg > simsFlanaganLeg_;

};


} // namespace low_thrust_direct_methods
} // namespace tudat

#endif // SIMSFLANAGAN_H
