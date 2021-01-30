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

#ifndef TUDAT_SIMS_FLANAGAN_H
#define TUDAT_SIMS_FLANAGAN_H

#include <tudat/simulation/simulation.h>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "tudat/astro/low_thrust/simsFlanaganModel.h"
#include "tudat/simulation/optimisation.h"
#include "tudat/astro/low_thrust/shape_based/shapeBasedMethod.h"

namespace tudat
{
namespace low_thrust_trajectories
{


//! Transform thrust model as a function of time into Sims Flanagan thrust model.
std::vector< double > convertToSimsFlanaganThrustModel( std::function< Eigen::Vector3d( const double ) > thrustModelWrtTime,
                                                        const double maximumThrust,
                                                        const double timeOfFlight, const int numberSegmentsForwardPropagation,
                                                        const int numberSegmentsBackwardPropagation );

//! Convert a shape-based thrust profile into a possible initial guess for Sims-Flanagan.
std::function< Eigen::Vector3d( const double ) > getInitialGuessFunctionFromShaping(
        std::shared_ptr< shape_based_methods::ShapeBasedMethod > shapeBasedLeg,
        const int numberSegmentsSimsFlanagan,
        const double timeOfFlight,
        std::function< double( const double ) > specificImpulseFunction,
        std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );



class SimsFlanagan : public low_thrust_trajectories::LowThrustLeg
{
public:

    //! Constructor.
    SimsFlanagan(
            const Eigen::Vector6d& stateAtDeparture,
            const Eigen::Vector6d& stateAtArrival,
            const double centralBodyGravitationalParameter,
            const double initialSpacecraftMass,
            const double maximumThrust,
            const std::function< double ( const double ) > specificImpulseFunction,
            const int numberSegments,
            const double timeOfFlight,
            std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings ) :
        LowThrustLeg( stateAtDeparture, stateAtArrival, timeOfFlight, initialSpacecraftMass, true ),
        centralBodyGravitationalParameter_( centralBodyGravitationalParameter ),
        maximumThrust_( maximumThrust ),
        specificImpulseFunction_( specificImpulseFunction ),
        numberSegments_( numberSegments ),
        optimisationSettings_( optimisationSettings )
    {
        // Calculate number of segments for both the forward propagation (from departure to match point)
        // and the backward propagation (from arrival to match point).
        numberSegmentsForwardPropagation_ = ( numberSegments_ + 1 ) / 2;
        numberSegmentsBackwardPropagation_ = numberSegments_ / 2;

        // Convert the thrust model proposed as initial guess into simplified thrust model adapted to the Sims-Flanagan method.
        if ( optimisationSettings_->initialGuessThrustModel_.first.size( ) != 0 )
        {
            if ( optimisationSettings_->initialGuessThrustModel_.first.size( ) !=
                 static_cast< unsigned int >( 3 * numberSegments ) )
            {
                throw std::runtime_error(
                            "Error when providing an initial guess for Sims-Flanagan, size of the thrust model initial guess unconsistent with number of segments" );
            }
            else
            {
                initialGuessThrustModel_.first = optimisationSettings_->initialGuessThrustModel_.first;
            }

        }
        else
        {
            initialGuessThrustModel_.first = std::vector< double >( );
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

        // Create Sims-Flanagan leg from the best optimisation individual.
        simsFlanaganModel_ = std::make_shared< SimsFlanaganModel >(
                    stateAtDeparture_, stateAtArrival_,
                    centralBodyGravitationalParameter_, initialMass_, maximumThrust_, specificImpulseFunction_,
                    timeOfFlight_, throttles );

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
        return simsFlanaganModel_->getTotalDeltaV( );
    }

    //! Compute current cartesian state.
    Eigen::Vector6d computeCurrentStateVector( const double currentTime );

    //! Compute state history.
    void getTrajectory(
            std::vector< double >& epochsVector,
            std::map< double, Eigen::Vector6d >& propagatedTrajectory )
    {
        simsFlanaganModel_->propagateTrajectory( epochsVector, propagatedTrajectory );
    }

    Eigen::Vector3d computeCurrentThrustForce( double time,
                                          std::function< double ( const double ) > specificImpulseFunction,
                                          std::shared_ptr<numerical_integrators::IntegratorSettings< double > > integratorSettings );

    //! Compute direction thrust acceleration in cartesian coordinates.
    Eigen::Vector3d computeCurrentThrustAccelerationDirection(
            double currentTime );

    //! Compute magnitude thrust acceleration.
    double computeCurrentThrustAccelerationMagnitude(
            double currentTime );

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
    std::shared_ptr< SimsFlanaganModel > getOptimalSimsFlanaganModel( )
    {
        return simsFlanaganModel_;
    }


    //! Retrieve acceleration map (thrust and central gravity accelerations).
    basic_astrodynamics::AccelerationMap retrieveLowThrustAccelerationMap(
            const simulation_setup::SystemOfBodies& bodies,
            const std::string& bodyToPropagate,
            const std::string& centralBody,
            const std::function< double ( const double ) > specificImpulseFunction,
            const std::shared_ptr< numerical_integrators::IntegratorSettings< double > > integratorSettings );


    //! Define appropriate translational state propagator settings for the full propagation.
    std::pair< std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > >,
    std::shared_ptr< propagators::TranslationalStatePropagatorSettings< double > > > createLowThrustTranslationalStatePropagatorSettings(
            const std::string& bodyToPropagate,
            const std::string& centralBody,
            const basic_astrodynamics::AccelerationMap& accelerationModelMap,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave );


protected:

    int convertTimeToLegSegment( double currentTime )
    {
        return simsFlanaganModel_->convertTimeToLegSegment( currentTime );
    }

private:

    double centralBodyGravitationalParameter_;

    //! Maximum allowed thrust.
    double maximumThrust_;

    std::function< double ( const double ) > specificImpulseFunction_;

    //! Number of segments into which the leg is subdivided.
    int numberSegments_;

    //! Initial guess for the optimisation.
    //! The first element contains the thrust throttles corresponding to the initial guess for the thrust model.
    //! The second element defines the bounds around the initial time (in percentage).
    std::pair< std::vector< double >, double > initialGuessThrustModel_;

    //! Optimisation settings.
    std::shared_ptr< simulation_setup::OptimisationSettings > optimisationSettings_;

    //! Fitness vector of the optimisation best individual.
    std::vector< double > championFitness_;

    //! Design variables vector corresponding to the optimisation best individual.
    std::vector< double > championDesignVariables_;

    //! Number of segments for the forward propagation from departure to match point.
    int numberSegmentsForwardPropagation_;

    //! Number of segments for the backward propagation from arrival to match point.
    int numberSegmentsBackwardPropagation_;

    //! Sims-Flanagan leg corresponding to the best optimisation output.
    std::shared_ptr< SimsFlanaganModel > simsFlanaganModel_;

};


} // namespace low_thrust_trajectories
} // namespace tudat

#endif // TUDAT_SIMS_FLANAGAN_H
