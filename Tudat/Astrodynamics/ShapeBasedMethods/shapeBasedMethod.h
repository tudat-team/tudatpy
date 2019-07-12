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

#ifndef SHAPEBASEDMETHOD_H
#define SHAPEBASEDMETHOD_H

#include <Tudat/SimulationSetup/tudatSimulationHeader.h>
#include <math.h>
#include <vector>
#include <Eigen/Dense>
#include <map>

namespace tudat
{
namespace shape_based_methods
{


// Base class for shape based methods.
class ShapeBasedMethod
{
public:

    //! Empty constructor.
    ShapeBasedMethod( ){ };

    //! Default destructor.
    virtual ~ShapeBasedMethod( ) { }

    //! Compute current thrust acceleration in cartesian coordinates.
    virtual Eigen::Vector3d computeCurrentThrustAccelerationVector( const double independentVariable ) = 0;

    //! Compute deltaV.
    virtual double computeDeltaV( ) = 0;

    //! Compute current cartesian state.
    virtual Eigen::Vector6d computeCurrentStateVector( const double independentVariable ) = 0;

    //! Function to compute the shaped trajectory and the propagation fo the full problem.
    virtual void computeShapedTrajectoryAndFullPropagation( simulation_setup::NamedBodyMap& bodyMap,
            basic_astrodynamics::AccelerationMap& accelerationMap,
            const std::string& centralBody,
            const std::string& bodyToPropagate,
            std::function< double ( const double ) > specificImpulseFunction,
            const std::shared_ptr<numerical_integrators::IntegratorSettings<double> > integratorSettings,
            std::pair< std::shared_ptr< propagators::PropagationTerminationSettings >,
                                                    std::shared_ptr< propagators::PropagationTerminationSettings > > terminationSettings,
            std::map< double, Eigen::VectorXd >& fullPropagationResults,
            std::map< double, Eigen::VectorXd >& shapingMethodResults,
            std::map< double, Eigen::VectorXd>& dependentVariablesHistory,
            propagators::TranslationalPropagatorType propagatorType = propagators::cowell,
            const std::shared_ptr< propagators::DependentVariableSaveSettings > dependentVariablesToSave =
            std::shared_ptr< propagators::DependentVariableSaveSettings >( ) ) = 0;



protected:

private:


};


} // namespace shape_based_methods
} // namespace tudat

#endif // SHAPEBASEDMETHOD_H
